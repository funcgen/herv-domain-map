#!/usr/bin/env python3
"""Check complete InterProScan JSON against FASTA/XML; optionally restore an omitted root trailer."""
import argparse,hashlib,json,os,tempfile
from pathlib import Path
from xml.etree import ElementTree as E
def fasta(path):
    records={};name=None;parts=[]
    for line in path.open():
        if line.startswith('>'):
            if name is not None:
                if name in records:raise ValueError('Duplicate FASTA identifier: '+name)
                records[name]=''.join(parts)
            name=line[1:].split()[0];parts=[]
        else:parts.append(line.strip())
    if name is not None:
        if name in records:raise ValueError('Duplicate FASTA identifier: '+name)
        records[name]=''.join(parts)
    return records
def validate(path,fasta_path,xml_path=None,repair=False,report=None):
    original=path.read_bytes();text=original.decode('utf-8');candidate=original;repaired=False
    try:document=json.loads(text)
    except json.JSONDecodeError as error:
        if not repair or error.pos!=len(text):raise
        candidate=original+b'\n]}\n';document=json.loads(candidate);repaired=True
    expected=fasta(fasta_path);seen={}
    for record in document['results']:
        seq=record['sequence']
        if hashlib.md5(seq.encode()).hexdigest()!=record['md5']:raise ValueError('JSON sequence/MD5 mismatch')
        for xref in record['xref']:
            key=xref['id']
            if key in seen or expected.get(key)!=seq:raise ValueError('JSON identifier/sequence mismatch: '+key)
            seen[key]=seq
    if seen!=expected:raise ValueError('JSON does not represent the complete FASTA identifier set')
    if xml_path:
        xml_seen={}
        for _,protein in E.iterparse(xml_path,events=('end',)):
            if protein.tag.rsplit('}',1)[-1]!='protein':continue
            seq=next(x.text.strip() for x in protein if x.tag.rsplit('}',1)[-1]=='sequence')
            for x in protein:
                if x.tag.rsplit('}',1)[-1]=='xref':
                    key=x.get('id')
                    if key in xml_seen or expected.get(key)!=seq:raise ValueError('XML identifier/sequence mismatch: '+str(key))
                    xml_seen[key]=seq
            protein.clear()
        if xml_seen!=expected:raise ValueError('XML does not represent the complete FASTA identifier set')
    if repaired:
        fd,name=tempfile.mkstemp(dir=path.parent,prefix=path.name+'.validated-')
        try:
            with os.fdopen(fd,'wb') as output:output.write(candidate)
            os.replace(name,path)
        finally:Path(name).unlink(missing_ok=True)
    result={'JSON_valid':True,'identifiers_verified':len(seen),'unique_sequence_records':len(document['results']),'XML_checked':xml_path is not None,'root_trailer_repaired':repaired,'original_SHA256':hashlib.sha256(original).hexdigest(),'final_SHA256':hashlib.sha256(candidate).hexdigest(),'repair_policy':'only append missing top-level results-array/root-object closing tokens after full record validation'}
    if report:report.write_text(json.dumps(result,indent=2)+'\n')
    return result
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('json',type=Path);p.add_argument('--fasta',required=True,type=Path);p.add_argument('--xml',type=Path);p.add_argument('--repair-omitted-root',action='store_true');p.add_argument('--report',type=Path);a=p.parse_args()
    print(json.dumps(validate(a.json,a.fasta,a.xml,a.repair_omitted_root,a.report),indent=2))
if __name__=='__main__':main()

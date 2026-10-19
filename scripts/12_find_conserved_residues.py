#!/usr/bin/env python3
import argparse,csv
from lxml import etree as E
NS='https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/schemas';N={'i':NS}
FIELDS=['Protein','Domain','Site Description','Conserved Residues']
def parse_interproscan_xml_grouped(xml_path):
 for _,p in E.iterparse(str(xml_path),events=('end',),tag='{'+NS+'}protein'):
  ids=[x.get('id') for x in p.findall('i:xref',N)]
  for m in p.findall('i:matches/i:rpsblast-match',N):
   sig=m.find('i:signature',N);desc=sig.get('desc','unknown') if sig is not None else 'unknown'
   for site in m.findall('.//i:rpsblast-site',N):
    residues=[x.get('residue')+x.get('start') for x in site.findall('i:site-locations/i:site-location',N) if x.get('residue') and x.get('start')]
    if residues:
     for k in ids:yield dict(zip(FIELDS,[k,desc,site.get('description','unknown'),';'.join(residues)]))
  p.clear()
  while p.getprevious() is not None:del p.getparent()[0]
def main():
 p=argparse.ArgumentParser();p.add_argument('-i','--input',required=True);p.add_argument('-o','--output',required=True);a=p.parse_args()
 with open(a.output,'w') as f:
  w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(parse_interproscan_xml_grouped(a.input))
if __name__=='__main__':main()

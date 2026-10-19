#!/usr/bin/env python3
"""Frozen release selection. Native domtbl scores are authoritative; legacy filters are preserved."""
from pathlib import Path
import re,csv,json,sys,collections,argparse
P=Path(__file__).resolve().parents[1];C=['Query','Gene','Domain','DomainType','Score','E-value','Domain_E-value','Coverage','Conservedness','Subfamily','Locus','Name','Line'];EX=['domain_bits','independent_E','hmm_from','hmm_to','ali_from','ali_to','model_length','query_length','Family','new_profile','disposition']
config={}
def locus(q):return re.sub(r'_[0-9]+\|.*$','',q)
def label(c):return 'Full-length' if c>=.95 else 'High' if c>=.8 else 'Moderate' if c>=.6 else 'Low'
def rowkey(r):return (-r['Score'],-r['domain_bits'],r['independent_E'],-r['Coverage'],r['Domain'],r['Query'],r['ali_from'],r['hmm_from'],r['ali_to'],r['hmm_to'])
def filter_file(species,stage,path):
 mapping={a[0]:a[1].lower() for l in CLASSES.read_text().splitlines() if len(a:=l.split())>=2};best={};review={};counts=collections.Counter();byprofile=collections.defaultdict(collections.Counter)
 out=OUT/f'{species}_{stage}.selected.tsv';out.parent.mkdir(exist_ok=True);allfile=None;writer=None
 if stage=='C':
  allfile=(OUT/f'{species}_{stage}.all_passing.tsv').open('w');writer=csv.DictWriter(allfile,fieldnames=C+EX,delimiter='\t',lineterminator='\n');writer.writeheader()
 for line in path.open():
  if line.startswith('#') or not line.strip():continue
  a=line.split();counts['raw_rows']+=1;name='INT_galadriel' if a[0]=='galadriel' else a[0];pr=byprofile[name];pr['reported_rows']+=1
  seqE,seqbits,ce,ie,dbits=map(float,[a[6],a[7],a[11],a[12],a[13]]);new=name in config['model_score_thresholds']
  if seqE>MAXE or ce>MAXD:continue
  counts['pass_original_E']+=1;pr['pass_original_E']+=1
  if new:
   th=config['model_score_thresholds'][name]
   if seqbits<th['sequence_bits'] or dbits<th['domain_bits'] or ie>1e-5:continue
  counts['pass_score_guards']+=1;pr['pass_score_guards']+=1
  q=a[3];loc=locus(q);fam=loc.split('_pos_')[0];cov=round((int(a[16])-int(a[15])+1)/int(a[2]),3);dt=name.split('_')[0];
  if cov<MINC:continue
  disposition='review_outside_scope' if new and fam not in config['new_model_primary_scope'][species] else 'primary';r=dict(zip(C,[q,mapping.get(name,'other'),name,dt,seqbits,seqE,ce,cov,label(cov),re.sub(r'-int.*','',q),loc,loc+'_'+dt,line.strip()]));r.update(dict(zip(EX,[dbits,ie,int(a[15]),int(a[16]),int(a[17]),int(a[18]),int(a[2]),int(a[5]),fam,new,disposition])))
  if writer:writer.writerow(r)
  if disposition=='review_outside_scope':
   key=(loc,name)
   if key not in review or rowkey(r)<rowkey(review[key]):review[key]=r
   counts['outside_scope_passing_rows']+=1;continue
  key=(loc,dt)
  if key not in best or rowkey(r)<rowkey(best[key]):best[key]=r
 if allfile:allfile.close()
 rows=sorted(best.values(),key=lambda r:(r['Locus'],r['DomainType']))
 for suffix,fields in [('selected.tsv',C),('evidence.tsv',C+EX)]:
  with (OUT/f'{species}_{stage}.{suffix}').open('w') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows({k:r[k] for k in fields} for r in rows)
 if stage=='C':
  with (OUT/f'{species}_outside_scope_review.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=C+EX,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(sorted(review.values(),key=lambda r:(r['Locus'],r['Domain'])))
 counts['selected_records']=len(rows);counts['selected_loci']=len(set(r['Locus'] for r in rows));counts['new_profile_selected_records']=sum(r['new_profile'] for r in rows);(OUT/f'{species}_{stage}.filter_counts.json').write_text(json.dumps(dict(counts),indent=2))
 with (OUT/f'{species}_{stage}.profile_filter_counts.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=['Profile','reported_rows','pass_original_E','pass_score_guards'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(dict(Profile=k,**{a:v.get(a,0) for a in ['reported_rows','pass_original_E','pass_score_guards']}) for k,v in sorted(byprofile.items()))
 print(species,stage,json.dumps(dict(counts)),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description='Select one deterministic domain per internal region and domain type from HMMER domtblout')
 ap.add_argument('-i','--input',type=Path,required=True);ap.add_argument('-o','--output',type=Path,required=True)
 ap.add_argument('--domain-classes',type=Path,required=True);ap.add_argument('--profile-thresholds',type=Path)
 ap.add_argument('--primary-scope',type=Path);ap.add_argument('--species',choices=['mouse','human'])
 ap.add_argument('--max-evalue',type=float,default=1e-5);ap.add_argument('--max-domain-evalue',type=float,default=1e-5);ap.add_argument('--min-coverage',type=float,default=0)
 args=ap.parse_args();CLASSES=args.domain_classes;OUT=args.output.parent;OUT.mkdir(parents=True,exist_ok=True);MAXE=args.max_evalue;MAXD=args.max_domain_evalue;MINC=args.min_coverage
 new={'GAG_MERVL','GAG_ERVL_Hs'};present=set()
 with args.input.open() as f:
  for l in f:
   if l.strip() and not l.startswith('#') and l.split()[0] in new:present.add(l.split()[0])
 if present and not (args.profile_thresholds and args.primary_scope and args.species):ap.error('ERV-L models require --profile-thresholds, --primary-scope and --species; E-values alone are insufficient')
 config={'model_score_thresholds':{},'new_model_primary_scope':{args.species:[]}}
 if args.profile_thresholds:
  with args.profile_thresholds.open() as f:
   for r in csv.DictReader(f,delimiter='\t'):
    if float(r['independent_E_max'])!=1e-5:ap.error('This selector pins independent E <=1e-5')
    config['model_score_thresholds'][r['Profile']]={'sequence_bits':float(r['sequence_bits']),'domain_bits':float(r['domain_bits'])}
 if args.primary_scope:config['new_model_primary_scope']=json.loads(args.primary_scope.read_text())
 if not present<=set(config['model_score_thresholds']):ap.error('Missing threshold for a new profile')
 filter_file(args.species or 'custom','C',args.input)
 generated=OUT/f"{args.species or 'custom'}_C.selected.tsv"
 if generated.resolve()!=args.output.resolve():generated.replace(args.output)

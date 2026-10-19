#!/usr/bin/env python3
"""Group domain peptides by exact profile and run native hmmalign."""
import argparse,os,subprocess,tempfile,collections
from pathlib import Path
from Bio import SeqIO
def main():
 p=argparse.ArgumentParser();p.add_argument('fasta');p.add_argument('hmm_db');p.add_argument('output_dir');a=p.parse_args();out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
 groups=collections.defaultdict(list)
 for r in SeqIO.parse(a.fasta,'fasta'):groups[r.id.rsplit('|',3)[2]].append(r)
 bindir=os.environ.get('HMMER_BIN','');tool=lambda n:str(Path(bindir)/n) if bindir else n
 with tempfile.TemporaryDirectory() as tmp:
  for profile,records in sorted(groups.items()):
   if '/' in profile or '..' in profile:raise ValueError('Unsafe profile name')
   fa=Path(tmp)/'group.fa';model=Path(tmp)/'model.hmm';SeqIO.write(records,fa,'fasta')
   with model.open('w') as f:subprocess.run([tool('hmmfetch'),a.hmm_db,profile],stdout=f,check=True)
   subprocess.run([tool('hmmalign'),'-o',str(out/(profile+'.sto')),str(model),str(fa)],check=True)
if __name__=='__main__':main()

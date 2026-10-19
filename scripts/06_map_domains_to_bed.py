#!/usr/bin/env python3
"""Map domains through the frozen zero-based internal BED, not FASTA header starts."""
import csv,argparse
from pathlib import Path
def main():
 p=argparse.ArgumentParser();p.add_argument('-i','--input',required=True);p.add_argument('-o','--output',required=True);p.add_argument('--internal-bed',required=True);a=p.parse_args()
 regions={}
 with open(a.internal_bed) as f:
  for row in csv.reader(f,delimiter='\t'):
   if len(row)>=6:regions[row[3]]=row
 out=[]
 with open(a.input) as f:
  for r in csv.DictReader(f,delimiter='\t'):
   b=regions[r['Locus']];start,end=int(b[1]),int(b[2]);fields=r['Line'].split();u,v=int(fields[17]),int(fields[18]);offset=int(r['Query'].rsplit('|',1)[1].split('-')[0])-1
   left=offset+3*(u-1);right=offset+3*v
   st,en=(start+left,start+right) if b[5]=='+' else (end-right,end-left)
   assert start<=st<en<=end,(r['Query'],st,en)
   name='|'.join([r['Query'],r['Domain'],r['Gene'],r['Conservedness'],r['Coverage']]);out.append([b[0],st,en,name,min(1000,max(0,int(float(r['Score'])))),b[5]])
 with open(a.output,'w') as f:csv.writer(f,delimiter='\t',lineterminator='\n').writerows(sorted(out,key=lambda x:(x[0],x[1],x[2],x[3])))
if __name__=='__main__':main()

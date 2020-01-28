import os, re
import pyximport; pyximport.install()
import fraction_id_calc as fic

afa = './data/KS_precomputed_hmmalign.afa'

d = {}
cur = ''
seq = ''
with open(afa, 'r') as f:
    lines = f.read().splitlines() 
    for ln in lines:
        if re.match("^>", ln) is not None:
            if cur == '':
                cur = ln[1:]
                continue
            else:
                d[cur] = seq
                seq = ''
                cur = ln[1:]
        else:
            seq += ln
d[cur] = seq
            
fidlist = []
for i1, s1 in enumerate(sorted(d)):
    ln = s1+','
    f = []
    for i2, s2 in enumerate(sorted(d)):
        if(i2 > i1):
            f.append(str(fic.ficalc(d[s1], d[s2])))
    ln += ','.join(reversed(f))
    fidlist.append(ln)

outfile = './precomputed_fracid.csv'
with open(outfile, 'w') as o:
    for ln in reversed(fidlist):
        o.write(ln+"\n")

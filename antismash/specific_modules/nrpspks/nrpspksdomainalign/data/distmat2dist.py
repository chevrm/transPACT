import os, sys
import re
import logging
import numpy as np
import subprocess

def reformat_pairwise_dict_distmat(dist_raw, dist_out):
    o = open(dist_out,"w")
    with open(dist_raw) as f:
        lines = f.read().splitlines()
        p = re.compile("^\s*0\.00\s+")
        for ln in reversed(lines):
            if re.match("^\s*0\.00", ln) is not None:
                newln = p.sub('', ln)
                s = re.split("\s+", newln)
                n = s.pop()
                seqname = s.pop()
                i = []
                for score in s:
                    ident = 1 - (float(score)/100)
                    i.append(ident)
                seqname = seqname.replace("----", "|")
                o.write(seqname+',')
                o.write(','.join(str(j) for j in i)+"\n")
    o.close()
    return 1

reformat_pairwise_dict_distmat(sys.argv[1], 'out.dist')

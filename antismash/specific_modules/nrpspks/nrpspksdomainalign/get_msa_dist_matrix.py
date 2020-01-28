import os
import re
import logging
import numpy as np
import subprocess
from Bio import SeqIO
import pyximport; pyximport.install()
import fraction_id_calc as fic

def _map_old_domain_id(domain_names, domain_seqs, indexed_domain_per_cluster):
    domain_dict = dict(zip(domain_names, domain_seqs))
    new_domain_id_1 = [indexed_domain_per_cluster[c].keys() for c in indexed_domain_per_cluster.keys()]
    new_domain_id = [id for sublist in new_domain_id_1 for id in sublist]

    new_domain_dict = {}
    for id in new_domain_id:
        id_info = id.split("|")[0:9]
        id_info[1] = "c"
        id_new1 = "|".join(id_info)
        id_new = ">" + id_new1
        if id_new in domain_dict.keys():
            new_domain_dict[">" + id] = domain_dict[id_new]

    domain_names_new = new_domain_dict.keys()
    domain_seqs_new = new_domain_dict.values()

    return domain_names_new, domain_seqs_new

def _generate_raw_fasta(domain_names_new, domain_seqs_new, data_dir, training_seq_filename, out_dir):
    ori_raw_filename = os.path.join(data_dir, training_seq_filename)
    out_filename = os.path.join(out_dir, 'add_' + training_seq_filename)
    ori_raw = [l for l in open(ori_raw_filename, 'r').read().split('\n') if l != '']
    for n, s in zip(domain_names_new, domain_seqs_new):
        ori_raw.append(n)
        ori_raw.append(s)

    out = open(out_filename, 'w')
    for l in ori_raw:
        out.write('%s\n' % (l))
    out.close()
    return (out_filename)

def _run_mafft(seq_filename):
    logging.info("MAFFT: generating pairwise distance matrix of multiple sequence alignment of all domains")
    outfile = re.sub(r'\.txt', r'.faa', seq_filename)
    outfile_dist = '.'.join([seq_filename, 'hat2'])
    cmd = ' '.join(["mafft", "--retree", "1", "--distout", seq_filename, ">", outfile])
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    err = process.communicate()[1]
    retcode = process.returncode
    if retcode != 0:
        logging.debug('%s\n' % (err))

    return outfile_dist

def _reformat_pairwise_dict_mafft(in_dist_file, out_dist_file):
    seq_nr = int(open(in_dist_file, 'r').read().split('\n')[1])

    s = 3
    e = 3 + seq_nr
    KS_id = open(in_dist_file, 'r').read().split('\n')[s:e]
    KS_id_adj = [id.split('=')[1] for id in KS_id]

    input = [i for i in open(in_dist_file, 'r').read().split() if i != ''][(s + 2 * seq_nr):]
    input_adj = [0.5 * float(i) for i in input]

    len_per_col = list(reversed(range(0, seq_nr)))
    value_matrix = np.zeros((seq_nr, seq_nr), dtype=np.float)
    entry1 = input_adj[0:len_per_col[0]]
    value_matrix[1:, 0] = entry1

    for i in range(1, len(len_per_col)):
        start = sum(len_per_col[:i])
        end = sum(len_per_col[:i + 1])
        entry = input_adj[start:end]
        matrix_row_index = seq_nr - len_per_col[i]
        value_matrix[matrix_row_index:, i] = entry

    out = open(out_dist_file, 'w')
    for i in range(0, len(len_per_col)):
        row_i = value_matrix[i, 0:i]
        row_i_adj = ','.join([str(r) for r in row_i])
        out.write('%s,%s\n' % (KS_id_adj[i], row_i_adj))

    out.close()


def run_get_msa_dist_matrix(domain_names, domain_seqs, data_dir, training_seq_filename, out_dir,
                            indexed_domain_per_cluster):
    new_names, new_seqs = _map_old_domain_id(domain_names, domain_seqs, indexed_domain_per_cluster)
    seq_out = _generate_raw_fasta(new_names, new_seqs, data_dir, training_seq_filename, out_dir)
    dist_raw_filename = _run_mafft(seq_filename=seq_out)
    if os.path.getsize(dist_raw_filename) > 0:
        out_dist_file = re.sub(r'\.txt\.hat2', r'.dist', dist_raw_filename)
        _reformat_pairwise_dict_mafft(in_dist_file=dist_raw_filename, out_dist_file=out_dist_file)
        return out_dist_file
    else:
        logging.debug("Distance matrix is not generated successfully by MAFFT")

def make_new_seq_aln(new_faa, new_alignment, new_names, new_seqs):
    f = open(new_faa,"w")
    for name, seq in zip(new_names, new_seqs):
        name = name.replace("|", "----")
        f.write(name+"\n"+seq+"\n")
    f.close()
    cmd = ' '.join(["mafft", "--retree", "1", "--distout", new_faa, ">", new_alignment])
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    err = process.communicate()[1]
    retcode = process.returncode
    if retcode != 0:
        logging.debug('%s\n' % (err))
    return new_alignment

def dist_from_aln_distmat(alignment, out):
    cmd = ' '.join(["distmat", alignment, "-protmethod", "0", "-outfile", out]) # 0=uncorrected, 1=jukes-cantor
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    err = process.communicate()[1]
    retcode = process.returncode
    if retcode != 0:
        logging.debug('%s\n' % (err))
    return 1

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
                seqname = seqname.replace("----", "|")
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

def precomp_order(dist_file):
    seqorder = []
    with open(dist_file) as f:
        lines = f.read().splitlines()
        for line in lines:
            seqname, rest = line.split(',', 1)
            seqorder.append(seqname)
    return seqorder

def run_individual_distmat(p_order, aln_file, new_names, out_dir, precomp_dist):
    aln_dict = {}
    t = os.path.join(out_dir, 'tmp.afa')
    o = os.path.join(out_dir, 'tmp.dist')
    fa = SeqIO.parse(open(aln_file),'fasta')
    for rec in fa:
        name = rec.id
        name = name.replace('----', '|')
        aln_dict[name] = str(rec.seq)
    dist_new = []
    for query in new_names:
        print "query loop"
        done = 0
        query = query.replace('>','')
        dline = query
        for i in p_order:
            if done % 1000 == 0:
                print "\t"+str(done)
            qnew = query.replace("|", "----")
            inew = i.replace("|", "----")
            tfh = open(t,"w")
            tfh.write('>'+qnew+"\n"+aln_dict[query]+"\n>"+inew+"\n"+aln_dict[i]+"\n")
            tfh.close()
            dist_from_aln_distmat(t, o)
            with open(o) as f:
                lines = f.read().splitlines()
                p = re.compile("^\s*0\.00\s+")
                for ln in reversed(lines):
                    if re.match("^\s*0\.00", ln) is not None:
                        newln = p.sub('', ln)
                        s = re.split("\s+", newln)
                        if len(s) == 3:
                            ident = 1 - (float(s[0])/100)
                            dline = dline+','+str(ident)
                            break
            done+=1
        dist_new.append(dline)
        p_order.append(query)
    dist = os.path.join(out_dir, 'full.dist')
    dout = open(dist,"w")
    with open(precomp_dist) as p:
        lines = p.read().splitlines()
        for ln in lines:
            dout.write(ln+"\n")
    for d in dist_new:
        dout.write(d+"\n")
    dout.close()
    return dist

def run_cython_dist(p_order, aln_file, new_names, out_dir, precomp_dist):
    aln_dict = {}
    t = os.path.join(out_dir, 'tmp.afa')
    o = os.path.join(out_dir, 'tmp.dist')
    fa = SeqIO.parse(open(aln_file),'fasta')
    for rec in fa:
        aln_dict[rec.id] = str(rec.seq)
    dist_new = []
    for query in new_names:
        query = query.replace('>','')
        dline = query
        for i in p_order:
            fid = fic.ficalc(aln_dict[query], aln_dict[i])
            dline = dline+','+str(fid)
        dist_new.append(dline)
        p_order.append(query)
    dist = os.path.join(out_dir, 'full.dist')
    dout = open(dist,"w")
    with open(precomp_dist) as p:
        lines = p.read().splitlines()
        for ln in lines:
            dout.write(ln+"\n")
    for d in dist_new:
        dout.write(d+"\n")
    dout.close()
    return dist

def run_hmmalign(domains, out_dir, hmm):
    tmpfaa = os.path.join(out_dir, 'tmphmm.faa')
    tmpstk =  os.path.join(out_dir, 'tmphmm.stk')
    qaln = os.path.join(out_dir, 'query_aln.afa')
    aln = {}
    for name, seq in domains:
        with open(tmpfaa, "w") as f:
            f.write(name+"\n"+seq+"\n")
        hmmcmd = ' '.join(['hmmalign', '--amino', '--informat', 'FASTA', '-o', tmpstk, hmm, tmpfaa])
        proc = subprocess.Popen(str(hmmcmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        proc.communicate()
        retcode = proc.returncode
        if retcode == 1:
            print "Alignment failed\n"+hmmcmd+"\n"
        with open(tmpstk, 'r') as infile:
            for line in infile:
                if line.startswith("#"):
                    continue
                elif re.match("^[\s\/]", line[0]) is None:
                    header = line.split(' ')[0]
                    algn = line.split(' ')[-1]
                    algn = ''.join([pos for pos in algn if (pos.isupper() or pos == '-')])
                    if header in aln.keys():
                        aln[header] += algn
                    else:
                        aln[header] = algn
    with open(qaln, "w") as out:
        for h in aln:
            out.write('>'+h+"\n"+aln[h]+"\n")
    return qaln

def combine_aln(a1, a2, out):
    with open(out, 'w') as o:
        with open(a1, 'r') as i1:
            for line in i1:
                o.write(line)
        with open(a2, 'r') as i2:
            for line in i2:
                o.write(line)
    
def run_get_msa_dist_matrix_precomputed(domain_names, domain_seqs, data_dir, training_seq_filename, out_dir,indexed_domain_per_cluster, ref_alignment, precomputed_dist):
    new_names, new_seqs = _map_old_domain_id(domain_names, domain_seqs, indexed_domain_per_cluster)    
    hmm = os.path.join(data_dir, 'KS_rawseq_pred_training_transATPKS.hmm')
    newaln = run_hmmalign(zip(new_names, new_seqs), out_dir, hmm)
    combined_alignment = os.path.join(out_dir, 'combined.afa')
    combine_aln(ref_alignment, newaln, combined_alignment)
    p_order = precomp_order(precomputed_dist)
    full_dist = run_cython_dist(p_order, combined_alignment, new_names, out_dir, precomputed_dist)
    return full_dist

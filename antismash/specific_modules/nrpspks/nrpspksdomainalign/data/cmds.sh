#!/bin/sh

## Edit both KS_*afa and KS_*txt

rm KS_rawseq_pred_training_transATPKS.h*
hmmbuild KS_rawseq_pred_training_transATPKS.hmm KS_rawseq_pred_training_transATPKS.afa
hmmpress KS_rawseq_pred_training_transATPKS.hmm
cp KS_rawseq_pred_training_transATPKS.afa ks_training.afa
perl -p -i -e 's/\|/----/g' ks_training.afa
distmat -protmethod 0 -outfile ks_training.distmat ks_training.afa
python distmat2dist.py ks_training.distmat
mv out.dist ks_training.fracid

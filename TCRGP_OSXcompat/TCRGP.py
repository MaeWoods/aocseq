cd Downloads/TCRGP-master\ 3

python
import pickle
import ast
import csv
import numpy as np
from matplotlib import pyplot as plt
plt.style.use('ggplot')

from tcrgp309 import *

tcrs_vdj_all = file2dict('data/vdjdb_conf1.tsv',['Species','Gene','Epitope'],['CDR3','V','J','Reference','Meta'])
epis_vdj_all = list(tcrs_vdj_all['HomoSapiens']['TRB'].keys())
tcrs_vdj = {}
for epi in epis_vdj_all:
    n_s=len(tcrs_vdj_all['HomoSapiens']['TRB'][epi])
    if n_s>=50:
        print(epi+': '+str(n_s)+' samples')
        tcrs_vdj[epi]=tcrs_vdj_all['HomoSapiens']['TRB'][epi]
        for i in range(n_s):
            meta = ast.literal_eval(tcrs_vdj[epi][i][4])
            sub_id = meta['subject.id']
            reference = tcrs_vdj[epi][i][3]
            tcrs_vdj[epi][i][3] = reference+'_'+sub_id
            tcrs_vdj[epi][i] = tcrs_vdj[epi][i][:-1]
            
epis_vdj = list(tcrs_vdj.keys())

tcrs_vdj_conds = file2dict('data/vdjdb_conf1.tsv',['Species','Gene','Epitope'],['Epitope gene','Epitope species'])
print('{:22s} {:13s} {:s}'.format('Epitope','Epitope gene', 'Epitope species'))
for epi in epis_vdj:
    row = tcrs_vdj_conds['HomoSapiens']['TRB'][epi][0]
    print('{:22s} {:13s} {:s}'.format(epi,row[0],row[1]))
    
cdrs = create_cdr_dict(alignment='imgt',species=['human'])


control_file = 'data/human_pairseqs_v1_parsed_seqs_probs_mq20_clones_random_nbrdists.tsv'
store_fields=['va_reps','vb_reps','cdr3a','cdr3b']
controls=[]
with open(control_file, newline='') as tsvfile:
    reader = csv.DictReader(tsvfile,delimiter='\t')
    for row in reader:
        entry = [row[s] for s in store_fields]
        cA = cdrs['human']['A'][entry[0].split(';')[0]]
        cB = cdrs['human']['B'][entry[1].split(';')[0]]
        if '*' not in ''.join(cB)+''.join(cA)+entry[2]+entry[3]:
            controls.append(entry)
            
n_controls = len(controls)

subsmat = subsmatFromAA2('HENS920102')
pc_blo = get_pcs(subsmat,d=21)

_,_,_,_ = loso('training_data/examples/vdj_human_ATDALMTGY.csv','human','ATDALMTGY',pc_blo,cdr_types=[[],['cdr3']],m_iters=500,lr=0.005,nZ=0,mbs=0,va='va',vb='vb',cdr3a=None,cdr3b='cdr3b',epis='epitope',subs='subject')


_,_,_,_ = Ploso('/Users/maewoodsphd/TakaLab/Rscripts/TCRsMatched/vdj_human_10.csv','human','YVLDHLIVV',pc_blo,cdr_types=[[],['cdr3']],m_iters=500,lr=0.005,nZ=0,mbs=0,va='va',vb='vb',cdr3a=None,cdr3b='cdr3b',epis='epitope',subs='subject',myfile='/Users/maewoodsphd/TakaLab/Rscripts/TCRsMatched/Example10.txt')


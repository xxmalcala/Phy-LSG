

import os
from Bio import SeqIO
import numpy as np
import random
from collections import defaultdict
import math

def get_seqs(fasta_list, taxon_list):
    seq_db = {i:[] for i in taxon_list}
    for fas in fasta_list:
        all_seqs = {i.id[:10]:f'{i.seq}' for i in SeqIO.parse(fas,'fasta')}
        for k in seq_db.keys():
            if k in all_seqs.keys():
                seq_db[k].append(all_seqs[k])
            else:
                seq_db[k].append('-'*len(list(all_seqs.values())[-1]))
    bad_taxa = []
    for k, v in seq_db.items():
        missing_genes = [list(set(i)) for i in v].count(['-'])
        if missing_genes > len(fasta_list)*0.5:
            bad_taxa.append(k)
    for taxon in bad_taxa:
        del seq_db[taxon]
    return seq_db


def random_gene_order(seq_db, proportions = [.2, .4, .6], replicates = 5):
    num_genes = [n for n in range(len(list(seq_db.values())[0]))]
    gene_proportion = [math.floor(i) for i in np.array(proportions)*len(num_genes)]
    for n in range(1,replicates+1):
        rep_seqs = {}
        for p in gene_proportion:
            genes_to_concat = np.random.choice(num_genes, p, replace=False)
            for k, v in seq_db.items():
                rep_seqs[k] = list(np.array(v)[genes_to_concat])
            with open(f'Rep_{n}.{math.ceil(p/len(num_genes)*100)}percent.RandomGene.fas','w+') as w:
                for k, v in rep_seqs.items():
                    w.write(f'>{k}\n{"".join(v)}\n')


def random_positions(seq_db, project_name, proportions = [.2, .4, .6], replicates = 5):
    concat_matrix = {k:''.join(v) for k, v in seq_db.items()}
    # with open(f'{project_name}.AllConcat.fas','w+') as w:
        # for k, v in concat_matrix.items():
            # w.write(f'>{k}\n{v}\n')
    total_positions = [n for n in range(len(list(concat_matrix.values())[0]))]
    pos_proportion = [math.floor(i) for i in np.array(proportions)*len(total_positions)]
    for n in range(1,replicates+1):
        rep_seqs = {}
        for p in pos_proportion:
            pos_to_concat = np.random.choice(total_positions, p, replace=False)
            for k, v in concat_matrix.items():
                rep_seqs[k] = ''.join(np.array(list(v))[pos_to_concat])
            with open(f'Rep_{n}.{math.ceil(p/len(total_positions)*100)}percent.RandomPos.fas','w+') as w:
                for k, v in rep_seqs.items():
                    w.write(f'>{k}\n{"".join(v)}\n')

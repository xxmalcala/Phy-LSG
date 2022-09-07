#!/usr/bin/env python3


'''
blah blah
'''

# Author: Xyrus Maurer-Alcala
# Contact:
# Last Modified:
# usage: python

# Dependencies:
# Python3, BioPython

import os
import shutil
import sys

import pandas as pd

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def check_start_stop(seq, aln, q_start, q_end, gcode = 1):
    # Hunts for start and valid stop codons in the transcript by extending the
    # search beyond the reported BLAST alignment (but just by 10% of the
    # alignment length).
    max_dist = round(aln*0.1)-round(aln*0.1)%3
    max_start = max(q_start-max_dist, q_start%3)
    max_end = min(q_end+max_dist, len(seq)-(len(seq)-q_end)%3)

    # Looks for start codons in the expanded 5' end of the transcript.
    fprime_utr = [f'{Seq(seq[n:n+3]).translate(gcode)}' for n in range(max_start-1, q_start, 3)]
    stop_cdns = [pos for pos, char in enumerate(fprime_utr) if char == '*']
    start_cdns = [pos for pos, char in enumerate(fprime_utr) if char == 'M']
    if start_cdns:
        if stop_cdns:
            try:
                new_start = min(filter(lambda x: x > stop_cdns[0],start_cdns))*3-1
            except:
                new_start = q_start-1
        else:
            new_start = max_start+start_cdns[0]*3-1
    else:
        new_start = q_start-1

    # Looks for valid stop codons in the expanded 3' end of the transcript.
    try:
        new_end = q_end+3+[f'{Seq(seq[n:n+3]).translate(gcode)}' for n in range(q_end, max_end, 3)].index('*')*3
    except ValueError:
        new_end = q_end - (q_end-new_start)%3

    # The "best" ORF that can be identified given the BLASTX report and evidence.
    fin_orf = seq[new_start:new_end]

    return [fin_orf, f'{Seq(fin_orf).translate(gcode)}']


def call_orf(work_dir, taxon_code, tsv, fasta_file, gcode = 1):
    orf_fas_ntd = f'{work_dir}OG_Hits/{taxon_code}.OGhits.NTD.fas'
    orf_fas_aa = orf_fas_ntd.replace('NTD.fas','AA.fas')
    orf_rename_tsv = tsv.replace("hits.tsv","hits.NameConversion.tsv")

    final_orfs = {}
    orf_db = pd.read_table(tsv, header=None, index_col=False)
    orfs_to_call = {i.id:[f'{i.seq}'] for i in SeqIO.parse(fasta_file,'fasta')}

    for index, row in orf_db.iterrows():
        # Captures the useful part of the BLAST report: query, alignment length,
        # query's alignment  start position, query's alignment end position.
        q, h, aln, q_start, q_end = row[[0,1,3,6,7]]
        og = f'OG5_{h.split("_")[-1]}'
        orig_seq = orfs_to_call[q][0]
        final_orfs[q] = [og]
        if q_start < q_end:
            # Snippet below is for some quality control (easily validated).
            # initial_orf = orig_seq[q_start-1:q_end]
            final_orfs[q] += check_start_stop(orig_seq, aln, q_start, q_end, gcode)
        else:
            orfs_to_call[q] += [og]
            qs = len(orig_seq)-q_start+1
            qe = len(orig_seq)-q_end+1
            orig_seq_rc = f'{Seq(orig_seq).reverse_complement()}'
            # initial_orf = orig_seq_rc[qs-1:qe]
            final_orfs[q] += check_start_stop(orig_seq_rc, aln, qs, qe, gcode)

    for k, v in final_orfs.items():
        if 'Cov' in k:
            new_name = f'{k.split("_Len_")[0]}_Len_{len(v[-2])}_Cov_' \
                f'{k.split("_Cov_")[1]}_{v[0]}'
        else:
            new_name = f'{k.split("_Len_")[0]}_Len_{len(v[-2])}_{v[0]}'

        final_orfs[k] += [
            new_name,
            f'>{new_name}\n{v[-2]}\n',
            f'>{new_name}\n{v[-1].rstrip("*")}\n']

    with open(orf_fas_ntd,'w+') as ntd:
        for v in final_orfs.values():
            ntd.write(v[-2])

    with open(orf_fas_aa,'w+') as aa:
        for v in final_orfs.values():
            aa.write(v[-1])

    with open(orf_rename_tsv,'w+') as w:
        w.write('Original Name\tFinal Name\n')
        for k,v in final_orfs.items():
            w.write(f'{k}\t{v[-3]}\n')

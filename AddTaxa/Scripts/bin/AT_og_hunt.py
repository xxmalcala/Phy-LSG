#!/usr/bin/env python3

'''
blah blah
'''

# Author: Xyrus Maurer-Alcala
# Contact: maurerax@gmail.com or xyrus.maurer-alcala@izb.unibe.ch
# Last Modified:
# usage: python

# Dependencies:
# Python3, BioPython, diamond

import subprocess
import sys

from Bio import SeqIO


def diamond_og_db(db_dir, work_dir, taxon_code, fasta_file,
        threads, is_genome = True):
    dmnd_db = f'{db_dir}db_OG/OGsout.dmnd'
    og_hits_tsv = f'{work_dir}OG_Hits/SpreadSheets/{taxon_code}.OGhits.tsv'

    if is_genome:
        dmnd_cmd = f'diamond blastp -k 1 --quiet --very-sensitive ' \
            f'--subject-cover 50 -p {threads} -d {dmnd_db} -q {fasta_file} ' \
            f'-f 6 -o {og_hits_tsv}'
    else:
        dmnd_cmd = f'diamond blastx -k 1 --quiet --very-sensitive ' \
            f'--subject-cover 50 -p {threads} -d {dmnd_db} -q {fasta_file} ' \
            f'-f 6 -o {og_hits_tsv}'

    dmnd_call = subprocess.call(dmnd_cmd, shell=True,
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    return og_hits_tsv


# need the translation step PRIOR to this! Also, merge this with the translation
# steps?

# parse the OG hits for the transcriptomes, then get the putative ORFs..., then
# write out.

# skip and just assign and translate for genomic data?

def rename_og_hits(db_dir, work_dir, taxon_code, fasta_file,
        threads, is_genome = True):

    og_hits_tsv = diamond_og_db(db_dir, work_dir, taxon_code, fasta_file,
        threads, is_genome = True)

    og_hit_renaming = {}
    for hit in open(og_hits_tsv).readlines():
        orig_name = hit.split("\t")[0]
        og_num = "_".join(hit.split("\t")[1].split("_")[-2:])
        og_hit_renaming[orig_name] = f'{orig_name}_{og_num}'

    if is_genome:
        og_hits_aa_fasta = f'{work_dir}OG_Hits/{taxon_code}.OGhits.AA.fas'
        og_hits_ntd_fasta = og_hits_aa_fasta.replace(".AA.fas",".NTD.fas")
        ntd_fasta_file = fasta_file.replace("AA.fas","NTD.fas")

        with open(og_hits_aa_fasta, "w+") as w:
            for seq_rec in SeqIO.parse(fasta_file, 'fasta'):
                if seq_rec.description in og_hit_renaming:
                    w.write(f'>{og_hit_renaming[seq_rec.description]}\n' \
                        f'{seq_rec.seq}\n')

        with open(og_hits_ntd_fasta, "w+") as w:
            for seq_rec in SeqIO.parse(ntd_fasta_file, 'fasta'):
                if seq_rec.description in og_hit_renaming:
                    w.write(f'>{og_hit_renaming[seq_rec.description]}\n' \
                        f'{seq_rec.seq}\n')

        return og_hits_tsv, og_hits_ntd_fasta, og_hits_aa_fasta
    else:
        return og_hits_tsv

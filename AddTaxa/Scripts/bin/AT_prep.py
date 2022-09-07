#!/usr/bin/env python3

'''
blah blah
'''

# Author: Xyrus Maurer-Alcala
# Contact: maurerax@gmail.com or xyrus.maurer-alcala@izb.unibe.ch
# Last Modified:
# usage: python

# Dependencies:
# Python3, BioPython

import os
import readline
import shutil
import subprocess
import sys

from Bio import SeqIO
from collections import defaultdict
from pathlib import Path


### Evaluate genetic code!!!!


def check_taxon_code(taxon_code):
    valid_t_code = False
    valid_size = (2,2,4)

    while not valid_t_code:
        eval_code = taxon_code.split('_')
        if len(eval_code) == 3:
            if tuple(len(c) for c in eval_code) == valid_size:
                valid_t_code = True
                return taxon_code
            else:
                pass
        else:
            pass
        print(f'\n[Error] Taxon code: {taxon_code} is invalid.')
        print(f'Here is an example of a valid taxon code: "Op_me_Hsap"')
        taxon_code = input(f'\nEnter a new valid 10 character taxon code:  ')


def prep_folders(work_dir, taxon_code, fasta_file, is_genome = True):
    if is_genome:
        new_work_dir = f'{work_dir}/{taxon_code}_WGS_CDS/'
    else:
        new_work_dir = f'{work_dir}/{taxon_code}_Transcripts/'
        Path(f'{new_work_dir}rDNA/').mkdir(parents=True, exist_ok=True)
        Path(f'{new_work_dir}Euk_Eval/SpreadSheets').mkdir(parents=True, exist_ok=True)

    Path(f'{new_work_dir}Original/').mkdir(parents=True, exist_ok=True)
    Path(f'{new_work_dir}OG_Hits/SpreadSheets').mkdir(parents=True, exist_ok=True)
    Path(f'{new_work_dir}Translations/').mkdir(parents=True, exist_ok=True)

    copy_fas = f'{new_work_dir}Original/{fasta_file.split("/")[-1]}'
    shutil.copyfile(fasta_file, copy_fas)

    return new_work_dir


def get_gene_name(seq_name, source):
    if source != "other" and source != 'jgi':
    # for GenBank/RefSeq sourced CDSs!
        if "locus_tag=" in seq_name and "protein_id=" in seq_name:
            acc = seq_name.split("protein_id=")[1].split("]")[0]
            gname = seq_name.split("locus_tag=")[1].split("]")[0]
        # for Ensembl sourced data
        elif "gene:" in seq_name:
            acc = seq_name.split()[0]
            gname = seq_name.split("gene:")[1].split()[0]
        elif "gene=" in seq_name:
            acc = seq_name.split()[0]
            gname = seq_name.split("gene=")[1].split()[0]
        else:
        # Just to check out annoying annotations, like "pseudo", but is (hopefully) just temporary!
            # print(seq_name)
            return None, None
    elif source == 'jgi':
        acc = seq_name.split()[0].split("|")[2]
        gname = acc
    return gname, acc


# Check to make sure that the CDSs being compared are complete, with both
# valid start and stop codons, and without potential frameshifts.
def check_valid_seq(seq):
    ends = ["TGA","TAG","TAA"]
    if len(seq) % 3 != 0:
        return None
    elif seq[:3] != "ATG":
        return None
    elif seq[-3:] not in ends:
        return None
    else:
        return "Valid"


# Parses the FASTA files and returns just the largest isoform for each gene.
# This is not always the best option, but simplifies phylogenomic analyses.
def capture_isoforms(work_dir, taxon_code, fasta_file, source):
    reduced_fasta = f'{taxon_code}.LargeIso.fas'
    working_fas = f'{work_dir}Original/{reduced_fasta}'
    seq_len_db = defaultdict(int)
    seq_acc_db = defaultdict(str)

    for seq_rec in SeqIO.parse(fasta_file,"fasta"):
        gene_name, acc_name = get_gene_name(seq_rec.description, source)
        if not gene_name:
            continue
        if len(seq_rec.seq) > seq_len_db[gene_name]:
            if check_valid_seq(f'{seq_rec.seq.upper()}'):
                seq_len_db[gene_name] = len(seq_rec.seq)
                seq_acc_db[gene_name] = f'>{taxon_code}_XX_{acc_name}\n{seq_rec.seq.upper()}\n'

    with open(working_fas,"w+") as w:
        w.write("".join(seq_acc_db.values()))

    with open(working_fas.replace(".fas",".SeqCodes.tsv"), "w+") as w:
        w.write("Original_Gene_Name\tNew_Name\n")
        for k, v in seq_acc_db.items():
            new_name = v.split("\n")[0].strip(">")
            w.write(f'{k}\t{new_name}\n')

    shutil.copyfile(working_fas, f'{work_dir}Translations/{taxon_code}.NTD.fas')

    return working_fas

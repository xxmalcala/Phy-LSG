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
import shutil
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq



def valid_genetic_codes(genetic_code, check_valid = False):
    gen_code = {"universal":1, "ciliate":6, "blepharisma":4, "euplotes":10,
        "peritrich":30, "mesodinium":29, "ctg":12, "condylostoma":28,
        "1":1, "4":4, "6":6, "10":10, "12":12, "28":28, "29":29, "30":30}
    if check_valid:
        print("Below are the currently supported genetic codes:")
        print("\n".join(gen_code.keys()))
        sys.exit()
    else:
        pass
    if genetic_code.lower() in gen_code:
        return gen_code[genetic_code.lower()]
    else:
        print("[Error]: Invalid genetic code provided. Currently supported"
            " genetic codes are listed below:")
        print("\n".join(gen_code.keys()))
        sys.exit()

def translate_genome(work_dir, taxon_code, fasta_file, genetic_code):
    translated = []
    in_frame_stop = 0
    aa_fasta = f'{work_dir}Translations/{taxon_code}.AA.fas'

    for seq_rec in SeqIO.parse(fasta_file, "fasta"):
        aa = f'{seq_rec.seq.translate(table=genetic_code)}'.rstrip("*")
        if "*" in aa:
            in_frame_stop += 1
            aa = aa.replace("*","X")
        translated.append(f'>{seq_rec.description}\n{aa}\n')

    if in_frame_stop > 0:
        print(f'Warning: Found {in_frame_stop} sequences with in-frame stop codons.')

    with open(aa_fasta, 'w+') as w:
        w.write(''.join(translated))

    return aa_fasta


def translate_ORF(orf, genetic_code):
    return f'{Seq(orf).translate(table=genetic_code)}'

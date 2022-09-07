#!/usr/bin/env python3

##__Updated__: 2021_04_28
##__Author__: Xyrus Maurer-Alcala
##__email__: maurerax@gamil.com

from datetime import datetime
import logging
import os
import shutil
import subprocess
import time

from pathlib import Path


def trim_guid_aln(guid_aln, post_guid_folder, phylo_folder, col_keep=95):
    trim_col = f'{(100-col_keep)/100}'
    renamed_aln = f'{guid_aln.split(".ForTrim")[0]}.{col_keep}trim'

    trim_al_fas_cmd = f'trimal -in {post_guid_folder}/{guid_aln} -out ' \
        f'{phylo_folder}/FASTA/{renamed_aln}.fas -gt {trim_col} -fasta'

    trim_al_phy_cmd = f'trimal -in {post_guid_folder}/{guid_aln} -out ' \
        f'{phylo_folder}/Phylip/{renamed_aln}.phy -gt {trim_col} -phylip'

    trim_al_nex_cmd = f'trimal -in {post_guid_folder}/{guid_aln} -out ' \
        f'{phylo_folder}/Nexus/{renamed_aln}.nex -gt {trim_col} -nexus'


    trimal_fas_call = subprocess.call(trim_al_fas_cmd, shell=True,
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    trimal_phy_call = subprocess.call(trim_al_phy_cmd, shell=True,
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    trimal_nex_call = subprocess.call(trim_al_nex_cmd, shell=True,
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    shutil.copy(
        f'{phylo_folder}/Phylip/{renamed_aln}.phy',
        f'{phylo_folder}/{renamed_aln}.phy')

def prep_for_phylo(params, col_keep=95):

    aln_count = 0
    for aln_file in os.listdir(params["post_guid_folder"]):
        if "ForTrimAl" in aln_file:
            trim_guid_aln(aln_file,
                params["post_guid_folder"],
                params["phylo_folder"],
                col_keep=95)
            aln_count += 1

    curr_time = datetime.now().strftime("%m/%d/%Y %I:%M %p")
    logging.info(
        f'+---------------------------------+\n'\
        f'|          Trim-Al Steps          |\n'\
        f'+---------------------------------+'
        )

    logging.info(f'   Removed columns with ≥95% gap positions for {aln_count} ' \
        f'alignments, finished: {curr_time}')

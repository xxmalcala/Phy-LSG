#!/usr/bin/env python3

##__Updated__: 2021_05_18
##__Author__: Xyrus Maurer-Alcala
##__email__: maurerax@gamil.com or xyrus.maurer-alcala@izb.unibe.ch

"""This script attempts tor resume a PhyloToL run. This will require the
configuration file. This script accounts for "checkpoints" based on the
prescence/absence of files in the user-defined "RunName" folder (inside the
Results folder).

Note (1): this is a work in progress and is largely a conceptual update to
the current release of PhyloToL (see Ceron-Romero et al. 2019, MBE).

Note (2): More updates are to come. This includes better logging, database
hanlding, and user friendly inputs.

If you use this, please cite Ceron-Romero et al. 2019, MBE"""

from datetime import datetime
import logging, glob, shutil, subprocess, sys, time

from pathlib import Path
from Bio import SeqIO

def og_for_guidance(pre_guid_folder, post_guid_folder):
    pre_guid_fas = []
    post_guid_fas = []

    for fas in glob.glob(f'{pre_guid_folder}/*fasta'):
        pre_guid_fas.append(fas.split("/")[-1].split(".")[0])

    for fas in glob.glob(f'{post_guid_folder}/*Guid.fas'):
        post_guid_fas.append(fas.split("/")[-1].split(".")[0])

    ogs_to_fin = [f'{fas}.PreGuidance.fasta' for fas in pre_guid_fas if fas not in post_guid_fas]

    return list(set(ogs_to_fin))

# in prep and not implemented yet...
def resume_phylobuild():
    # check if build-phylo is on
    # if trees are being made, then just do this step!
    pass

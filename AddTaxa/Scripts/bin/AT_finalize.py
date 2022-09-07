#!/usr/bin/env python3

'''
blah blah
'''

# Author: Xyrus Maurer-Alcala
# Contact: maurerax@gmail.com or xyrus.maurer-alcala@izb.unibe.ch
# Last Modified:
# usage: python

# Dependencies:
# Python3, CD-HIT-EST

import shutil
import subprocess
import sys


def store_final_data(db_dir, work_dir, taxon_code, tsv, ntd, aa):
    r2g_ntd = f'{db_dir.split("/Databases")[0]}/ReadyToGo/NTD/'
    r2g_aa = r2g_ntd.replace("/NTD/","/AA/")
    r2g_tsv = r2g_ntd.replace("/NTD/","/TSV/")
    shutil.copy(ntd, r2g_ntd)
    shutil.copy(aa, r2g_aa)
    shutil.copy(tsv, r2g_tsv)

def final_cluster(tsv, ntd, aa):
    pass

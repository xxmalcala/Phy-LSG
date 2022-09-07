#!/usr/bin/env python3

'''
Copies the files in the "ReadyToGo" folder to their respective places in the
DataFiles folder for PhyloToL.
'''

# Author: Xyrus Maurer-Alcala
# Contact: maurerax@gmail.com or xyrus.maurer-alcala@izb.unibe.ch
# Last Modified:
# usage: python

# Dependencies:
# Python3, shutil

import os
import shutil


def copy_files():

    r2g_dir = f'{os.path.realpath(__file__).split("AddTaxa")[0]}' \
		f'AddTaxa/ReadyToGo/'
    datafiles_dir = f'{os.path.realpath(__file__).split("AddTaxa")[0]}' \
        f'DataFiles/AdditionalTaxa/'
    print(r2g_dir)
    for fas in os.listdir(f'{r2g_dir}NTD'):
        if fas.endswith("fas"):
            shutil.copy(f'{r2g_dir}NTD/{fas}',f'{datafiles_dir}ReadyToGo_NTD/')
    for fas in os.listdir(f'{r2g_dir}AA'):
        if fas.endswith("fas"):
            shutil.copy(f'{r2g_dir}AA/{fas}',f'{datafiles_dir}ReadyToGo_AA/')
    for tsv in os.listdir(f'{r2g_dir}TSV'):
        if tsv.endswith("tsv"):
            shutil.copy(f'{r2g_dir}TSV/{tsv}',f'{datafiles_dir}ReadyToGo_TSV/')

if __name__ == '__main__':
    copy_files()

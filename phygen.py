#!/usr/bin/env python3

##__Updated__: 2021_04_28
##__Author__: Xyrus Maurer-Alcala
##__email__: maurerax@gamil.com; xmaurer-alcala@amnh.org


import argparse, logging, os, sys

from datetime import datetime

from bin import pg_prep as prep
from bin import pg_iter_guid as iter_guid
from bin import pg_trim_guid as trimmer

from argparse import RawTextHelpFormatter,SUPPRESS


# add a "continue" argument!!!
def check_args():
	parser = argparse.ArgumentParser(description =
        f'\nTo UPDATE...\n',
        usage=SUPPRESS, formatter_class=RawTextHelpFormatter)

	required_arg_group = parser.add_argument_group('Required Arguments')
	required_arg_group.add_argument("--config", "-c", action = 'store',
		help = ' Configuration file.\n')

	optional_arg_group = parser.add_argument_group('Optional Arguments')
	optional_arg_group.add_argument("--make_config", "-make_config",
		action = "store_true", help = " Create a new config file.\n\n")
	optional_arg_group.add_argument("--resume", "-r",
		action = "store_true", help = ' Resumes run from last checkpoint.\n')

	# optional_arg_group = parser.add_argument_group('Other Options')
	# optional_arg_group.add_argument('--no-phylo','-np', action='store_false',
	# help = ' Do NOT run the phylogenetic tree building steps\n')

	if len(sys.argv[1:]) == 0:
		print (parser.description)
		sys.exit()
	args = parser.parse_args()

	if args.make_config:
		make_slim_config()

	return args


def make_slim_config():
	config_info = f'//\n\nThis file lists the current supported parameters '\
		f'in PhyloToL.\n\nSomething here...\n\n//\n\n### Data Files/Preparation\n' \
		f'\nDataFilesPATH = [Absolute PATH here]\n\n' \
		f'# Number of threads to use.\nThreads = 10\n\n# Folder name with the ' \
		f'reference OGs.\n# Allows for different sets of OGs (e.g. plastid vs '\
		f' nuclear gene OGs).\nRefOGdb = [Name of OG folder; OGs must start ' \
		f'with "OG"]\n\n# A name for the current run of PhyloToL.\nRunName ' \
		f'= [Unique Identifier]\n\n# List of OG names to be run (ensure that ' \
		f'OG names start with "OG").\nOGList = test.txt\n\n# List of taxa to ' \
		f'include in this run. If you want to include all taxa: TaxaForPhyloToL' \
		f' = all\nTaxaForPhyloToL = [all; or name of text file with taxa (one ' \
		f'per line)]\n\n# "BLAST" e-value cut-off of taxa/OGs to include... ' \
		f'Currently UNSUPPORTED.\n# blastCutOff = 1e-20\n\n### GUIDANCE2 ' \
		f'Options\n# Max number of Guidance2 iterations.\nGuidIter = 5\n\n' \
		f'# Keep intermediate Guidance2 files (these are often large).\n' \
		f'Guid_Intermediate = [y/n]\n\n# Cut-off thresholds. Somewhat ' \
		f'arbitrary, but defaults are optimized for DEEP evolution.\nseqCutoff' \
		f' = 0.3\ncolCutoff = 0.4\nresCutoff = 0.0\n\n# Remove sequences that ' \
		f'are very short after Guidance2.\nMinSeqLen_PostG = 10\n\n### TrimAl' \
		f' Options\n# Maximum proportion of gaps per column (default is 95%).\n' \
		f'GapPerc = 95\n\n### Phylogeny Construction Options\n# Automate the0' \
		f' tree building steps.\nBuildPhylo = [y/n]\n\n# Specify the ' \
		f'phylogenetic tree building program [iqtree2 or raxml-ng].\nPhyloProg' \
		f' = iqtree2\n\n# Number of rapid bootstraps to perform (if undesired,' \
		f' default = n).\nBootstraps = [y/n]\n\n# Specify the model to use. ' \
		f' If "auto", iqtree2 default is "MFP" and raxml-ng default is "LG+G"\n' \
		f'Model = auto'

	with open("new_config.txt", "w+") as w:
		w.write(config_info)

	sys.exit(0)


if __name__ == '__main__':
	args = check_args()

	start_time = datetime.now()
	params, depends = prep.prepare_for_phylotol(args.config, args.resume)

	iter_guid.run_guidance_steps(params, args.resume)

	trimmer.prep_for_phylo(params, col_keep=95)

	if params["BuildPhylo"].lower() != "n":
		from bin import phylotol_phylo_build as tree
		tree.build_tree(params)

	curr_time = datetime.now().strftime("%m/%d/%Y %I:%M %p")
	elapsed = datetime.now() - start_time

	print('\n#----------------------------#')
	print(f'\nPhyloToL Run, {params["RunName"]}, finished: {curr_time}')
	print(f'{int(elapsed.seconds/3600)} hours, {int(elapsed.seconds/60)} '\
		f'minutes, {elapsed.seconds} seconds elapsed.\n')
	# print("It is possible to automate the phylogeny building steps with " \
	# 	"scripts in the 'Util' folder.\n")

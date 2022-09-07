#!/usr/bin/env python3

'''
blah blah
'''

# Author: Xyrus Maurer-Alcala
# Contact: maurerax@gmail.com or xyrus.maurer-alcala@izb.unibe.ch
# Last Modified: 2021-05-05
# usage: python

# Dependencies:
# Python3, BioPython, Barrnap, CD-HIT-EST, Diamond

import argparse
import logging
import os
import sys
import time

from bin import AT_prep
from bin import AT_translate

from datetime import datetime

# from bin import phylotol_prep as prep
from argparse import RawTextHelpFormatter,SUPPRESS

### Add checkpoints for transcriptome!!!

### need to add option to skip the isoform step! -- just check the "source" arg

def check_args():
	parser = argparse.ArgumentParser(description =
        f'\nPrepares/formats genomic or transcriptomic data for PhyloToL\n',
        formatter_class=RawTextHelpFormatter, usage=SUPPRESS) # usage=SUPPRESS,
	parser._action_groups.pop()

	# Genome subparser and options
	subparsers = parser.add_subparsers(dest="command")
	genome = subparsers.add_parser("genome", help="genome",
		formatter_class=RawTextHelpFormatter)
	genome._action_groups.pop()

	g_required = genome.add_argument_group("Required Arguments")
	g_optional = genome.add_argument_group("Optional Arguments")

	g_required.add_argument("--fasta", "-fin", action = "store",
		help = " FASTA file of CDSs.\n", dest="fasta")

	g_required.add_argument("--taxon", "-tc", action = "store",
		help = " 10-character taxonomic code (e.g. Op_me_Hsap).\n\n",
		dest="taxon", type = str)

	g_required.add_argument("--genetic_code", "-gc", nargs='?',
		action = "store", default = "1", dest = "genetic_code",
		help = " Genetic code/translation table (NCBI-Based)." \
		" Default is table 1 (universal code).\n\n")

	g_optional.add_argument("--source", "-s", action = "store", dest='source',
		choices=('other', 'NCBI', 'Ensembl', 'jgi'), default = 'ncbi',
		help = " Source of the data. Note that only the largest\n" \
		 " isoforms will be kept for data from NCBI and Ensembl.\n\n")

	g_optional.add_argument("--threads", "-p", action = "store", dest='threads',
		type = int, default = 4,
		help = " Number of CPU threads to use. Default = 4.\n\n")

	g_optional.add_argument("--eval_gen_code", "-e", action = "store_true",
		help = " Attempts to assess the genetic code, based on the " \
		"presence/absence of stop codon reassignments.\n\n")

	g_optional.add_argument("--list_gen_code", action = "store_true",
		help = " Lists the supported genetic codes then quit.\n")

	# Transcriptome subparser and options
	transcriptome = subparsers.add_parser('transcriptome', help='transcriptome',
		formatter_class=RawTextHelpFormatter)
	transcriptome._action_groups.pop()

	t_required = transcriptome.add_argument_group("Required Arguments")
	t_optional = transcriptome.add_argument_group("Optional Arguments")

	t_required.add_argument("--fasta", "-fin", action = "store",
		help = " FASTA file of CDSs.\n\n", dest="fasta")

	t_required.add_argument("--taxon", "-tc", action = "store",
		help = " 10-character taxonomic code (e.g. Op_me_Hsap).\n\n",
		dest="taxon", type = str)

	t_required.add_argument("--genetic_code", "-gc", action = "store",
		default = 1, dest="genetic_code",
		help = " Genetic code/translation table (NCBI-Based).\n")

	t_optional.add_argument("--threads", "-p", action = "store", dest='threads',
		type = int, default = 4,
		help = " Number of CPU threads to use. Default = 4.\n\n")

	t_optional.add_argument("--min_len", "-l", action = "store", dest='length',
		type = int, default = 200,
		help = " Minimum transcript length. Default is 200bp.\n\n")

	t_optional.add_argument("--eval_gen_code", "-e", action = "store_true",
		help = " Uses alignment based assessment of stop-codon reassignments.\n" \
		" Note: Processing will stop until a genetic code\n is provided, following" \
		" assessment.\n\n")

	t_optional.add_argument("--list_gen_code", action = "store_true",
		help = " Lists the supported genetic codes then quit.\n\n")


	if len(sys.argv) < 3:
		print("[Error] AddingTaxa.py COMMANDS [options].")
		print("For help: 'AddingTaxa.py genome -h'")
		# add custom message here?
		# parser.print_help()
	    # parser.print_help(sys.stderr)
		sys.exit(1)

	return parser.parse_args()


def prep_fasta(args):
	args.tc = AT_prep.check_taxon_code(args.taxon)
	if args.command == "genome":
		args.wd = AT_prep.prep_folders(
					args.wd,
					args.tc,
					args.fasta)

		args.w_fas = AT_prep.capture_isoforms(
						args.wd,
						args.tc,
						args.fasta,
						args.source)
	elif args.command == "transcriptome":
		args.wd = AT_prep.prep_folders(
					args.wd,
					args.tc,
					args.fasta,
					False)


def filter_transcripts(args):
	from bin import AT_filt

	size_filt_fas = AT_filt.size_filt_seqs(args.wd,
						args.tc,
						args.fasta,
						args.length)

	rDNA_filt_fas = AT_filt.capture_rRNA(
						args.wd,
						args.tc,
						size_filt_fas,
						args.threads)

	args.euk_fas, args.prok_fas = AT_filt.prok_euk_comp(
									args.wd,
									args.tc,
									rDNA_filt_fas,
									args.threads)


def translate_seqs(args):
	if args.command == "genome":
		args.aa_fas = AT_translate.translate_genome(
						args.wd,
						args.tc,
						args.w_fas,
						args.gc)
	elif args.command == "transcriptome":
		from bin import AT_call_orf
		AT_call_orf.call_orf(
			args.wd,
			args.tc,
			args.og_tsv,
			args.euk_fas,
			args.gc)


def assign_ogs(args):
	from bin import AT_og_hunt
	if args.command == "genome":
		args.og_tsv, args.og_ntd, args.og_aa = (
			AT_og_hunt.rename_og_hits(
				args.db_dir,
				args.wd,
				args.tc,
				args.aa_fas,
				args.threads))
	elif args.command == "transcriptome":
		args.og_tsv = (
			AT_og_hunt.diamond_og_db(
			args.db_dir,
			args.wd,
			args.tc,
			args.euk_fas,
			args.threads,
			False))

def finalize_adding_taxa(args):
	from bin import AT_finalize
	if args.command == 'genome':
		AT_finalize.store_final_data(
			args.db_dir,
			args.wd,
			args.tc,
			args.og_tsv,
			args.og_ntd,
			args.og_aa)


if __name__ == '__main__':
	args = check_args()

	if args.list_gen_code:
		AT_translate.valid_genetic_codes(args.genetic_code, True)
	else:
		args.gc = AT_translate.valid_genetic_codes(args.genetic_code)

	args.db_dir = f'{os.path.realpath(__file__).split("AddTaxa")[0]}' \
		f'AddTaxa/Databases/'

	if args.command == "genome":
		args.wd = f'{os.path.realpath(__file__).split("AddTaxa")[0]}' \
			f'AddTaxa/Genomes'
		prep_fasta(args)
		translate_seqs(args)
		assign_ogs(args)
		finalize_adding_taxa(args)

	else:
		args.wd = f'{os.path.realpath(__file__).split("AddTaxa")[0]}' \
			f'AddTaxa/Transcriptomes'
		prep_fasta(args)
		filter_transcripts(args)
		assign_ogs(args)
		translate_seqs(args)
	# overall_time_start = time.time()

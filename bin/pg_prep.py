#!/usr/bin/env python3

##__Updated__: 2021_04_28
##__Author__: Xyrus Maurer-Alcala'
##__email__: maurerax@gamil.com or xyrus.maurer-alcala@izb.unibe.ch

"""This script prepares the data for the remainder of PhyloToL. This includes
checking the user provided configuration file, as well as the currently
supported set of dependencies.

Note (1): this is a work in progress and is largely a conceptual update to
the current release of PhyloToL (see Ceron-Romero et al. 2019, MBE).

Note (2): More updates are to come. This includes making the "pre-guidance"
files, as well as handling variable databases, taxa, and taxa from outside
the user-defined reference database.

If you use this, please cite Ceron-Romero et al. 2019, MBE"""

import glob
import logging
import os
import shutil
import sys

from datetime import datetime
from pathlib import Path
from Bio import SeqIO

# Checks that the required/optional dependencies are available and in
# the environment PATH if needed.
def check_dependencies():
    d = {"trimAl": shutil.which("trimal"),
        "MAFFT":shutil.which("mafft"),
        "GUIDANCE2":shutil.which("guidance.pl"),
        "RAxML-NG":shutil.which("raxml-ng"),
        "IQTree2":shutil.which("iqtree2")}
    return d

# Parses and stores information from the parameter file.  These will be
# passed along the relevant functions later.
def parse_parameters(config_file):
    param = {}

    for i in open(config_file).readlines():
        if "=" in i and "#" not in i:
            param[i.split()[0]] = i.split()[-1].rstrip("/")

    param["run_folder"] = f'{param["DataFilesPATH"].rsplit("DataFiles")[0]}' \
        f'Results/{param["RunName"]}'
    param["RefOGfolder"] = f'{param["DataFilesPATH"]}/{param["RefOGdb"]}'
    param["guid_folder"] = f'{param["run_folder"]}/Guidance_Steps'
    param["pre_guid_folder"] = f'{param["run_folder"]}/PreGuidance'
    param["post_guid_folder"] = f'{param["run_folder"]}/PostGuidance'
    param["phylo_folder"] = f'{param["run_folder"]}/Phylogeny'

    return param

# If a log/run has been stared already on the same day, then update the
# suffix to ensure subsequent runs do not overwrite current logs files.
def check_log_folder_exist(path):
    filename, extension = os.path.splitext(path)
    counter = 2

    while os.path.exists(path):
        path = f'{filename}.Run_{counter}{extension}'
        counter += 1

    return path

# Logs the information from the parameter file and saves a copy of the presence
# or absence of the current dependencies (e.g. third party software).
def param_depend_log(param):
    # Parse the parameter file and check for dependencies
    dependencies = check_dependencies()

    curr_time = datetime.now().strftime("%m/%d/%Y %I:%M %p")

    # Initializes the log with the date and increments the number of the suffix
    # if a run has already been initiated on the same day.
    current_date_log = check_log_folder_exist(
        f'{param["RunName"]}.{datetime.now().strftime("%b_%d")}.log')

    logging.basicConfig(filename=current_date_log,
        format='%(message)s',
        level=logging.DEBUG)

    # Saves the relevant information to the log file.
    logging.info(f'Thanks for using PhyloToL!')
    logging.info(f'\nPhyloToL Run, {param["RunName"]}, started: {curr_time}')
    print(f'\nPhyloToL Run, {param["RunName"]}, started: {curr_time}')
    logging.info(
        f'\n+---------------------------------+\n'\
        f'|     Parameters for PhyloToL     |\n'\
        f'+---------------------------------+'
        )

    for k, v in param.items():
        if "_folder" not in k:
            logging.info(f'   {k}: {v}')

    logging.info(
        f'\n+---------------------------------+\n'\
        f'|      Checking Dependencies      |\n'\
        f'+---------------------------------+'
        )

    for k, v in dependencies.items():
        if v:
            logging.info(f'   {k}: check')
        else:
            logging.warning(f'   {k}: missing')

    if param["BuildPhylo"] != "n":
        if param["PhyloProg"].lower() == "iqtree2":
            presence = ["IQTree2", dependencies["IQTree2"].lower()]
        elif param["PhyloProg"].lower() == "raxml-ng":
            presence = ["RAxML-NG", dependencies["IQTree2"].lower()]
        else:
            curr_time_2 = datetime.now().strftime("%m/%d/%Y %I:%M %p")
            logging.error(f'\n[Error] No valid phylogenetic tree building ' \
                f' option was specified,\ndespite the specifying to build ' \
                f'trees with the "BuildPhylo" option.\n')
            logging.error(f'Aborting PhyloToL run, {param["RunName"]}: '\
                f'{curr_time_2}\n')
            print("\n[Error] Phylogeny construction set to occur, but no "\
                "supported phylogeny building program specified.")
            print(f'Aborting PhyloToL run, {param["RunName"]}: {curr_time_2}\n')
            sys.exit(1)

        if presence[-1] == "missing":
            curr_time_2 = datetime.now().strftime("%m/%d/%Y %I:%M %p")
            logging.error(f'\n[Error] {presence[0]} is missing from the PATH.')
            logging.error(f'Aborting PhyloToL run, {param["RunName"]}: ' \
                f'{curr_time_2}\n')
            print(f'\n[Error] {presence[0]} is missing from the PATH.')
            print(f'Aborting PhyloToL run, {param["RunName"]}: {curr_time_2}\n')
            sys.exit(1)
        else:
            pass

    return param, dependencies

# Prepares the folder(s) for the current PhyloToL run. If the RunName
# in the config file has not been updated between individual runs, then
# this will abort the run.
def prep_folders(out_folder, pre_guid, post_guid, phylo_f, guid_steps):
    # checks to see if a folder with the current RunName exists in the
    # DatFiles folder and aborts if True.
    if os.path.isdir(f'{out_folder}'):
        curr_time = datetime.now().strftime("%m/%d/%Y %I:%M %p")
        print(f'\nThe folder {out_folder.split("/")[-1]} already exists.')

        if input(f"\nWould you like to overwrite these files? [y/n]   ").lower() != "y":

            print("\nPlease change the RunName in the config file and run again.")
            print(f'\nQuitting PhyloToL to ensure that the previous '
                f'{out_folder.split("/")[-1]} run\nis NOT overwritten and lost.\n')

            logging.warning(f'\n{"#"*30}\nERROR: RunName, '
                f'{out_folder.split("/")[-1]} already exists.\n'
                f'\nAborted PhyloToL run at {curr_time}\n{"#"*30}')
            sys.exit(1)

        else:
            logging.warning(f'\n{"#"*30}\nWARNING: RunName, '
                f'{out_folder.split("/")[-1]} already exists.\n\n'
                f'User decision: Overwrite {curr_time}\n{"#"*30}')

            print(f'\nRemoving prior data for PhyloToL run: '
                f'{out_folder.split("/")[-1]}')

            try:
                shutil.rmtree(f'{pre_guid}')
            except OSError as e:
                print(f'[ERROR]: {e.filename} - {e.strerror}')
            try:
                shutil.rmtree(f'{post_guid}')
            except OSError as e:
                print(f'[ERROR]: {e.filename} - {e.strerror}')
            try:
                shutil.rmtree(f'{guid_steps}')
            except OSError as e:
                print(f'[ERROR]: {e.filename} - {e.strerror}')
            try:
                shutil.rmtree(f'{phylo_f}')
            except OSError as e:
                print(f'[ERROR]: {e.filename} - {e.strerror}')

    Path(f'{pre_guid}').mkdir(parents=True, exist_ok=True)
    Path(f'{post_guid}/Removed_Seqs').mkdir(parents=True, exist_ok=True)
    Path(f'{guid_steps}').mkdir(parents=True, exist_ok=True)
    Path(f'{phylo_f}/Phylip').mkdir(parents=True, exist_ok=True)
    Path(f'{phylo_f}/Nexus').mkdir(parents=True, exist_ok=True)
    Path(f'{phylo_f}/FASTA').mkdir(parents=True, exist_ok=True)
    Path(f'{phylo_f}/Trees').mkdir(parents=True, exist_ok=True)

# Grabs all the sequences to run through PhyloToL using the lists of
# orthologous gene families (OGs) and taxa in the configuration file.
def preguidance_files(data_files_folder, pre_guid_folder, og_list,
    ref_og_folder, taxon_list):

    ### set option to grab the NUC files too?!
    taxon_folder = f'{data_files_folder}/AdditionalTaxa/ReadyToGo_AA'
    og_list_file = f'{data_files_folder}/{og_list}'

    taxa_to_run = []
    ogs_for_run = {}

    og_fasta_files = glob.glob(os.path.join(ref_og_folder, "*.fa*"))
    taxon_fasta_files = glob.glob(os.path.join(taxon_folder, "*.fa*"))

    # Notes the valid taxa.
    if taxon_list != "all":
        taxon_list_file = f'{data_files_folder}/{taxon_list}'
        taxa_to_run = [txn.rstrip() for txn in open(taxon_list_file).readlines() if "_" in txn]

    # Initially grabs data from the reference OGs.
    if og_list == "all":
        for og_fasta in og_fasta_files:
            og_name = og_fasta.split("/")[-1].split(".fa")[0]
            ogs_for_run[og_name] = []

            for seq_rec in SeqIO.parse(og_fasta,"fasta"):
                if taxa_to_run:
                    if seq_rec.description[:10] in taxa_to_run:
                        ogs_for_run[og_name].append(seq_rec)
                else:
                    ogs_for_run[og_name].append(seq_rec)


    # Grabs data from the approved list of reference OGs.
    else:
        ogs_for_run = {i.rstrip():[] for i in open(og_list_file).readlines() if i.startswith("OG")}
        for og_fasta in og_fasta_files:
            og_name = og_fasta.split("/")[-1].split(".fa")[0]
            if og_name in ogs_for_run.keys():
                for seq_rec in SeqIO.parse(og_fasta,"fasta"):
                    if taxa_to_run:
                        if seq_rec.description[:10] in taxa_to_run:
                            ogs_for_run[og_name].append(seq_rec)
                    else:
                        ogs_for_run[og_name].append(seq_rec)


    # Collects sequences from taxa that were added after the reference OG
    # creation.  These are run through the Adding-Taxa steps of PhyloToL.
    for taxon_fasta in taxon_fasta_files:
        if taxa_to_run:
            taxon_name = taxon_fasta.split("/")[-1][:10]
            if taxon_name in taxa_to_run:
                for seq_rec in SeqIO.parse(taxon_fasta,"fasta"):
                    og_name = "_".join(seq_rec.description.split("_")[-2:])
                    if og_name in ogs_for_run.keys():
                        ogs_for_run[og_name].append(seq_rec)
        else:
            for seq_rec in SeqIO.parse(taxon_fasta,"fasta"):
                og_name = "_".join(seq_rec.description.split("_")[-2:])
                if og_name in ogs_for_run.keys():
                    ogs_for_run[og_name].append(seq_rec)

    # Filters out sequences that are too short/long
    for k, v in ogs_for_run.items():
        ref_ogs = [len(i.seq) for i in v if i.id[6:10].islower()]
        avg = sum(ref_ogs)/len(ref_ogs)
        keep = [i for i in v if 0.5*avg <= len(i.seq) <= 1.5*avg]
        ogs_for_run[k] = keep

    # Creates the Pre-GUIDANCE2 fasta files for each OG. BIG issue if potential
    # in-frame stop codons are translated. Replaces with ambiguous "X" and
    # reports the frequency of occurrence.
    ogs_with_inframe_stop = []
    seqs_with_inframe_stop = 0
    for k, v in ogs_for_run.items():
        in_frame_stops = 0
        pre_guid_fasta = f'{pre_guid_folder}/{k}.PreGuidance.fasta'
        with open(pre_guid_fasta,"w+") as w:
            for seq_rec in v:
                if 'Cov' in seq_rec.id:
                    temp = seq_rec.id.replace('Cov_','TEST').replace('Cov','TEST')
                    c = float(temp.split('TEST')[1].split('_')[0])
                    if c < 10:
                        continue
                if '*' in f'{seq_rec.seq}':
                    in_frame_stops += 1
                    seqs_with_inframe_stop += 1
                    ogs_with_inframe_stop.append(k)
                w.write(f'>{seq_rec.description}\n{seq_rec.seq}\n'.replace('*','X'))

    curr_time = datetime.now().strftime("%m/%d/%Y %I:%M %p")
    logging.info(
        f'\n+---------------------------------+\n'\
        f'|           Pre-Guidance          |\n'\
        f'+---------------------------------+'
        )
    logging.info(f'   Generated {len(ogs_for_run)} Pre-Guidance files,' \
        f' finished: {curr_time}')
    if ogs_with_inframe_stop:
        logging.info(f'      [Warning]: {seqs_with_inframe_stop} protein ' \
            f'sequences, in {len(set(ogs_with_inframe_stop))} orthologous '\
            f' gene families, were found with in-frame stop codons "*".\n' \
            f'      Replaced with ambiguous amino acid character ("X") for now.')


# Runs all the above functions and passes along the parameters/dependencies
# to the rest of PhyloToL's steps.
def prepare_for_phylotol(config_file, resume = False):
# config_file = 'phylotol_parameters.txt'
    param = parse_parameters(config_file)
    depend = param_depend_log(param)

    if not resume:
        prep_folders(
            param["run_folder"],
            param["pre_guid_folder"],
            param["post_guid_folder"],
            param["phylo_folder"],
            param["guid_folder"]
        )

    preguidance_files(
        param["DataFilesPATH"],
        param["pre_guid_folder"],
        param["OGList"],
        param["RefOGfolder"],
        param["TaxaForPhyloToL"]
    )

    return param, depend



    # parse tsv to get seq names? -- use that to get names, then fucking roll!
    # can use pandas I guess, but have to rework the adding taxa steps!

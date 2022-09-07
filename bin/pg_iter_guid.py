#!/usr/bin/env python3

##__Updated__: 2021_04_15
##__Author__: Xyrus Maurer-Alcala
##__email__: maurerax@gamil.com or xyrus.maurer-alcala@izb.unibe.ch

"""This script uses GUIDANCE2 to filter out non-homologous sequences from a
given alignment. The resulting output is free to be "cleaned" of gap-rich
columns in the next step(s), based on user preference.

Note (1): this is a work in progress and is largely a conceptual update to
the current release of PhyloToL (see Ceron-Romero et al. 2019, MBE).

Note (2): More updates are to come. This includes better logging, database
hanlding, and user friendly inputs.

If you use this, please cite Ceron-Romero et al. 2019, MBE"""

from datetime import datetime
import logging
import os
import shutil
import subprocess
import time

from pathlib import Path
from Bio import SeqIO

# Renames the sequences in the fasta file for Guidance2 as it cannot
# correctly handle longer names.
def pre_guidance(fasta_file, pre_guid_folder, guid_folder, guidance_iter = 1):
    # Need to create temporary fasta files to track the Guidance2 Iterations
    guid_fasta_file = f'{fasta_file.split(".fas")[0]}.Guid{guidance_iter}.fas'
    pre_guid_fasta = f'{pre_guid_folder}/{fasta_file}'

    seq_codes = {}
    numbered_seqs = []
    seq_count = 1

    # Rename/Number the sequences in each iteration of Guidance2
    for seq_rec in SeqIO.parse(pre_guid_fasta,'fasta'):
        num_name = f'>{seq_count}{seq_rec.description[2:10]}'
        numbered_seqs.append(f'{num_name}\n{seq_rec.seq}')
        seq_codes[num_name.strip('>')] = seq_rec.description
        seq_count += 1

    # Write out the temporary files for the current iteration of Guidance2.
    # The aim is to protect the original FASTA files.
    with open(f'{guid_folder}/{guid_fasta_file}','w+') as w:
        w.write('\n'.join(numbered_seqs))
    return f'{guid_folder}/{guid_fasta_file}', seq_codes


# Checks that the number of protein sequences is enough for Guidance2
# to be run. Also, checks for sequences removed by Guidance2, so they
# can be modified accordingly.
def check_fasta_len(fasta_file, tossed_seq = False):
    grep_cmd = f'grep ">" {fasta_file} | wc -l'
    grep_call = subprocess.Popen(grep_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        shell=True)

    data, error = grep_call.communicate()

    if tossed_seq:
        # print(int(data))
        if int(data) > 0:
            return ('new-iteration', data)
        else:
            return ('finished', '0')
    else:
        if int(data) > 3:
            if int(data) <= 200:
                return ('run','genafpair')
            else:
                return ('run','auto')
        else:
            return ('below',None)


# Performs the Guidance2 iterations, following the parameters set in
# PhyloToL v3. See Cerón-Romero et al. 2019 MBE for more information.
def guidance_cmds(pre_guid_fasta, guid_folder, guidance_iter = 1,
    seq_cutoff = 0.3, col_cutoff = 0.4, threads = 10):

    temp_outdir = f'{guid_folder}/{pre_guid_fasta.split("/")[-1].strip(".fas")}'

    # change to a shutil check?
    # guidance_path = ('/Users/xxma/Desktop/Programs/guidance.v2.02/www/'
    #     f'Guidance/guidance.pl')
    guidance_path = shutil.which('guidance.pl')

    run_guidance, mafftalg = check_fasta_len(pre_guid_fasta)

    if run_guidance == 'below':
        return 'below'

    guid_cmd = (f'perl {guidance_path} ' \
        f'--seqFile {pre_guid_fasta} ' \
        f'--msaProgram MAFFT ' \
        f'--seqType aa ' \
        f'--outDir {temp_outdir} ' \
        f'--seqCutoff {seq_cutoff} ' \
        f'--colCutoff {col_cutoff} ' \
        f'--bootstraps 5 ' \
        f'--proc_num {threads} ' \
        f'--outOrder as_input ' \
        f'--MSA_Param \\\-"-{mafftalg} ' \
        f'--maxiterate 200 '
        f'--thread {threads} '\
        f'--bl 62 --anysymbol"')

    print(f'Running guidance iteration [{guidance_iter}]')

    # print(guid_cmd)

    iter_time = datetime.now().strftime("%m/%d/%Y %I:%M %p")
    logging.info(f'      Guidance iteration [{guidance_iter}] '
        f'started: {iter_time}')

    guid_call = subprocess.call(guid_cmd, shell=True,
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    return temp_outdir


# Copies the "clean" sequences into a new fasta file for the next iteration
# of Guidance2
def prep_next_guidance_iter(post_guid_folder, fas_for_guid, guid_iter, g_status):
    post_guid_fasta = (f'{post_guid_folder}/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names')

    renamed_fasta = (f'{fas_for_guid.rsplit("Guid",1)[0]}Guid'
        f'{guid_iter}.fas')

    shutil.copyfile(post_guid_fasta, renamed_fasta)

    return renamed_fasta


# Ensures no data loss. Saves the removed non-homologous sequences, from the
# current fasta file passed through Guidance, into a new FASTA file.
def collect_removed_seqs(post_guid_folder, removed_seq_fas, seq_codes, og_name):
    seq_folder_out = f'{post_guid_folder}/Removed_Seqs/'

    with open(f'{seq_folder_out}{og_name}.RemovedSeqs.fas','a+') as w:
        for i in SeqIO.parse(f'{removed_seq_fas}.With_Names','fasta'):
            original_name = seq_codes[i.description]
            w.write(f'>{original_name}\n{i.seq}\n')


def parse_columns(msa_file, seq_codes, column_scores, min_col_scr = 0.4, min_aa = 10):
    ### remove columns that are unscored and those below threshold
    good_cols = []
    post_guidance_seqs = []
    renamed_seqs = []
    short_seqs = []

    for i in open(column_scores).readlines()[1:-1]:
        if float(i.split()[-1]) >= min_col_scr:
            good_cols.append(int(i.split()[0])-1)

    for i in SeqIO.parse(msa_file,'fasta'):
        new_name = f'>{seq_codes[i.description]}\n'
        prot_seq = ''.join([aa for pos, aa in enumerate(f'{i.seq}') if pos in good_cols])
        if len(prot_seq) - prot_seq.count("-") > min_aa:
            post_guidance_seqs.append(f'{new_name}{prot_seq}\n')
            renamed_seqs.append(f'{new_name}{i.seq}\n')
        else:
            short_seqs.append(f'{new_name}{i.seq}\n')

    return post_guidance_seqs, renamed_seqs, short_seqs


# Ensure that the original sequence names in the "finsihed" alignment are
# restored. Also removes columns marked for removal by Guidance2.
def finalize_guidance_files(msa_file, temp_guid_folder, og_name, seq_codes,
        post_guid_folder, min_col_scr, min_aa):

    col_scrs = f'{temp_guid_folder}/MSA.MAFFT.Guidance2_res_pair_col.scr'
    # guid_msa = f'{temp_guid_folder}/MSA.MAFFT.aln.With_Names'

    post_guid = f'{post_guid_folder}/{og_name}.PostGuid.fas'
    for_trimming = f'{post_guid.replace(".fas",".ForTrimAl.fas")}'

    post_guid_seqs, renamed_seqs, short_seqs = parse_columns(
                                                msa_file,
                                                seq_codes,
                                                col_scrs,
                                                min_col_scr,
                                                min_aa)
    if len(renamed_seqs) < 4:
        print("[ERROR] Fewer than 4 sequences remain. No Post-Guidance files" \
            " will be created.\n")
    else:
        if short_seqs:
            removed_seqs = f'{post_guid_folder}/Removed_Seqs/{og_name}.RemovedSeqs.fas'
            with open(removed_seqs, "a+") as w:
                w.write("".join(short_seqs))


        with open(f'{for_trimming}',"w+") as w:
            w.write("".join(post_guid_seqs))

        with open(f'{post_guid}',"w+") as w:
            w.write("".join(renamed_seqs))


# Manages the iterations of Guidance2 as well as evaluates its outputs.
def iter_guidance(fasta_file, pre_guid_folder, guid_folder, post_guid_folder,
        max_iter = 10, seq_cutoff = 0.3, col_cutoff = 0.4, threads = 10, min_aa = 10):

    og_name = f'{fasta_file.split("/")[-1].split(".")[0]}'

    fas_for_guid, seq_codes = pre_guidance(
                                fasta_file,
                                pre_guid_folder,
                                guid_folder)

    print(f'\n#{"--"*4} {og_name} {"--"*4}#')

    for guid_iter in range(1,max_iter+1):
        temp_guid_folder = guidance_cmds(
                            fas_for_guid,
                            guid_folder,
                            guid_iter,
                            seq_cutoff,
                            col_cutoff,
                            threads)

        postguid_seqs_tossed = f'{temp_guid_folder}/Seqs.Orig.fas.FIXED.Removed_Seq'

        clean_by_guidance = f'{temp_guid_folder}/Seqs.Orig.fas.FIXED.' \
            f'Without_low_SP_Seq.With_Names'

        guid_status_toss, tossed = check_fasta_len(postguid_seqs_tossed, True)
        guid_status = check_fasta_len(clean_by_guidance)

        if guid_status[0] == "below":
            print(f'\n[WARNING]: Too few sequences are present/remain for ' \
                f'homology assessment with Guidance.\nYou may wish to reduce'
                f' the stringency for the "seqCutoff" in the parameter file.')

            fin_time = datetime.now().strftime("%m/%d/%Y %I:%M %p")
            logging.info(f'      Too few sequences in current iteration for '
                f'{fasta_file.split("/")[-1]}: {fin_time}\n')

            break

        elif guid_iter >= max_iter:
            msa_file = f'{temp_guid_folder}/MSA.MAFFT.aln.With_Names'


            if int(tossed) > 0:
                fin_msa_file = msa_file.replace('aln.','aln.clean.')
                bad_seqs = [i.id for i in SeqIO.parse(f'{postguid_seqs_tossed}.With_Names','fasta')]
                good_seqs = [f'>{i.id}\n{i.seq}\n' for i in SeqIO.parse(msa_file,'fasta') if i.id not in bad_seqs]
                with open(fin_msa_file,'w+') as w:
                    w.write(''.join(good_seqs))
                logging.info(f'      Guidance iteration [{guid_iter}]:'
                    f' {int(tossed)} bad sequence(s)')
            else:
                fin_msa_file = msa_file

            finalize_guidance_files(
                fin_msa_file,
                temp_guid_folder,
                og_name,
                seq_codes,
                post_guid_folder,
                col_cutoff,
                min_aa)

            fin_time = datetime.now().strftime("%m/%d/%Y %I:%M %p")
            logging.info(f'      Guidance iterations finished: {fin_time}\n')

            break

        elif guid_status_toss == "finished":

            msa_file = f'{temp_guid_folder}/MSA.MAFFT.aln.With_Names'

            finalize_guidance_files(
                msa_file,
                temp_guid_folder,
                og_name,
                seq_codes,
                post_guid_folder,
                col_cutoff,
                min_aa)

            fin_time = datetime.now().strftime("%m/%d/%Y %I:%M %p")
            logging.info(f'      Guidance iterations finished: {fin_time}\n')

            break

        elif guid_status_toss == 'new-iteration':
            logging.info(f'      Guidance iteration [{guid_iter}]:'
                f' {int(tossed)} bad sequence(s)')

            fas_for_guid = (f'{fas_for_guid.rsplit("Guid",1)[0]}Guid'
                f'{guid_iter+1}.fas')

            collect_removed_seqs(post_guid_folder,
                postguid_seqs_tossed,
                seq_codes,
                og_name)

            shutil.copyfile(clean_by_guidance, fas_for_guid)



# Removes the intermediate folders/files made by/for Guidance2.
def remove_intermediates(guid_folder):
    for f in os.listdir(guid_folder):
        if not f.endswith('.fas'):
            try:
                shutil.rmtree(f'{guid_folder}/{f}')
            except OSError as e:
                try:
                    file_to_rm = Path(f'{guid_folder}/{f}')
                    file_to_rm.unlink()
                except:
                    print(f'[ERROR]: {e.filename} - {e.strerror}')
        else:
            try:
                file_to_rm = Path(f'{guid_folder}/{f}')
                file_to_rm.unlink()
            except:
                print('ERROR')

# Manages the entire Guidance2 run/steps.
def run_guidance_steps(params, resume = False):
    pre_guid_files = []

    if resume:
        from bin import phylotol_resume as resumer
        pre_guid_files = resumer.og_for_guidance(
                            params["pre_guid_folder"],
                            params["post_guid_folder"]
                            )
        if pre_guid_files:
            remove_intermediates(params["guid_folder"])
            logging.info(
                f'\n+---------------------------------+\n'\
                f'|     Resuming GUIDANCE2 Steps    |\n'\
                f'+---------------------------------+'
                )
    else:
        for fas in os.listdir(params["pre_guid_folder"]):
            if "PreGuid" in fas:
                pre_guid_files.append(fas)

        logging.info(
            f'\n+---------------------------------+\n'\
            f'|         GUIDANCE2 Steps         |\n'\
            f'+---------------------------------+'
            )

    print(
        f'\n+---------------------------------+\n'\
        f'|         GUIDANCE2 Steps         |\n'\
        f'+---------------------------------+'
        )


    for pre_guid_fas in pre_guid_files:
        logging.info(f'   {pre_guid_fas.split("/")[-1]}:')

        iter_guidance(
            pre_guid_fas,
            params["pre_guid_folder"],
            params["guid_folder"],
            params["post_guid_folder"],
            int(params["GuidIter"]),
            float(params["seqCutoff"]),
            float(params["colCutoff"]),
            int(params["Threads"]),
            int(params["MinSeqLen_PostG"]))

        if params["Guid_Intermediate"] == "n":
            remove_intermediates(params["guid_folder"])

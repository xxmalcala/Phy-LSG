#!/usr/bin/env python3

'''
blah blah
'''

# Author: Xyrus Maurer-Alcala
# Contact: maurerax@gmail.com or xyrus.maurer-alcala@izb.unibe.ch
# Last Modified:
# usage: python

# Dependencies:
# Python3, BioPython, CD-HIT-EST, Barrnap

import subprocess
from pathlib import Path
from Bio import SeqIO


def clust_transcripts(work_dir, taxon_code, fasta_file):
    post_clust_fas = f'{work_dir}Original/{taxon_code}.cdHit.fas'

    cd_hit_cmd = f'cd-hit-est -G 0 -c 0.97 -aS 1.0 -aL .005 -i {fasta_file}' \
        f' -o {post_clust_fas}'

    cd_hit_call = subprocess.call(cd_hit_cmd, shell=True,
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    return post_clust_fas


def size_filt_seqs(work_dir, taxon_code, fasta_file, min_len):
    post_clust_fas = clust_transcripts(work_dir, taxon_code, fasta_file)
    size_filt_fas = f'{work_dir}Original/{taxon_code}.{min_len}bp.fas'

    seq_codes = {}
    renamed_transcripts = []
    transcript = 1
    for seq_rec in SeqIO.parse(post_clust_fas,'fasta'):
        if len(seq_rec.seq) >= int(min_len):
            slen = len(seq_rec.seq)
            if "_cov_" in seq_rec.description:
                scov = f'{float(seq_rec.description.split("_cov_")[1].split("_")[0]):.2f}'
                new_name = f'{taxon_code}_XX_Trans_{transcript}_Len_{slen}_Cov_{scov}'
            else:
                new_name = f'{taxon_code}_XX_Trans_{transcript}_Len_{slen}'
            seq_codes[seq_rec.description] = new_name
            renamed_transcripts.append(f'>{new_name}\n{seq_rec.seq}\n')
            transcript += 1

    with open(size_filt_fas,"w+") as w:
        w.write("".join(renamed_transcripts))

    with open(size_filt_fas.replace(".fas",".SeqCodes.tsv"), "w+") as w:
        w.write("Original_Gene_Name\tNew_Name\n")
        for k, v in seq_codes.items():
            w.write(f'{k}\t{v}\n')
    try:
        Path(post_clust_fas).unlink()
    except:
        print(f'[Error] Could not find {post_clust_fas}')

    return size_filt_fas


def capture_rRNA(work_dir, taxon_code, size_filt_fas, threads):
    rDNA_free_fas = f'{work_dir}rDNA/{taxon_code}.rDNA_Free.fas'
    bact_rRNA = f'{work_dir}rDNA/{taxon_code}.rDNA.Bact.fas'

    barrnap_files = {"bac": bact_rRNA,
        "arc":bact_rRNA.replace("Bact.fas","Arch.fas"),
        "euk":bact_rRNA.replace("Bact.fas","Euk.fas"),
        "mito":bact_rRNA.replace("Bact.fas","Mito.fas")}

    for clade, fasta in barrnap_files.items():
        barrnap_cmd = f'barrnap --quiet ' \
            f'--kingdom {clade} ' \
            f'--lencutoff 0.6 ' \
            f'--threads {threads} ' \
            f'--outseq {fasta} ' \
            f'{size_filt_fas}'

        barrnap_call = subprocess.call(barrnap_cmd, shell=True,
            stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    rDNA_seqs = []
    for fasta in barrnap_files.values():
        for seq_rec in SeqIO.parse(fasta,"fasta"):
            rDNA_seqs.append(seq_rec.description.split(":")[-2])

    with open(rDNA_free_fas, "w+") as w:
        for seq_rec in SeqIO.parse(size_filt_fas,"fasta"):
            if seq_rec.description not in rDNA_seqs:
                w.write(f'>{seq_rec.description}\n{seq_rec.seq}\n')

    return rDNA_free_fas


def diamond_prok_euk(work_dir, taxon_code, rDNA_filt_fas, threads):
    bact_hit_tsv = f'{work_dir}Euk_Eval/SpreadSheets/{taxon_code}.BactHits.tsv'
    euk_hit_tsv = f'{work_dir}Euk_Eval/SpreadSheets/{taxon_code}.EukHits.tsv'

    dmnd_db = f'{work_dir.split("Transcriptomes/")[0]}Databases/db_BvsE/'

    dmnd_bact = f'diamond blastx -k 1 --quiet -d {dmnd_db}Bact.dmnd ' \
        f'-p {threads} -f 6 -o {bact_hit_tsv} -q {rDNA_filt_fas}'

    dmnd_euk = dmnd_bact.replace("Bact","Euk")

    dmnd_bact_call = subprocess.call(dmnd_bact, shell=True,
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    dmnd_euk_call = subprocess.call(dmnd_euk, shell=True,
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    return bact_hit_tsv, euk_hit_tsv


def parse_prok_euk_hits(work_dir, taxon_code, rDNA_filt_fas, threads):
    prok_hit, euk_hit = diamond_prok_euk(work_dir, taxon_code, rDNA_filt_fas, threads)
    seq_hit_db = {}

    for line in open(prok_hit).readlines():
        eval = float(line.split("\t")[-2])
        seq_hit_db[line.split("\t")[0]] = [eval]

    for line in open(euk_hit).readlines():
        if line.split("\t")[0] in seq_hit_db.keys():
            eval = float(line.split("\t")[-2])
            seq_hit_db[line.split("\t")[0]].append(eval)
        else:
            seq_hit_db[line.split("\t")[0]] = ["NA",eval,"Euk"]

    for k, v in seq_hit_db.items():
        if len(v) == 1:
            seq_hit_db[k] += ["NA","Prok"]
        elif len(v) == 2:
            if v[1] == 0:
                seq_hit_db[k].append("Euk")
            elif v[0] == 0:
                seq_hit_db[k].append("Prok")
            else:
                prok_euk_ratio = v[0]/v[1]
                euk_prok_ratio = v[1]/v[0]
                if prok_euk_ratio >= 100:
                    seq_hit_db[k].append("Euk")
                elif euk_prok_ratio >= 1000:
                    seq_hit_db[k].append("Prok")
                else:
                    seq_hit_db[k].append("NA")

    return seq_hit_db


def prok_euk_comp(work_dir, taxon_code, rDNA_filt_fas, threads):
    categorized_hits = parse_prok_euk_hits(work_dir, taxon_code, rDNA_filt_fas, threads)
    euk_hit_fas = f'{work_dir}Euk_Eval/{taxon_code}.EukUndHits.fas'
    prok_hit_fas = f'{work_dir}Euk_Eval/{taxon_code}.ProkHits.fas'

    prok_seqs = []
    euk_seqs = []

    for seq_rec in SeqIO.parse(rDNA_filt_fas,'fasta'):
        if seq_rec.description not in categorized_hits:
            euk_seqs.append(f'>{seq_rec.description}\n{seq_rec.seq}\n')
        elif categorized_hits[seq_rec.description][-1] != "Prok":
            euk_seqs.append(f'>{seq_rec.description}\n{seq_rec.seq}\n')
        else:
            prok_seqs.append(f'>{seq_rec.description}\n{seq_rec.seq}\n')

    with open(euk_hit_fas,"w+") as w:
        w.write("".join(euk_seqs))

    with open(prok_hit_fas,"w+") as w:
        w.write("".join(prok_seqs))

    return euk_hit_fas, prok_hit_fas

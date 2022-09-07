

import os
import sys

from Bio import SeqIO
from collections import defaultdict
from glob import glob

from ete3 import Tree

def parse_tree(msa_file, tree_file):
    try:
        t = Tree(tree_file)
    except:
        try:
            t = Tree(tree_file.replace(".fas",""))
            t.set_outgroup(t.get_midpoint_outgroup())
        except:
            print(tree_file)
            sys.exit()
    out_name = tree_file.split('/')[-1].split('.tre')[0]
    minor_clades = []
    taxa_db = defaultdict(list)
    bad_seqs = [i.id for i in SeqIO.parse(msa_file,'fasta') if f'{i.seq}'.count('-')/len(i.seq) > 0.5]
    for i in t.get_leaves():
        if i.name not in bad_seqs:
            taxa_db[i.name[:10]].append(i)
            minor_clades.append(i.name[:5])
    plogs = []
    for k, v in taxa_db.items():
        if len(v) > 1:
            pass
            for node in v:
                iter_parents(out_name, node.name, node, child = None, mjr_c = [], taxon_list = [], recursion_depth = 0)
        else:
            plogs.append(v[0].name)
    plogs += eval_paralogs(out_name, t, msa_file)
    return plogs


def iter_parents(out_name, seq_name, node, child = None, mjr_c = [], taxon_list = [], recursion_depth = 0):
    # Need a method to ignore "self"-hits
    for taxon_seq in node.get_leaves():
        if taxon_seq.name not in taxon_list:
            taxon_list.append(taxon_seq.name)
            mjr_c.append(taxon_seq.name[:2])
    if mjr_c:
        if len(set(mjr_c)) == 1:
            child = node
            try:
                parent = node.up
            except:
                parent = None
            if parent:
                recursion_depth += 1
                iter_parents(out_name,seq_name, parent, child, mjr_c, taxon_list, recursion_depth)
        else:
            if recursion_depth < 1:
                with open(f'{out_name}_CladeSize.txt','a') as w:
                    w.write(f'{seq_name}\t1\n')
            else:
                with open(f'{out_name}_CladeSize.txt','a') as w:
                    w.write(f'{seq_name}\t{len(child.get_leaves())}\n')


def eval_paralogs(out_name, t, msa_file):
    taxa_seq_db = defaultdict(list)
    plogs = []
    for line in open(f'{out_name}_CladeSize.txt').readlines():
        taxa_seq_db[line.split('\t')[0][:10]].append(line.rstrip().split('\t'))
    for k, v in taxa_seq_db.items():
        v.sort(key=lambda x: -int(x[-1]))
        # print(k)
        # print(v)
        if v[0][1] > v[1][1]:
            plogs.append(v[0][0])
        else:
            seqs = [i[0]for i in v if i[1] == v[0][1]]
            seq_dist = []
            for s in seqs:
                node = t.search_nodes(name=s)[0]
                phylo_dist = node.get_distance(node.up)
                seq_dist.append((node.name, phylo_dist))
            plogs.append(min(seq_dist, key=lambda x: x[1])[0])
    return plogs


def parse_msa(plogs, msa_file):
    seqs_to_keep = [f'>{i.id}\n{i.seq}\n' for i in SeqIO.parse(msa_file,'fasta') if i.description in plogs]
    with open(f'Pruned/PlogClean.{msa_file.split("/")[-1]}','w+') as w:
        w.write(''.join(seqs_to_keep))


def select_paralogs(msa_file, tree_file):
    plogs = parse_tree(msa_file, tree_file)
    parse_msa(plogs, msa_file)


def concatenate_msas(taxon_list, project_name):
    taxa_to_merge = {i.rstrip():[] for i in open(taxon_list).readlines() if i.rstrip() != ''}
    gf_order = {}
    cur_pos = 0
    prune_dir = 'Pruned/'
    for f in os.listdir(prune_dir):
        if f.endswith('fas'):
            gf_name = f'TempOG_{f.split("PlogClean.")[1].split(".PostGuid")[0]}'
            x = {i.id[:10]:f'{i.seq}' for i in SeqIO.parse(f'{prune_dir}{f}','fasta')}
            slen = len(list(x.values())[1])
            for t in taxa_to_merge.keys():
                if t in x.keys():
                    taxa_to_merge[t].append(x[t])
                else:
                    taxa_to_merge[t].append('-'*slen)
            gf_order[gf_name] = [cur_pos+1, cur_pos+slen]
            cur_pos += slen
    with open(f'{project_name}.Concat.fas','w+') as w:
        for k, v in taxa_to_merge.items():
            w.write(f'>{k}\n{"".join(v)}\n')
    with open(f'{project_name}.Partitions.nex','w+') as w:
        w.write('#nexus\nbegin sets;\n')
        for k, v in gf_order.items():
                w.write(f'\tcharset {k}aa = {v[0]}-{v[1]};\n')
        w.write('end;')

def get_files(msa_dir, tree_dir):
    # fix this to be flexible!
    file_pair = [(f'{msa_dir}/{f}',f'{tree_dir}/{f.replace("fas","treefile")}') for f in os.listdir(msa_dir) if f.endswith('.fas')]
    for pair in file_pair:
        select_paralogs(pair[0], pair[1])

#    concatenate_msas(taxon_list, project_name)

if __name__ == '__main__':
    if len(sys.argv[1:]) != 4:
        print('usage: python new_concat.py msa-dir tree-dir taxon-list project-name')
        sys.exit(1)
    else:
        msa_dir = sys.argv[1]
        tree_dir = sys.argv[2]
        if len(sys.argv[1:]) == 4:
            taxon_list = sys.argv[3]
            project = sys.argv[4]
    os.mkdir('Pruned')
    get_files(msa_dir, tree_dir)
    concatenate_msas(taxon_list, project)

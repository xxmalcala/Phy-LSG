

import os
import sys

from collections import defaultdict
from ete3 import Tree


def parse_tree(tree_file):
    t = Tree(tree_file)
    all_taxa = list(set([node.name[:10] for node in t.get_leaves()]))
    mnr_clades = list(set([taxon[:5] for taxon in all_taxa]))
    mjr_clades = list(set([taxon[:2] for taxon in all_taxa if taxon[:2] != 'EE']))
    if len(mnr_clades) < 4:
        print(f'Ignoring {tree_file} as it contains fewer than 4 "minor clades"')
        return None, None
    else:
        return t, mjr_clades

def adjust_tree_root(tree_file):
    t, mjr_clades = parse_tree(tree_file)
    tree = Tree(tree_file)
    if not tree:
        return None, None
    priority = ['BaZa', 'Op','Pl','Am','Ex','Sr']
    clade_sizes = {i:[None,0] for i in priority}
    for node in tree.iter_descendants("postorder"):
        mjr_c = [i.name[:2] for i in node.get_leaves() if i.name[:2] in mjr_clades]
        if len(mjr_c) > 1:
            if mjr_c.count('Ba') + mjr_c.count('Za') >= len(mjr_c)-1 and len(mjr_c) > clade_sizes['BaZa'][1]:
                clade_sizes['BaZa'] = [node, len(mjr_c)]
                break
            else:
                for clade in priority[1:]:
                    if mjr_c.count(clade) >= len(mjr_c)-1 and len(mjr_c) > clade_sizes[clade][1]:
                        clade_sizes[clade] = [node, len(mjr_c)]
    for k, v in clade_sizes.items():
        if v[0]:
            tree.set_outgroup(v[0])
            return tree

def check_same_taxon(node, reps = 0):
    taxon_node_names = list(set([leaf.name[:10] for leaf in node.up.get_leaves()]))
    if len(taxon_node_names) == 1:
        parent = node.up
        return parent, True, reps + 1
    else:
        return node.up, False, reps


def check_sisters(tree_file):
    tree = adjust_tree_root(tree_file)
    gene_fam_name = tree_file.split('/')[-1].split('.')[0]
    avg_blen = sum([i.dist for i in tree.get_leaves()])/len(tree.get_leaves())
    verbose_summary = defaultdict(list)
    for node in tree.get_leaves():
        taxon_eval = [node, True, 0]
        query_taxon = node.name[:10]
        mjr = 'same-major'
        mnr = 'same-minor'
        qblen_type = 'long'
        while taxon_eval[1]:
            taxon_eval = check_same_taxon(taxon_eval[0], taxon_eval[2])
        sister_seqs = [taxon.name for taxon in taxon_eval[0].get_leaves()]
        sister_taxa = list(set([i[:10] for i in sister_seqs if i[:10] != query_taxon]))
        mjr_c = list(set([i[:2] for i in sister_taxa]))
        mnr_c = list(set([i[:5] for i in sister_taxa]))
        if len(set(mjr_c)) > 1:
            mjr = 'non-monophyletic'
        elif len(set(mjr_c)) == 1 and mjr_c[0] != query_taxon[:2]:
            mjr = mjr_c[0]
        if len(set(mnr_c)) > 1:
            mnr = 'non-monophyletic'
        elif len(set(mnr_c)) == 1 and mnr_c[0] != query_taxon[:5]:
            mnr = mnr_c[0]
        if node.dist < avg_blen:
            qblen_type = 'short'
        if len(sister_taxa) < 10:
            verbose_summary[query_taxon].append(
                [gene_fam_name, node.name, mjr, mnr, ','.join(sister_taxa),
                node.dist, avg_blen, qblen_type])
    return verbose_summary


def check_many_trees(tree_folder):
    comp_summary = defaultdict(list)
    for tree_file in os.listdir(tree_folder):
        gene_fam_name = tree_file.split('.')[0]
        detailed_summary = check_sisters(f'{tree_folder}/{tree_file}')
        for k, v in detailed_summary.items():
            comp_summary[k] += v
    return comp_summary

def summarize_results(comp_summary, blen_filt = True, verbose = True):
    priority = ['BaZa', 'Op','Pl','Am','Ex','Sr']
    short_summary = {}
    for k, v in comp_summary.items():
        unique_gfs = []
        total_seqs = len(v)
        category = []
        short_blen = []
        for i in v:
            unique_gfs.append(i[0])
            if (blen_filt and i[3][:2] and i[-1] == 'short') or (not blen_filt):
                category.append(i[3][:2])
                short_blen.append(i[-1])
        short_summary[k] = [category.count('sa'), category.count('Op'),
            category.count('Pl'), category.count('Am'), category.count('Ex'),
            category.count('Sr'), category.count('EE'), category.count('Ba'),
            category.count('Za'), category.count('no'), total_seqs,
            f'{category.count("sa")/total_seqs:.3f}']
    with open('Brief_Summary.Phylogenetic_Sisters.tsv','w+') as w:
        w.write('Taxon\tSame-Minor-Clade\tOp\tPl\tAm\tEx\tSr\tEE\tBa\tZa\t' \
                'Non-Monophyletic\tTotal-Tips\tProp. Same\n')
        for k, v in short_summary.items():
            temp = '\t'.join([f'{i}' for i in v])
            w.write(f'{k}\t{temp}\n')

if __name__ == '__main__':
    if len(sys.argv[1:]) == 1:
        tree_folder = sys.argv[1]
    else:
        print('Usage: python count_sisters.py [FOLDER-WITH-TREES]')
        sys.exit()
    trees_summary = check_many_trees(tree_folder)
    summarize_results(trees_summary, blen_filt = False, verbose = True)

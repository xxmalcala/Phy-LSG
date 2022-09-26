# WARNING

This is part of an on-going project to improve accesibility to phylogenomic approaches for "naïve" users to assess "young" lineage-specific genes. As such, this is an **UNFINISHED**, albeit functional, pipeline.

## Description
**Phy-LSG** is a phylogenomic pipeline based on the existing [PhyloToL](https://github.com/Katzlab/PhyloTOL) ([Cerón-Romero et al. 2019, *MBE*](https://academic.oup.com/mbe/article/36/8/1831/5486329)) approach aimed to tackle "young" lineage-specific genes. 

In addition to generating robust **M**ulti-**S**equence **A**lignments (MSAs), support scripts are available to assess possible contamination (SCRIPT NAME) and ability to update databases with additional taxa (transcriptomes or genomes). 

## Dependencies
+ [Barrnap](https://github.com/weizhongli/cdhit)
+ [CUB](https://github.com/xxmalcala/CUB) (may become integrated...)
+ [CD-HIT](https://github.com/weizhongli/cdhit)
+ [DIAMOND](https://github.com/bbuchfink/diamond)
+ [GUIDANCE2](http://guidance.tau.ac.il/source.php)
  - [BioPerl](https://bioperl.org/)
+ [Python 3.6+](https://www.python.org/downloads/)
  - [ETE3](http://etetoolkit.org/)
  - [BioPython](https://biopython.org/wiki/Download)
+ [TrimAl](https://github.com/inab/trimal)
+ [Vsearch](https://github.com/torognes/vsearch)

## Usage
Generate a new configuration file:
```
$ python phygen.py --make_config
```
Run PhyGen on an existing configuration file:
```
$ python phygen.py --config my-config.txt
```
Resume an existing PhyGen run:
```
$ python phygen.py --resume --config my-config.txt
```

## Outputs include:
+ Folder of "pruned" homologs ready for phylogenetic reconstruction.
+ Post-Guidance MSAs with and without informed column removal.
  - Provides opportunity to produce corresponding nucleotide alignments for downstream applications (e.g. HyPhy).
+ Log file summarizing steps, including the numer of sequences removed.
+ FASTA file for each gene family with the removed putative non-homologous sequences.

## Planned Updates
- [ ] Improve multi-threading
- [ ] Phylogeny Reconstruction support (e.g. IQTree2)
- [ ] rDNA database for unusual (i.e. Foraminifera) rDNA removal
- [ ] Conda packaging!

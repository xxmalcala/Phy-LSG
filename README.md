# WARNING

This is part of an on-going "pet" project to improve accesibility to phylogenomic approaches for "naïve" users. As such, this is an **UNFINISHED**, albeit functional, pipeline. It remains unpublished (in its current state) and it would be best to use the existing [PhyloToL](https://github.com/Katzlab/PhyloTOL) pipeline for the time being!

## Description
**PhyGen** is a phylogenomic pipeline based on the existing [PhyloToL](https://github.com/Katzlab/PhyloTOL) approach, with a suite of quality of life changes to improve modularity, reporting, and data curation. 

In addition to generating robust **M**ulti-**S**equence **A**lignments (MSAs), support scripts are available to assess possible contamination (SCRIPT NAME) and ability to update databases with additional taxa (transcriptomes or genomes). 

## Dependencies
+ [Barrnap](https://github.com/weizhongli/cdhit)
+ [CD-HIT](https://github.com/weizhongli/cdhit)
+ [DIAMOND](https://github.com/bbuchfink/diamond)
+ [GUIDANCE2](http://guidance.tau.ac.il/source.php)
  - [BioPerl](https://bioperl.org/)
+ [Python 3.6+](https://www.python.org/downloads/)
  - [ETE3](http://etetoolkit.org/)
  - [BioPython](https://biopython.org/wiki/Download)
+ [TrimAl](https://github.com/inab/trimal)

## Usage
Generate a new configuration file:
```
$ python phygen.py -something
```
Run PhyGen on an existing configuration file:
```
$ python phygen.py -c my-config.txt
```
Resume an existing PhyGen run:
```
$ python phygen.py -r -c my-config.txt
```

## Planned Updates
- [ ] Improve multi-threading
- [ ] Shareable ortholog database
- [ ] Conda packaging!
- [ ] Include phylogeny building support (e.g. IQTree2)
- [ ] Automate gene-species-tree reconciliation tools (ASTRAL-III or GeneRAX)
- [ ] Composition based contaminant screening

Microbial Resistance, Plasmid and Virulence (RPV) identification 
================================================================

## Table of Contents

  * [Introduction](#introduction)
  * [Description](#description)
  * [Dependencies](#dependencies)
  * [Installation](#installation)
    * [Environment](#environment)
    * [Databases](#databases)
  * [Usage](#usage)
  * [References](#reference)


## **Introduction**
--------------------
Monitoring the evolution of resistance, plasmid and virulence of bacterial strains is essential in epidemiology. This tool allows the identification of resistance, plasmid and virulence genes. Several tools and databases are available. For this, we rely on three public and regularly updated databases. [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/) for resistance data, [PlasmidFinder](https://bitbucket.org/genomicepidemiology/plasmidfinder_db/src/master/) for plasmid data and [VirulenceFactors](http://www.mgc.ac.cn/VFs/main.htm)(VFDB) database for virulence data. 


## **Description**
------------------
The program has been developed for use in a command line under a linux system. \
So some basics are necessary. You will have a simple guide to make the program work without being a linux expert

RVP is a python script that takes as input:
- A database in [fasta](https://en.wikipedia.org/wiki/FASTA_format) format
- One or more nucleotide sequences of the assambled strains (the contigs in [fasta](https://en.wikipedia.org/wiki/FASTA_format) format)

And it returns a table mentioning the genes identified for each strain.

## **Dependencies**
------------

- Linux 
- Python >= 3.7
- Databases : ResFinder_db, PlasmidFinder_db, VFDB
- Install requierements : Python packages
- blat


## **Installation**
-------------------

### **Environment**
  - Install python virtual environnement
  ```
  python3 -m venv /path/to/new/virtual/environment
  ```

  - Install python packages
  ```
  pip install -r requierements.txt
  ```
  - Install RPV
  ```
  git clone 
  ```

### **Databases**
Create a directory for downloading last databases updates.

  - ResFinder database 
  ```
  git clone https://aatchon@bitbucket.org/genomicepidemiology/resfinder_db.git
  cd resfinder_db
  bash RPV/scr/rmv_dup_seq.sh named_temporary_database.db named_finaly_database.db
  ```

  - PlasmiFinder database
  ```
  git clone https://aatchon@bitbucket.org/genomicepidemiology/plasmidfinder_db.git
  cd resfinder_db
  bash RPV/scr/rmv_dup_seq.sh named_temporary_database.db named_finaly_database.db
  ```

  - Virulence Factor database
    - Go to [Download page](http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz)
    - Remove duplicate sequences
    ```
    gunzip VFDB_setA_nt.fas.gz
    awk '/^>/{f=!d[$1];d[$1]=1}f' VFDB_setA_nt.fas > named_finaly_database.db
    ```

## **Usage**
-------------

### `rpv.py`

**Usage**:

```console
$ rpv.py [OPTIONS] COMMAND [ARGS]...
```

**Commands**:

* `plasmid`: Search plasmid genes
* `resistance`: Search resistance genes
* `virulence`: Search virulence genes

### `rpv.py command`

**Usage**:

```console
$ rpv.py command [OPTIONS] DATABASE SEQUENCES...
```

**Arguments**:

* `DATABASE`: Database in fasta format  [required]
* `SEQUENCES...`: Sequence(s) in fasta format.  [required]

**Options**:

* `--output PATH`: Output file path of the result (.tsv).
* `--identity INTEGER RANGE`: Minimum percentage of sequence identity [default: 90]
* `--help`: Show this message and exit.

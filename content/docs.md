# `rpv.py`

**Usage**:

```console
$ rpv.py [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `plasmid`: Search plasmid genes
* `resistance`: Search resistance genes
* `virulence`: Search virulence genes

## `rpv.py plasmid`

Search plasmid genes

**Usage**:

```console
$ rpv.py plasmid [OPTIONS] DATABASE SEQUENCES...
```

**Arguments**:

* `DATABASE`: PlasmidFinder Database format fasta.  [required]
* `SEQUENCES...`: Séquences au format fasta.  [required]

**Options**:

* `--output PATH`: Chemin du fichier de sortie du résultat (.tsv).
* `--identity INTEGER RANGE`: Pourcentage minimum d'indentité de séquence  [default: 90]
* `--help`: Show this message and exit.

## `rpv.py resistance`

Search resistance genes

**Usage**:

```console
$ rpv.py resistance [OPTIONS] DATABASE SEQUENCES...
```

**Arguments**:

* `DATABASE`: ResFindder Database au format fasta.  [required]
* `SEQUENCES...`: Séquences au format fasta.  [required]

**Options**:

* `--output PATH`: Chemin du fichier de sortie du résultat (.tsv).
* `--identity INTEGER RANGE`: Pourcentage minimum d'indentité de séquence  [default: 90]
* `--help`: Show this message and exit.

## `rpv.py virulence`

Search virulence genes

**Usage**:

```console
$ rpv.py virulence [OPTIONS] DATABASE SEQUENCES...
```

**Arguments**:

* `DATABASE`: VFDB (Virulence Factor Database) au format fasta  [required]
* `SEQUENCES...`: Séquences au format fasta.  [required]

**Options**:

* `--output PATH`: Chemin du fichier de sortie du résultat (.tsv).
* `--identity INTEGER RANGE`: Pourcentage minimum d'indentité de séquence  [default: 90]
* `--help`: Show this message and exit.

#!/usr/bin/python
# -*- coding: utf-8 -*-

##  -------  Descriptions ------- ##
##
## Date : 22/02/2022
## Auteur : Atchon Alban
## email : alban.atchon@gmail.com
## github : https://github.com/ATCHON
## Version : 0.1.0
##
## ------------------------------- ##

__version__ = "0.1.0"


import shlex
import subprocess
from typing import List, Optional

import pandas as pd
import typer
from pathlib import Path

WK_DIR = Path(__file__).resolve().parent
TEMP_DIR = Path.cwd() / "temp"
TEMP_DIR.mkdir(exist_ok=True)
BLAT_EXE = WK_DIR / "content" / "src" / "blat"
FILTER_BLAST = WK_DIR / "content" / "src" / "filter_blast.py"


app = typer.Typer()


def version_callback(value: bool):
    if value:
        typer.echo(f"RPV Version: {__version__}")
        raise typer.Exit()


def run_blat(fasta, db, identity=int):
    """
    run_blat : Function running BLAT (Blast-Like Alignment Tool), for the identification of genes of interest.

    Args: fasta (Fasta file): Assembly sequence of the strain of interest.
          db (Multi-fasta file): Database => ResFinder database, PlasmidFinder database or VirulenceFactor database
    Returns: tsv, Résultat temporaire de l'identification des gènes via Blat.
    """

    blat_cmd = f"{BLAT_EXE} -fine -minIdentity={identity} -out=blast8 {db} {fasta} {TEMP_DIR / 'tmp.blat'}"
    subprocess.run(shlex.split(blat_cmd))

    flt_cmd = f"python3 {FILTER_BLAST} -r {db} {TEMP_DIR / 'tmp.blat'} -o {TEMP_DIR / 'tmp.tsv'}"
    subprocess.run(shlex.split(flt_cmd))
    return pd.read_csv(f"{TEMP_DIR / 'tmp.tsv'}", sep='\t')


def sort_dataframe(df):
    df = df[['pident'] + [col for col in df.columns if col != 'pident']]
    df = df[['gene'] + [col for col in df.columns if col != 'gene']]
    df = df[['Sample'] + [col for col in df.columns if col != 'Sample']]
    df = df.iloc[:, :-8]
    df = df.reset_index(drop=True)
    return df


def manage_resistance_df(df):
    df['gene'] = df['saccver'].str.split('_', expand=True)[0]
    df = df.drop(['saccver', 'descref'], axis=1)
    df = sort_dataframe(df)
    return df


def manage_virulence_df(df):
    df['gene'] = df['descref'].str.split(' ', expand=True, )[1].str.replace(r'(', '', regex=True).str.replace(r')', '',
                                                                                                              regex=True)
    df['description'] = df['descref'].str.split(' ', 2, expand=True)[2]
    df = df.drop(['saccver', 'descref'], axis=1)
    df = df[['description'] + [col for col in df.columns if col != 'description']]
    df = sort_dataframe(df)
    return df


def manage_plasmid_df(df):
    df['gene'] = df['saccver'].str.split('_', expand=True)[0]
    df['description'] = df['descref'].str.split('_', 2, expand=True)[2]
    df = df.drop(['saccver', 'descref'], axis=1)
    # df = df[['description'] + [col for col in df.columns if col != 'description']]
    df = sort_dataframe(df)
    return df


# @app.command('run')
def main(database, sequences, identity):
    """
    main : Execute search on database and add sample name on dataframe

    Args:
        identity: Sequences identity percent.
        database : Database of nucleotide sequences in fasta format without duplicate sequences.
        sequences : Nucleotide sequences from read assemblies.
    Returns: dataframe
    """
    df = None
    for nb_strain, sequence in enumerate(sequences, start=1):
        seq_name = Path(sequence).name
        typer.echo(f"Search genes for sample {seq_name} .... N°{nb_strain}")
        df_temp = run_blat(sequence, database, identity)
        df_temp['Sample'] = str(seq_name.split('.')[0])
        df = pd.concat([df, df_temp]) if df is not None else df_temp
        typer.echo('--' * 20)
    return df


@app.command(help="Search resistance genes")
def resistance(database: Path = typer.Argument(..., help="ResFinder Database fasta format."),
               output: Optional[Path] = typer.Option(None, help="Output file path of the result (.tsv)."),
               identity: int = typer.Option(90, max=100, help="Minimum percentage of sequence identity."),
               sequences: List[Path] = typer.Argument(..., help="Fasta format sequences.")):
    """
    Search for resistance genes via the Resfinder database
    """
    df = main(database=database, sequences=sequences, identity=identity)
    if df.empty is not True:  # type: ignore
        df = manage_resistance_df(df)
        typer.echo(df)
        if output:
            df.to_csv(Path(output), sep='\t')
    else:
        typer.echo("No resistance genes have been identified")


@app.command(help="Search plasmid genes")
def plasmid(database: Path = typer.Argument(..., help="PlasmidFinder Database fasta format."),
            output: Optional[Path] = typer.Option(None, help="Output file path of the result (.tsv)."),
            identity: int = typer.Option(90, max=100, help="Minimum percentage of sequence identity."),
            sequences: List[Path] = typer.Argument(..., help="Fasta format sequences.")):
    """ 
    Search for plasmid genes via the PlasmidFinder database 
    """
    df = main(database=database, sequences=sequences, identity=identity)
    if df.empty is not True:  # type: ignore
        df = manage_plasmid_df(df)
        typer.echo(df)
        if output:
            df.to_csv(Path(output), sep='\t')
    else:
        typer.echo("No plasmid gene has been identified")


@app.command(help="Search virulence genes")
def virulence(database: Path = typer.Argument(..., help="VFDB (Virulence Factor Database) fasta format."),
              output: Optional[Path] = typer.Option(None, help="Output file path of the result (.tsv)."),
              identity: int = typer.Option(90, max=100, help="Minimum percentage of sequence identity."),
              sequences: List[Path] = typer.Argument(..., help="Fasta format sequences.")):
    """
    Search for virulence genes via the Virulence Factor database 
    """
    df = main(database=database, sequences=sequences, identity=identity)
    if df.empty is not True:  # type: ignore
        df = manage_virulence_df(df)
        typer.echo(df)
        if output:
            df.to_csv(Path(output), sep='\t')
    else:
        typer.echo("No virulence genes have been identified")


@app.callback()
def version_call(version: Optional[bool] = typer.Option(
                None, "--version", callback=version_callback)):
    typer.echo("")


if __name__ == "__main__":
    app()

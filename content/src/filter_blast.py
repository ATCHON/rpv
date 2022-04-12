#!/usr/bin/python
# -*- coding: utf-8 -*-
# Auteur : Benoit Valot / github :https://github.com/bvalot

"""Filter Blast result (format ) to return best unique position"""

import argparse
import sys

import numpy
from Bio import SeqIO

desc = "Filter Blast result (format 6) to return best unique position for query"
command = argparse.ArgumentParser(prog='filter_blast.py',
                                  description=desc, usage='%(prog)s [options] blast')
command.add_argument('-o', '--out', nargs="?",
                     type=argparse.FileType("w"), default=sys.stdout,
                     help='Return filter blast file, default=stdout')
command.add_argument('-i', '--identity', nargs="?",
                     type=int, default=90,
                     help='Minimum percentage identity to return result, default 90')
command.add_argument('-c', '--coverage', nargs="?",
                     type=int, default=80,
                     help='Minimum percentage coverage to return result, work only if fasta is supply, default 80')
command.add_argument('--overlap', nargs="?",
                     type=int, default=80,
                     help='Minimum percentage overlap between to alignment to merge it, default 80')
command.add_argument('--fraction', nargs="?",
                     type=int, default=2,
                     help='Maximum percentage fraction of length difference to compare identity between 2 alignments, default 2')
command.add_argument('-r', '--reference', nargs="?",
                     type=argparse.FileType("r"),
                     help='Add reference fasta file for coverage check')
command.add_argument('-q', '--query', nargs="?",
                     type=argparse.FileType("r"),
                     help='Add query fasta file for coverage check')
command.add_argument('blast', type=argparse.FileType("r"),
                     help='Tabular blast result (format 6)')
command.add_argument('-v', '--version', action='version',
                     version='%(prog)s 0.3.0')


def read_description(fasta):
    if fasta is None:
        return None
    return {
        seq.id: (seq.description, len(seq))
        for seq in SeqIO.parse(fasta, 'fasta')
    }


def fraction(len1, len2):
    return float(numpy.abs(len1 - len2)) / numpy.mean([len1, len2])


class Align:
    """A simple Align class"""

    def __init__(self, value):
        self.value = value
        self.saccver = self.get("saccver")
        self.qaccver = self.get("qaccver")
        self.qstart = int(self.get("qstart"))
        self.qend = int(self.get("qend"))
        self.length = int(self.get("length"))
        self.pident = float(self.get("pident"))

    def overlap(self, align, overlap):
        if self.qaccver != align.qaccver:
            return False
        over = len(self.maps().intersection(align.maps()))
        if over == 0 or (float(over) / min(self.length, align.length)) * 100 < overlap:
            return False
        return True

    def get(self, key):
        return self.value.get(key)

    def maps(self):
        return set(range(self.qstart, self.qend))

    def coverage(self, totlength):
        return (float(self.length) / totlength) * 100

    def __len__(self):
        return self.length


class Aligns:
    """A class to compile alignment"""

    def __init__(self, overlap, fraction):
        self.overlap = overlap
        self.fraction = fraction
        self.queries = {}
        # self.alignment = []

    def add_align(self, align):
        """Choice best on similar length and pident"""
        als = self.queries.setdefault(align.qaccver, [])
        add = False
        for i, al in enumerate(als):
            if al.overlap(align, self.overlap):
                add = True
                if fraction(align.length, al.length) * 100 < self.fraction:
                    if align.pident > al.pident:
                        als[i] = align
                elif align.length > al.length:
                    als[i] = align
                # if align.length == al.length:
                #     if align.pident > al.pident:
                #         self.alignment[i] = align
                # elif align.length > al.length:
                #     self.alignment[i] = align
        if not add:
            als.append(align)

    def __iter__(self):
        for als in self.queries.values():
            yield from als

    def __len__(self):
        return sum(map(len, self.queries.values()))


# sourcery no-metrics skip: raise-specific-error
if __name__ == '__main__':
    """Performed job on execution script"""
    args = command.parse_args()
    output = args.out

    descref = read_description(args.reference)
    descquery = read_description(args.query)

    # header = ["chr", "database", "feature", "start", "stop", "ident", "strand", "other", "informations"]
    header = ["qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
              "sstart", "send", "evalue", "bitscore"]
    aligns = Aligns(args.overlap, args.fraction)
    for line in args.blast.readlines():
        if line[0] == "#":
            output.write(line)
            continue
        value = dict(zip(header, line.strip().split()))
        al = Align(value)

        ## identity
        if al.pident < args.identity:
            continue

        ## coverage
        lens = []
        if descref is not None:
            if descref.get(al.get('saccver')) is None:
                raise Exception("The sequence " + al.get('saccver') + " is not found in your fasta file")
            lens.append(descref.get(al.get('saccver'))[1])
        if descquery is not None:
            if descquery.get(al.get('qaccver')) is None:
                raise Exception("The sequence " + al.get('qaccver') + " is not found in your fasta file")
            lens.append(descquery.get(al.get('qaccver'))[1])

        if lens and al.coverage(min(lens)) < args.coverage:
            continue
        aligns.add_align(al)

    ## export results
    output.write("\t".join(header))
    if descref is not None:
        output.write("\tdescref")
    if descquery is not None:
        output.write("\tdescquery")
    output.write("\n")

    for al in aligns:
        towrite = [al.get(h) for h in header]
        if descref is not None:
            towrite.append(descref.get(al.get('saccver'))[0])
        elif descquery is not None:
            towrite.append(descquery.get(al.get('qaccver'))[0])
        else:
            towrite.append(".")
        output.write("\t".join(towrite) + "\n")

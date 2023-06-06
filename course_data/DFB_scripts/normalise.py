#!/usr/bin/env python

"""
A simple normalising function designed for use during part 3 of the WCSC de novo
practical session. Uses khmer to interleave the fastqs and jellyfish to count
the kmers.
Reads are retained if their median kmer frequency * random (0, 1) < C, currently
hard-coded as 10. The heuristic nature of this process is controlled by setting
random.seed(45)

input
-----
fq  two fastqs as positional arguments. Not gzipped. Intended to be
    trimmed_1P.fastq and trimmed_2P.fastq.

output
------
normalised_R1.fastq & normalised_R2.fastq
    These are the two normalised fastq files for downstream SPAdes analysis

"""
import argparse
import io
import os

from utils import generalUtilities as gU

from statistics import median
from random import random, seed

seed(45)
C = 10
nproc = int(gU.shell("nproc").stdout)

def median_proc(*read_pair):
    "If the median kmer frequency * random < C, keep the read"
    if random() * median([
        kmers.get(read_pair[seq_pos][start:][:31], 1)
        for seq_pos in (1, 5)
        for start in range(len(read_pair[seq_pos]) - 31)
        ]) < C:
        
        return (read_pair[:4], read_pair[4:])
    
if __name__ == '__main__':
    
    # Obtain fastq paired files
    parser = argparse.ArgumentParser()
    parser.add_argument('fq', type=str, nargs=2)
    args = parser.parse_args()
    
    outfastqs = ('normalised_R1.fastq',
                 'normalised_R2.fastq')
    
    print("Interleaving FASTQs")
    interleaved = "interleaved.fastq"
    if not gU.exists(interleaved):
        gU.interleave_fastq(args.fq, o=interleaved)

    print("Counting kmers with Jellyfish")
    jf = "norm.jf"
    if not gU.exists(jf):
        gU.jf_count(
            "--out-counter-len", "2",
            m=31, t=8, s="100M", c=1, o=jf, fastq=interleaved
        )
    
    print("Dumping Jellyfish data into dict")
    kmers = dict(
        (s, int(n))
        for n, s in gU.fasta_parser(
            io.StringIO(gU.jf_dump(jf=jf, L=1).stdout)
        )
    )
    
    print(f"Normalising randomly (C={C})")
    read_set = [*zip(*filter(
        None,
        gU.threaded(
            func=median_proc,
            data=[*gU.iter_zip(8, open(interleaved).readlines())],
            procs=nproc
        )
    ))]
    
    print("Writing output files")
    for reads, outfastq in zip(read_set, outfastqs):
        with open(outfastq, 'w') as out:
            out.write("\n".join(gU.chain(reads)))
            
    gU.remove(jf, interleaved)
    print('Normalisation complete')
    

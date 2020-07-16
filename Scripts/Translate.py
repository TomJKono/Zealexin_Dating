#!/usr/bin/env python
"""A simple script to take a set of unaligned nucleotide CDS sequences,
translate them to amino acid sequences, align them with clustal-omega, then
back-translate them into nucleotides. Requires Biopython and clustal-omega.
Writes the aligned sequences to stdout. Takes one argument:
    1) Input FASTA

Author:
Thomas Kono (konox006@umn.edu)
Minneosta Supercomputing Institute
Minneapolis, MN"""

try:
    import sys
    import subprocess
    import tempfile
    from Bio import SeqIO
    from Bio import SeqRecord
    from Bio import Seq
except ImportError:
    sys.stderr.write('This script requires Biopythonn.\n')
    exit(1)


def nuc_to_aa(s):
    """Translate the nucleotide sequences and write them to a temporary file.
    Return the path to the temporary file."""
    translated = tempfile.NamedTemporaryFile(
        prefix='Translated_',
        suffix='.fa',
        mode='w+t')
    nuc_seqs = {}
    for seq in SeqIO.parse(s, 'fasta'):
        s = SeqRecord.SeqRecord(
            seq.seq.translate(),
            id=seq.id,
            description='')
        SeqIO.write(s, translated, 'fasta')
    translated.seek(0)
    return translated


def main(fasta):
    """Main function."""
    # Load up the sequence dictionary
    orig_nuc = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    tl_seq = nuc_to_aa(fasta)
    print(tl_seq.read())
    tl_seq.close()
    return


try:
    in_fa = sys.argv[1]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    exit(2)

main(in_fa)

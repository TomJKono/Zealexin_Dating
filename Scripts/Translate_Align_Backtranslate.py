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
    from Bio.Align.Applications import ClustalOmegaCommandline
    from Bio.Seq import IUPAC
except ImportError:
    sys.stderr.write('This script requires Biopython.\n')
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


def align(s):
    """Align the translated sequences. Returns the path to the output file
    with the alignment."""
    aligned = tempfile.NamedTemporaryFile(
        prefix='Aligned_',
        suffix='.fa',
        mode='w+t')
    # Align with clustal-omega
    co_cmd = ClustalOmegaCommandline(
        infile=s.name,
        outfile=aligned.name,
        seqtype='protein',
        force=True,
        iterations=10,
        distmat_full=True,
        distmat_full_iter=True)
    co_cmd()
    aligned.seek(0)
    return aligned


def backtranslate(a, orig_seqs):
    """Back-translate the alignment using the original nucleotide sequences as
    a guide. Return the list of SeqRecords with the backtranslated seqs."""
    bt = []
    aln = SeqIO.parse(a.name, 'fasta')
    for prot_seq in aln:
        entry = ''
        codon = 0
        nuc = orig_seqs[prot_seq.id]
        for aa in prot_seq:
            if aa == '-':
                entry += '---'
            else:
                entry += str(nuc.seq[codon*3:(codon*3)+3])
                codon += 1
        bt.append(
            SeqRecord.SeqRecord(
                Seq.Seq(entry),
                id=prot_seq.id,
                description='')
            )
    # Close the alignment handle
    a.close()
    return bt


def main(fasta):
    """Main function."""
    # Load up the sequence dictionary
    orig_nuc = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    tl_seq = nuc_to_aa(fasta)
    aligned_seq = align(tl_seq)
    bt_seq = backtranslate(aligned_seq, orig_nuc)
    # Write it to stdout
    SeqIO.write(bt_seq, sys.stdout, 'fasta')
    # Close handles
    tl_seq.close()
    aligned_seq.close()
    return


try:
    in_fa = sys.argv[1]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    exit(2)

main(in_fa)

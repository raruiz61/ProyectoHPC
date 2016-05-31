#! /usr/bin/env python

import build_seq
from Bio import SeqIO
document = "_document"
for seq_record in SeqIO.parse("file.fasta", "fasta"):
    build_seq.build_seq(seq_record.seq)
    cmd.select(document,"all")
    cmd.save (seq_record.id+".pdb", document, -1, 'pdb')
    cmd.delete("all")

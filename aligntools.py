#!/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

import itertools

# yields a copy of seq, where elements not equal mask are replaced with elements from other_seq
# there should be the same number of elements in other_seq as there are unmasked elements in seq
def copy_masked(seq, other_seq, mask):
	it = iter(seq)
	for other_el in other_seq:
		while it.next() == mask:
			yield mask
		yield other_el
	try:
		while True:
			yield it.next()
	except StopIteration:
		pass

def apply_alignment(aligned_seq, new_seq):
	return ''.join(copy_masked(aligned_seq, new_seq, "-"))

def mapping_to_alignment(mapping, len1, len2, count_base=1):
	# note: by default, all indexes we talk about here are 1-based
	i1 = count_base
	i2 = count_base

	alignment = list()
	for res1, res2 in mapping:
		while i1 < res1:
			alignment.append( (i1, None) )
			i1 += 1
		while i2 < res2:
			alignment.append( (None, i2) )
			i2 += 1
		alignment.append( (res1, res2) )
		i1 += 1
		i2 += 1
	while i1 <= len1:
		alignment.append( (i1, None) )
		i1 += 1
	while i2 <= len2:
		alignment.append( (None, i2) )
		i2 += 1

	return alignment


def alignment_fasta(alignment, fasta1, fasta2):
	seq1, seq2 = zip(*alignment)
	seq1 = ''.join(( "-" if i is None else fasta1[i-1] for i in seq1 ))
	seq2 = ''.join(( "-" if i is None else fasta2[i-1] for i in seq2 ))
	return (seq1, seq2)


if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(description='Applies an existing alignment to an unaligned sequence')
	parser.add_argument('alignment', type=argparse.FileType('r'), help='MSA output in FASTA format')
	parser.add_argument('fasta', type=argparse.FileType('r'), nargs='+', help='FASTA file(s) with one sequence matching each sequence from the MSA')
	args = parser.parse_args()

	# read corresponding sequences together from all files mentioned, the aligned file first
	for entries in itertools.izip(*map(fasta2line, [args.alignment] + args.fasta)):
		head, seq = entries[0]
		for other_head, other_seq in entries[1:]:
			print apply_alignment(seq, other_seq)

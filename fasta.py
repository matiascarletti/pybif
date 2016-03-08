#!/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

def fasta_seqs(lines, strip_leading_bracket=False):
	head = ""
	for line in lines:
		if line.startswith(">"):
			if (len(head)):
				yield (head, seq)
			head = line.strip()
			if strip_leading_bracket:
				head = head[1:] # exclude leading '>' symbol
			seq = ""
		else:
			seq += line.strip()
	if (len(head)):
		yield (head, seq)

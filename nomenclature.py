#!/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

_aa1to3_map = None
_aa3to1_map = None


def _init_aa1to3_map():
	global _aa1to3_map;
	if _aa1to3_map is None:
		_aa1to3_map = {
			"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU",
			"F": "PHE", "G": "GLY", "H": "HIS", "I": "ILE",
			"K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
			"P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER",
			"T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
		}


def reverse_lookup(lookup):
	result = dict()
	for key in lookup:
		value = lookup[key]
		result[value] = key
	return result


def _init_aa3to1_map():
	global _aa3to1_map
	if _aa3to1_map is None:
		_init_aa1to3_map()
		_aa3to1_map = reverse_lookup(_aa1to3_map)


def aa3to1(aa_3):
	_init_aa3to1_map()
	return _aa3to1_map[aa_3]


def aa3_seq_to_fasta(aa_3_seq):
	return (aa3to1(aa) for aa in aa_3_seq)


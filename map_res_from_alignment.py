#!/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

import sys
import argparse
import operator
import string

parser = argparse.ArgumentParser(description='Given an .alignment file, map residue numbers between the two proteins.')
parser.add_argument('res', metavar='RES', type=str, nargs='+', help='One or more residue numbers')
parser.add_argument('--reverse', action='store_true', help='Reverse the mapping. Map residues numbers from the second protein to the first (the default is from the first to the second).')
args = parser.parse_args()

coords_started = False
get_cols = operator.itemgetter(0, 5) # native and align
res_set = set(map(string.strip, args.res))

for line in sys.stdin:
	# the residue and coordinate map goes on from the line 
	# 'alignment:' until the end of the file
	if line.startswith('alignment:'):
		coords_started = True
	elif coords_started:
		res_a, res_b = map(string.strip, get_cols(line.split('\t')))
		if (not args.reverse) and (res_a in res_set):
			print res_b
		elif (args.reverse) and (res_b in res_set):
			print res_a

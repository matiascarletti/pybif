#!/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

import operator
from collections import defaultdict
import sys

def pdb_res_int(s):
	try:
		return int(s)
	except ValueError:
		return int(s[:-1])

def parse_pdb_range(s):
	p = s[1:].split('-', 1)
	if len(p) == 1:
		return (pdb_res_int(s), pdb_res_int(s))
	else:
		return (pdb_res_int(s[0] + p[0]), pdb_res_int(p[1]))

def parse_seq_range(s):
	p = s.split('-', 1)
	if len(p) == 1:
		return (int(s), int(s))
	else:
		return (int(p[0]), int(p[1]))


def ecod_dom_pdb_ranges(dom):
	for chain_range in dom.split(","):
		chain, dom_range = chain_range.split(':', 1)
		try:
			yield chain, parse_pdb_range(dom_range)
		except ValueError:
			yield chain, None


def ecod_dom_seq_ranges(dom):
	for chain_range in dom.split(","):
		chain, dom_range = chain_range.split(':', 1)
		try:
			yield chain, parse_seq_range(dom_range)
		except ValueError:
			yield chain, None


def ecod_pdb_ranges_map(lines, verbose=False):
	get_cols = operator.itemgetter(4,6)
	db = defaultdict(lambda: defaultdict(list))

	for line in lines:
		if line.startswith('#'):
			continue

		pdb_id, dom = get_cols(line.strip().split('\t'))
		for chain, dom_range in ecod_dom_pdb_ranges(dom):
			if dom_range is not None:
				db[pdb_id][chain].append( dom_range )
			elif verbose:
				sys.stderr.write('Warning: ECOD line error, skipping domain: ' + line)

	return db

def domain_range_containing(dom_ranges, start_res, end_res):
	for dom_start_res, dom_end_res in dom_ranges:
		if start_res >= dom_start_res and end_res <= dom_end_res:
			return (dom_start_res, dom_end_res)

	return None

#!/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

# queries pdbaa.fasta.gz for sequences by pdb_id and chain

import gzip
import sys
import argparse
import os

from argparsetools import 
import fasta

def query_pdb_sequences(pdbaa_fasta, query_chains, query_pdbs = None):
	query_chains = set(query_chains)
	if query_pdbs == None:
		query_pdbs = set( pdb_id for pdb_id, chain in query_chains )
	for name, seq in fasta.fasta_seqs(pdbaa_fasta):
		pdb_infos = ( pdb_info.split(":", 2) for pdb_info in name[1:].split("|") )
		# first check if any of the PDBs appear in the query pdbs
		for pdb_id, chains, mol_type in pdb_infos:
			pdb_id = pdb_id.lower()
			if not pdb_id in query_pdbs:
				continue
			# if so, search for exact PDB-chain pair matches in query_chains
			for chain in chains.split(","):
				if (pdb_id, chain) in query_chains:
					query_chains.remove((pdb_id, chain))
					yield (pdb_id, chain, seq)
					if not query_chains:
						break
			if not query_chains:
				break
		if not query_chains:
				break

def main():
	parser = argparse.ArgumentParser(description='Get sequences from a PDBAA FASTA file given a list of PDB IDs')
	parser.add_argument('--pdbaa-fasta', type=optionally_compressed_filetype('r'), help='an (optionally gzip-compressed) PDBAA FASTA file. May also be specified by the PDBAA_FASTA environment variable. The PDBAA FASTA file is expected to have FASTA headers in the format >pdb_id:chain[,...]:polymer_type[|...]')
	parser.add_argument('pdb_chain', nargs='*', help='a list of PDB chain identifiers in the format pdb_id:chain; specify - to read a tab-separated table (rather than colons) from stdin. The PDB IDs are case-insensitive, but the chain identifiers are case-sensitive.')
	args = parser.parse_args()

	if args.pdbaa_fasta is None:
		pdbaa_fasta_filename = os.getenv('PDBAA_FASTA')
		if pdbaa_fasta_filename is None:
			parser.error('--pdbaa-fasta is not specified and environment variable PDBAA_FASTA is not defined.')
		args.pdbaa_fasta = opt_compressed_filetype('r')(os.getenv('PDBAA_FASTA'))

	queries = set(args.pdb_chain)
	if len(args.pdb_chain) == 1 and args.pdb_chain[0] == "-":
		queries = set( (tuple(line.strip().split("\t", 1)) for line in sys.stdin) )
	else:
		queries = set([ tuple(pdb_chain.split(":")) for pdb_chain in args.pdb_chain ])

	for pdb_id, chain, seq in query_pdb_sequences(args.pdbaa_fasta, queries):
		print ">" + pdb_id + ":" + chain
		print seq


if __name__ == "__main__":
	main()

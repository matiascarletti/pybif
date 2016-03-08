#!/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

import os
from lxml import etree
import aligntools
import pdb
import textwrap
import sys
import os

_fatcat = None

def get_fatcat():
	global _fatcat
	if _fatcat is None:
		_fatcat = os.getenv('FATCAT_EXEC', 'fatcat')
		if not os.path.isfile(_fatcat) or not os.access(_fatcat, os.X_OK):
			raise RuntimeError('FATCAT path %s is not executable' % (_fatcat,))
	return _fatcat

class FatcatRunError(Exception):
	pass

def fatcat_raw(file1, file2):
	fatcat_exec = get_fatcat()
	# the first three lines of output are ignored since they are just logging information from FATCAT and not XML
	command = '%s -file1 %s -file2 %s -printXML 2> /dev/null | tail -n+4' \
        % (fatcat_exec, file1, file2)
	xml_string = os.popen(command).read()
	if len(xml_string) == 0:
		raise FatcatRunError('Error running command: %s' % command)
	return xml_string

def fatcat_xml(file1, file2):
	xml_string = fatcat_raw(file1, file2)
	return etree.fromstring(xml_string)
		

#def fatcat(file1, file2):
#	root = fatcat_xml(file1, file2)
#	result = dict(root.attrib)
#	blocks = []
#	for block in root.iter('block'):
#		block = dict(block.attrib)
#		eqrs = []
#		for eqr in root.iter('eqr'):
#			eqrs.append(dict(eqr.attrib))
#		block["eqrs"] = eqrs
#	result["blocks"] = blocks
#	return result

def matrix_and_shift_from_fatcat_rigid(root):
	return matrix_and_shift_from_fatcat_block(next(root.iter("block")))

def matrix_and_shift_from_fatcat_block(block):
	matrix = next(block.iter("matrix"))
	shift = next(block.iter("shift"))
	# If
	# - u is the template vector
	# - v is the target vector
	# - M is the matrix listed
	# - s is the shift vector listed
	# then (u'*M+s) is the superposition of u against v.
	# So we return the transpose of M, which will allow to multiply u by M on the left to get M*u+s ~ v
	return (
		[ [ float(matrix.get("mat%d%d" % (column+1,row+1))) for column in range(3) ] for row in range(3) ],
		[ float(shift.get(key))  for key in ("x", "y", "z") ],
		)


# this assumes that there is one block in the fatcat alignment (which I currently understand
# to be true for rigid alignments), and that each structure has a single chain.	
def fatcat_pdb_num_mapping(root):
	block = next(root.iter("block"))
	for eqr in block.iter("eqr"):
		yield (eqr.get("pdbres1"), eqr.get("pdbres2"))

def threeway_alignment(args, no_apply=False):
	pdb1_filename, pdb2_filename, apply_transform_to_filename = args.pdb1, args.pdb2, args.pdb3
	xml_raw = fatcat_raw(pdb1_filename, pdb2_filename)
	root = etree.fromstring(xml_raw)
	matrix, shift_vec = matrix_and_shift_from_fatcat_rigid(root)
	pdb2_cas = list(pdb.ca_atoms(pdb.pdb_atoms(open(pdb2_filename, 'r'))))
	# vecs = ( atom["pos"] for atom in pdb2_cas )
	# cm_vec = [ float(sum(coord_vals))/len(coord_vals) for coord_vals in zip(*vecs) ]

	if ( not no_apply ):
		atoms = list(pdb.pdb_atoms(open(apply_transform_to_filename, 'r')))
		atoms = pdb.apply_transform(atoms, matrix, shift_vec_before=(0,0,0), shift_vec_after=shift_vec)
		for line in pdb.output_atoms(atoms):
			print line

	if args.output_alignment != None:
		with open(args.output_alignment, 'w') as out_file:
			out_file.write(xml_raw)

	if args.output_fasta != None:
		pdb1_cas = list(pdb.ca_atoms(pdb.pdb_atoms(open(pdb1_filename, 'r'))))
		pdb1_map = { atom["res_id"].strip(): index for index, atom in enumerate(pdb1_cas, 1) }
		pdb2_map = { atom["res_id"].strip(): index for index, atom in enumerate(pdb2_cas, 1) }
		# get mapping, but translate pdb numbering to sequential numbering
		mapping = [ (pdb1_map[res1], pdb2_map[res2]) for (res1, res2) in fatcat_pdb_num_mapping(root) ]
		alignment = aligntools.mapping_to_alignment(mapping, len(pdb1_cas), len(pdb2_cas))
		fasta1, fasta2 = aligntools.alignment_fasta( alignment, pdb.to_fasta(pdb1_cas), pdb.to_fasta(pdb2_cas) )
		with open(args.output_fasta, 'w') as out_file:
			out_file.write(textwrap.dedent('''
				>{name1}
				{seq1}
				>{name2}
				{seq2}{newline}
			''').strip().format(name1=pdb1_filename, seq1=fasta1, name2=pdb2_filename, seq2=fasta2, newline='\n'))


def twoway_alignment(args):
	args.pdb3 = args.pdb2
	threeway_alignment(args, no_apply=args.no_apply)

def main():
	import argparse

	# --- copied from mammoth - consider refactoring ---
	parser = argparse.ArgumentParser()

	parser.add_argument('--output-alignment', metavar='FILE', type=str, help='filename to write alignment data in this script\'s verbose custom format, which includes the transformation matrix and vectors, as well as the aligned coordinates side-by-side.')
	parser.add_argument('--output-fasta', metavar='FILE', type=str, help='filename to write alignment data in FASTA-format')

	subparsers = parser.add_subparsers(help='Sub-command help')

	parser_3way = subparsers.add_parser('threeway', help='Find transformation from pdb2 to pdb1, and apply to pdb3')
	parser_3way.add_argument('pdb1', type=str, help='Target model')
	parser_3way.add_argument('pdb2', type=str, help='Mobile model to use for alignment')
	parser_3way.add_argument('pdb3', type=str, help='Mobile model onto which the transformation will be applied')
	parser_3way.set_defaults(func=threeway_alignment)

	parser_2way = subparsers.add_parser('align', help='Find transformation from pdb2 to pdb1, and apply to pdb2')
	parser_2way.add_argument('pdb1', type=str, help='Target model')
	parser_2way.add_argument('pdb2', type=str, help='Mobile model')
	parser_2way.add_argument('--no-apply', action='store_true', help='Do not actually apply the transformation, just find it, and output it if --output-alignment is specified.')
	parser_2way.set_defaults(func=twoway_alignment)

	args = parser.parse_args()

	args.func(args)
	# --- end copied from mammoth --

if __name__ == "__main__":
	main()

#!/usr/bin/python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

## make mammoth structure alignments
# here: find transformation based on pdb1 and pdb2; apply to pdb3
# fixed: hetatm/atm: first siz chars are copied from original line
# rewritten from superimpose_part.py

import string
from glob import glob
import os
import tempfile
import shutil
import re
import sys
from collections import OrderedDict
import argparse
import json
import textwrap
import math

import aligntools
import nomenclature
import pdb

class MaxsubParseError(Exception):
	pass

def parse_maxsub_sup(input):
	exptl_coords = list()
	pred_coords = list()
	lines = input.readlines()
	chain_num = 1
	for line_num, line in enumerate(lines):
		if line.startswith("REMARK Transformation Matrix:"):
			matrix = [ map(float, string.split(row)[1:]) for row in lines[(line_num+1):(line_num+4)] ]
		elif line.startswith("REMARK Translation vector (Prediction):"):
			pred_vec = map(float, string.split(lines[line_num+1])[1:])
		elif line.startswith("REMARK Translation vector (Experiment):"):
			exptl_vec = map(float, string.split(lines[line_num+1])[1:])
		elif line.startswith("TER"):
			chain_num += 1
		elif line.startswith("ATOM  "):
			x = float(line[30:38])
			y = float(line[38:46])
			z = float(line[46:54])
			# chain = line[21] is blank for maxsub_sup
			res_num = line[22:26]
			aa = line[17:20]
			ins_code = line[26]
			res_id = res_num + ins_code
			if chain_num == 1:
				pred_coords.append((res_id, aa, x, y, z))
			elif chain_num == 2:
				exptl_coords.append((res_id, aa, x, y, z))
			else:
				raise MaxsubParseError('More than two chains in maxsub_sup.pdb file')

	return dict(matrix=matrix, pred_vec=pred_vec, exptl_vec=exptl_vec, pred_coords=pred_coords, exptl_coords=exptl_coords)

class RasmolParseError(Exception):
	pass

def parse_rasmol_tcl(rasmol_tcl_file):

	in_npairs = False
	pred_exptl_map = list()
	for line in rasmol_tcl_file:
		if line.startswith("# npairs ="):
			npairs = int(line[12:])
			in_npairs = True
		elif line.startswith("#") and in_npairs:
			# maps residues between experimental and prediciton models
			# the residue numbers used here are sequential (according to appearance in ATOM sequence)
			res_pred = int(line[2:6])
			res_exptl = int(line[6:10])
			pred_exptl_map.append((res_pred, res_exptl))

	if npairs != len(pred_exptl_map):
		raise RasmolParseError('rasmol.tcl files has npairs=%d, but %d aligned residues were listed' % (npairs, len(pred_exptl_map)))

	return pred_exptl_map

def apply_transform(model_file, matrix, pred_vec, exptl_vec, outfile, renumber_atoms = False):
	# TODO : make use of the pdb library
	atom_count = 0
	for line in model_file:
		if line[:6] in ['ATOM  ','HETATM']:
			atom_count += 1
			x = float(line[30:38])
			y = float(line[38:46])
			z = float(line[46:54])
			xt, yt, zt = apply_transform_to_vec((x,y,z), matrix, pred_vec, exptl_vec)

			if renumber_atoms:
				atom_number = str(atom_count).rjust(5)
			else:
				atom_number = line[6:11]

			outfile.write('%-6s%s%s%8.3f%8.3f%8.3f%s\n'\
				  %(line[0:6],atom_number,line[11:30],xt,yt,zt,line[54:-1]))
	
		elif line[:6] == 'ENDMDL':
			break

	outfile.write('TER\n')


def apply_transform_to_vec(vec,matrix,pred_vec,exptl_vec):
	ans = [0.0]*3
	for i in range(3):
		for j in range(3):
			ans[i] = ans[i] + matrix[i][j]*(vec[j]-exptl_vec[j])
		ans[i] = ans[i] + pred_vec[i]

	return ans

class MammothAlignError(Exception):
	pass

_mastodon = None

def get_mastodon():
	global _mastodon
	if _mastodon is None:
		_mastodon = os.getenv('MASTODON_EXEC', 'mastodon')
	return _mastodon


def align(pred_filename, exptl_filename, retain_files = False, last_residue_fix=True):
	model_count = 0
	atom_count = 0

	if not retain_files:
		pred_filename = os.path.abspath(pred_filename)
		exptl_filename = os.path.abspath(exptl_filename)
		orig_cwd = os.getcwd()
	try:
		if not retain_files:
			temp = tempfile.mkdtemp()
		try:
			if not retain_files:
				os.chdir(temp)

			mastodon = get_mastodon()

			stdout = os.popen('%s -p %s -e %s -r 1 2> /dev/null' \
				% (mastodon, pred_filename,exptl_filename)).read()
	
			if len(stdout) == 0:
				raise MammothAlignError('(pred=%s, exptl=%s): No output. Is %s executable? Are outputs empty?' % (pred_filename, exptl_filename, mastodon))

			try:
				trans_data = parse_maxsub_sup(open('maxsub_sup.pdb', 'r'))
			except IOError as e:
				raise MammothAlignError('(pred=%s, exptl=%s): IOError (%s) while opening maxsub_sup.pdb' % (pred_filename, exptl_filename, str(e)))

			try:
				pred_exptl_map = parse_rasmol_tcl(open('rasmol.tcl', 'r'))
			except IOError as e:
				raise MammothAlignError('(pred=%s, exptl=%s): IOError (%s) while opening rasmol.tcl' % (pred_filename, exptl_filename, str(e)))

			matrix, pred_vec, exptl_vec, pred_coords, exptl_coords = (trans_data[key] for key in ("matrix", "pred_vec", "exptl_vec", "pred_coords", "exptl_coords"))

			# this was a nice idea, but in cases where there are insertions codes (e.g. 1goaA), the mapping fails

			#try:
			#	exptl_map = { exptl_coords.keys()[res_exptl-1]: pred_coords.keys()[res_pred-1] for res_exptl, res_pred in pairs }
			#except IndexError as e:
			#	raise MammothAlignError('(pred=%s, exptl=%s): IndexError (%s) while matching the residue map from rasmol.tcl to the actual residues of the experimental structure in maxsub_sup.pdb' % (pred_filename, exptl_filename, str(e)))

			#try:
			#	pred_map = { pred_coords.keys()[res_pred-1]: exptl_coords.keys()[res_exptl-1] for res_exptl, res_pred in pairs }
			#except IndexError as e:
			#	raise MammothAlignError('(pred=%s, exptl=%s): IndexError (%s) while matching the residue map from rasmol.tcl to the actual residues of the predicted structure in maxsub_sup.pdb' % (pred_filename, exptl_filename, str(e)))

			tokens = [ re.split("\s+", line) for line in stdout.splitlines() if line.startswith(' PSI(end)') ][0]
			rmsd = tokens[8]
			n_ali = tokens[4]
	
			if last_residue_fix:
				pred_atoms = list(pdb.ca_atoms(pdb.pdb_atoms(open(pred_filename, 'r'))))
				if len(pred_atoms) == len(pred_coords)+1:
					# fix last residue being missing from mammoth alignment
					exptl_atoms = list(pdb.ca_atoms(pdb.pdb_atoms(open(exptl_filename, 'r'))))
					last_pred_ca = pred_atoms[-1]
					last_exptl_ca = exptl_atoms[-1]
					res_id1, aa1, x1, y1, z1 = ( last_pred_ca[key] for key in ("res_id", "res_name", "x", "y", "z") )
					x1, y1, z1 = (x1-pred_vec[0],y1-pred_vec[1],z1-pred_vec[2])
					res_id2, aa2, x2, y2, z2 = ( last_exptl_ca[key] for key in ("res_id", "res_name", "x", "y", "z") )
					x2, y2, z2 = apply_transform_to_vec((x2,y2,z2), matrix, (0,0,0), exptl_vec)

					pred_coords.append((res_id1, aa1, x1, y1, z1))
					exptl_coords.append((res_id2, aa2, x2, y2, z2))

					# check if the added residue continues an alignment stretch
					if pred_exptl_map[-1][0] == len(pred_atoms)-1:
						# if the pred model res is added to a stretch, try to align against the following exptl model res
						res_id2, aa2, x2, y2, z2 = exptl_coords[pred_exptl_map[-1][1]]
					elif pred_exptl_map[-1][1] == len(exptl_atoms)-1:
						# otherwise, the pred model res is not the last aligned residue, but exptl may be aligned against a previous residue in pred model
						res_id1, aa1, x1, y1, z1 = pred_coords[pred_exptl_map[-1][0]]

					distance_sq = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
					if distance_sq < 16: # angstrom squared
						pred_exptl_map.append( (len(pred_atoms), len(exptl_atoms)) )
						rmsd = float(rmsd)
						n_ali = int(n_ali)
						rmsd = '%.2f' % math.sqrt( ((rmsd**2 * n_ali) + distance_sq) / float(n_ali+1) )
						n_ali = str(n_ali + 1)

			return dict(matrix=matrix, pred_vec=pred_vec, exptl_vec=exptl_vec, 
					exptl_coords=exptl_coords, pred_coords=pred_coords, 
					pred_exptl_map=pred_exptl_map, rmsd=rmsd, n_ali=n_ali)
	
		finally:
			if not retain_files:
				shutil.rmtree(temp)
	finally:
		if not retain_files:
			os.chdir(orig_cwd)


def aligned_coords_to_fasta( coords, res_nums ):
	# here is another implementation:
	# seq = [ aa for res_id, aa, x, y, z in coord ] + ["-"]
	# pos = [ res_num-1 if res_num != None else len(seq)-1 for res_num in res_nums ]
	# return "".join(( nomenclature.aa_3to1(aa) for aa in operator.itemgetter(pos)(seq) ))
	result = list()
	for res_num in res_nums:
		if res_num == None: # non-aligned residue; treat as gap
			result.append('-')
		else:
			res_id, aa, x, y, z = coords[res_num-1]
			result.append(nomenclature.aa3to1(aa))

	return ''.join(result)

def gen_coord_alignment(pred_filename, exptl_filename, align_data, alignment = None):
	pred_coords, exptl_coords, pred_exptl_map = (align_data[key] for key in ("pred_coords", "exptl_coords", "pred_exptl_map"))
	if alignment == None:
		alignment = aligntools.mapping_to_alignment(pred_exptl_map, len(pred_coords), len(exptl_coords))

	# each item in pred_coords or exptl_coords is the tuple (res_id, aa, x, y, z)
	# we expect items to be indexed by the res_num (minus 1, since indexing in Python is zero-based)
	# we produce empty-string 4-tuples (well, 4-lists, because it's easier) in place of gaps
	alignment = '\n'.join(
			'\t'.join(
				(map(str, pred_coords[pred_res_num-1]) if pred_res_num != None else [""]*4) +
				(map(str, exptl_coords[exptl_res_num-1]) if exptl_res_num != None else [""]*4)
			) for pred_res_num, exptl_res_num in alignment
		)

	# note: dedent() and .strip() are used so the heredoc can be nicely aligned with rest of code
	return textwrap.dedent('''
		Alignment data
		# pred is the target, stationary model
		# exptl is superimposed to pred
		pred: {pred_filename}
		exptl: {exptl_filename}
		
		# transform rule is pred_transl_vec + matrix * (vec - exptl_transl_vec)
		# for each vector vec in the exptl model
		
		matrix: {matrix}
		pred_transl_vec: {pred_vec}
		exptl_transl_vec: {exptl_vec}
		
		# following is the per-residue alignment of the models
		n_ali: {n_ali}
		rmsd: {rmsd}
		alignment:
		{alignment}
	''').strip().format(**dict(align_data, **locals()))

def gen_fasta_alignment(pred_filename, exptl_filename, align_data, alignment = None):
	pred_coords, exptl_coords, pred_exptl_map = (align_data[key] for key in ("pred_coords", "exptl_coords", "pred_exptl_map"))
	if alignment == None:
		alignment = aligntools.mapping_to_alignment(pred_exptl_map, len(pred_coords), len(exptl_coords))

	pred_aligned_seq = aligned_coords_to_fasta(pred_coords, (pred_res_num for pred_res_num, exptl_res_num in alignment))
	exptl_aligned_seq = aligned_coords_to_fasta(exptl_coords, (exptl_res_num for pred_res_num, exptl_res_num in alignment))

	return textwrap.dedent('''
			>{pred_filename} {n_ali} {rmsd}
			{pred_aligned_seq}
			>{exptl_filename}
			{exptl_aligned_seq}
		''').strip().format(**dict(align_data, **locals()))

def threeway_alignment(args, no_apply=False):
	pred_filename, exptl_filename, apply_transform_to_filename = args.pdb1, args.pdb2, args.pdb3
	align_data = align(pred_filename, exptl_filename, retain_files=args.keep_files)

	matrix, pred_vec, exptl_vec, rmsd, n_ali = (align_data[key] for key in ("matrix", "pred_vec", "exptl_vec", "rmsd", "n_ali"))
	if ( not no_apply ):
		with open(apply_transform_to_filename, 'r') as apply_transform_to_file:
			apply_transform(apply_transform_to_file, matrix, pred_vec, exptl_vec, outfile=sys.stdout, renumber_atoms = False)

	if args.output_alignment != None:
		with open(args.output_alignment, 'w') as out_file:
			out_file.write(gen_coord_alignment(pred_filename, exptl_filename, align_data))
			out_file.write('\n') # ensure trailing line

	if args.output_fasta != None:
		with open(args.output_fasta, 'w') as out_file:
			out_file.write(gen_fasta_alignment(pred_filename, exptl_filename, align_data))
			out_file.write('\n') # ensure trailing line

def twoway_alignment(args):
	args.pdb3 = args.pdb2
	threeway_alignment(args, no_apply=args.no_apply)

def main():

	parser = argparse.ArgumentParser()

	parser.add_argument('--keep-files', action='store_true', help='Keep files (maxsub_sup.pdb, maxsub_sup2.pdb, rasmol.tcl) produced by MAMMOTH')
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
	# TODO : incorporate stuff from superimpose.py

def orig_version():
	if len(sys.argv) <=2:
		print '\n'
		print '-'*75
		print 'USAGE: %s <pdb1> <pdb2> <pdb3> > <superposition-pdb>'%sys.argv[0]
		print '\n will superimpose pdb2 onto pdb1 and apply transformation to pdb3'
		print '-'*75
		print '\n\n'
		assert 0==1

	# apparently, MAMMOTH outputs the matrix that one would apply to the experimental
	# structure, in order to get them superimposed.

	pred_filename, exptl_filename, apply_transform_to_filename = sys.argv[1:4] 

	align_data = align(pred_filename, exptl_filename, retain_files=True)
	matrix, pred_vec, exptl_vec, rmsd, n_ali = (align_data[key] for key in ("matrix", "pred_vec", "exptl_vec", "rmsd", "n_ali"))

	apply_transform(open(apply_transform_to_filename, 'r'), matrix, pred_vec, exptl_vec, outfile=sys.stdout, renumber_atoms = False)

	sys.stderr.write('%s -vs- %s: %s over %s residues\n' % (pred_filename, exptl_filename, rmsd, n_ali))

if __name__ == "__main__":
	main()

#!/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

import itertools
import sys
import math

import nomenclature
import filetools

# using a class is as fast as a dict, but uses about 1/9th of the memory
# can be used for HETATM as well as ATOM records
class PdbAtom(object):
	__slots__ = [
		'rec_name',
		'atom_num',
		'atom_name',
		'alt_loc',
		'res_name',
		'chain',
		'res_num',
		'ins_code',
		'x',
		'y',
		'z',
		'occupancy',
		'temp_factor',
		'element',
		'charge',
	]

	def __init__(self, line=None):
		if line is None:
			return

		self.rec_name = line[:6]
		self.atom_num = int(line[6:11])
		self.atom_name = line[12:16]
		self.alt_loc = line[16]
		self.res_name = line[17:20]
		self.chain = line[21]
		self.res_num = int(line[22:26])
		self.ins_code = line[26]
		self.x = float(line[30:38])
		self.y = float(line[38:46])
		self.z = float(line[46:54])

		l = len(line)
		if (l > 54):
			self.occupancy = float(line[54:60])
		if (l > 60):
			self.temp_factor = float(line[60:66])
		if (l > 76):
			self.element = line[76:78]
		if (l > 78):
			self.charge = line[78:80]

	def __getitem__(self, key):
		return getattr(self, key)

	def __setitem__(self, key, value):
		setattr(self, key, value)

	@property
	def res_id(self):
		# this includes the insertion code
		return '%4d%1s' % (self.res_num, self.ins_code)

	@res_id.setter
	def res_id(self, value):
		self.res_num = int(value[:4])
		self.ins_code = value[4]

	@property
	def pos(self):
		return (self.x, self.y, self.z)

	@pos.setter
	def pos(self, value):
		self.x, self.y, self.z = value

	def to_dict(self):
		return {key: getattr(self, key) for key in self.__slots__}

	def from_dict(self, d):
		for key in self.__slots__:
			if key in d:
				setattr(self, key, d[key])

	def from_atom(self, other):
		for key in self.__slots__:
			if hasattr(other, key):
				setattr(self, key, getattr(other, key))

	def __str__(self):
		# this is about two times quicker than "...".format(**self)
		return "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s" % (
			self.rec_name, self.atom_num, self.atom_name, self.alt_loc, self.res_name, self.chain, self.res_num,
			self.ins_code, self.x, self.y, self.z, self.occupancy, self.temp_factor, self.element, self.charge)


def pdb_atoms(lines, include_hetatms=False, allow_multi_model=False):
	# TODO : support ENDMDL, MODEL
	for line in lines:
		if not allow_multi_model and line.startswith('ENDMDL  '):
			break
		if not (line.startswith('ATOM  ') or (include_hetatms and line.startswith('HETATM'))):
			continue
		yield PdbAtom(line)


def ca_atoms(atoms):
	return (atom for atom in atoms if atom["atom_name"] == " CA ")


_backbone_atom_names = None
_backbone_atom_names_set = None
#_backbone_atom_order = None
_backbone_geometry = None
#_backbone_ncac_atom_names_set = None
#_backbone_co_atom_names_set = None
def _init_backbone_consts():
	global _backbone_atom_names, _backbone_atom_names_set, _backbone_geometry
#	global _backbone_atom_order, _backbone_co_atom_names_set, _backbone_ncac_atom_names_set
	_backbone_atom_names = [" N  ", " CA ", " C  ", " O  "]
	_backbone_atom_names_set = set(_backbone_atom_names)
#	_backbone_ncac_atom_names_set = set([" N  ", " CA ", " C  "])
#	_backbone_co_atom_names_set = set([" C  ", " O  "])
#	_backbone_atom_order = dict(((atom_name, index) for index, atom_name in enumerate(_backbone_atom_names)))
	_backbone_geometry = dict()
	# pairs tell distances in Angstroms
	_backbone_geometry[(" N  ", " CA ")] = 1.45
	_backbone_geometry[(" CA ", " C  ")] = 1.52
	_backbone_geometry[(" C  ", " N  ")] = 1.33
	_backbone_geometry[(" C  ", " O  ")] = 1.23
	

def backbone_atoms(atoms):
	if _backbone_atom_names is None:
		_init_backbone_consts()
	return (atom for atom in atoms if atom["atom_name"] in _backbone_atom_names)


#def ordered_backbone_atoms(atoms, single_chain=False):
#	if _backbone_atom_order is None:
#		_init_backbone_consts()
#	
#	return itertools.chain.from_iterable( ( 
#		sorted(res_atoms, key=lambda atom: _backbone_atom_order[atom["atom_name"]])
#		for res_id, res_atoms
#		in atoms_by_res(backbone_atoms(atoms), single_chain=single_chain) ) )


# from https://docs.python.org/2/libary/itertools.html
# later developed into http://stackoverflow.com/a/35137212/404861
def sliding_window(iterable, n):
    its = [ itertools.islice(iter, i, None) for i, iter in enumerate(itertools.tee(iterable, n)) ]
    return itertools.izip(*its)

def sliding_window_longest(iterable, n, fillvalue=None):
    its = [ itertools.islice(iter, i, None) for i, iter in enumerate(itertools.tee(iterable, n)) ]
    return itertools.izip_longest(*its, fillvalue=fillvalue)


def is_backbone_continuous(atoms, single_chain=False, relative_tolerance=0.1):
	if _backbone_geometry is None:
		_init_backbone_consts()

	def disconnected(atom1, atom2):
		distance = math.sqrt(atom_distance_sq(atom1, atom2))
		ideal_distance = _backbone_geometry[(atom1["atom_name"], atom2["atom_name"])]
		return (abs(distance - ideal_distance)/ideal_distance > relative_tolerance)

	for res1, res2 in sliding_window_longest( (list(res_atoms) for res_id, res_atoms in atoms_by_res(atoms, single_chain=single_chain)), n=2 ):
		res1 = { atom["atom_name"]: atom for atom in res1 }
		if res2 is not None:
			res2 = { atom["atom_name"]: atom for atom in res2 }
		print res1, res2
		if disconnected(res1[" N  "], res1[" CA "]) or disconnected(res1[" CA "], res1[" C  "]) or \
				disconnected(res1[" C  "], res1[" O  "]) or ( res2 is not None and disconnected(res1[" C  "], res2[" N  "]) ):
			return False

	return True


def atoms_by_res(atoms, single_chain=False, ignore_ins_code=False):
	if single_chain:
		if ignore_ins_code:
			return itertools.groupby(atoms, lambda (atom): atom["res_num"])
		else:
			return itertools.groupby(atoms, lambda (atom): atom["res_id"])
	else:
		if ignore_ins_code:
			return itertools.groupby(atoms, lambda (atom): atom["chain"] + str(atom["res_num"]))
		else:
			return itertools.groupby(atoms, lambda (atom): atom["chain"] + atom["res_id"])


def atoms_by_chain(atoms):
	return itertools.groupby(atoms, lambda (atom): atom["chain"])

def atom_distance_sq(atom1, atom2):
	return ( (atom1.x-atom2.x)**2 + (atom1.y-atom2.y)**2 + (atom1.z-atom2.z)**2 )

def chain_range(atoms, range_start=None, range_end=None, sequential_numbering=False):
	if sequential_numbering:
		# figure out sequential numbering for each atom
		if range_start is None:
			range_start = 1
		if range_end is None:
			range_end = sys.maxint

		for res_count, (res_num, res_atoms) in enumerate(atoms_by_res(atoms, single_chain=True), 1):
			if (res_count >= range_start) and (res_count <= range_end):
				for atom in res_atoms:
					yield atom
	else:
		# each atom already has correct residue number
		for atom in atoms:
			res_num = atom["res_num"]
			if (range_start is None or res_num >= range_start) and (range_end is None or res_num <= range_end):
				yield atom


def ranges(atoms, chain_ranges, sequential_numbering=False):
	for chain, chain_atoms in atoms_by_chain(atoms):
		if chain not in chain_ranges:
			continue
		range_res = chain_ranges[chain]
		if range_res is None:
			# yield all atoms for chain
			for atom in chain_atoms:
				yield atom
		else:
			range_start, range_end = range_res
			for atom in chain_range(chain_atoms, range_start, range_end, sequential_numbering):
				yield atom


def renumber_residues(atoms, residue_mapping=None, start_residue=1):
	res_count = start_residue-1
	for uniq_res_num, res_atoms in atoms_by_res(atoms):
		res_count += 1
		## first = True
		for atom in res_atoms:
			## if first:
			## 	print res_count, atom["res_name"]
			## 	first = False
			if residue_mapping is not None:
				residue_mapping.append(atom["res_num"], res_count)
			atom["res_num"] = res_count
			yield atom


def renumber_atoms(atoms):
	atom_count = 0
	for atom in atoms:
		atom_count += 1
		atom["atom_num"] = atom_count
		yield atom


def mutate_res(atoms, res_map, atom_maps, map_to_atom=False, delete_unmapped_atoms=False):
	for atom in atoms:
		modified = False
		res_name = atom["res_name"]
		if res_name in res_map:
			if res_name is None:
				# skip this residue type
				continue
			atom["res_name"] = res_map[res_name]
			modified = True
		if res_name in atom_maps:
			atom_map = atom_maps[res_name]
			atom_name = atom["atom_name"]
			if atom_name in atom_map:
				new_atom_name, new_element = atom_map[atom_name]
				if new_atom_name is None:
					# skip this atom
					continue
				else:
					atom["atom_name"] = new_atom_name
					atom["element"] = new_element
					modified = True
			elif delete_unmapped_atoms:
				# skip unmapped atom if requested
				continue

		if modified and map_to_atom:
			atom["rec_name"] = "ATOM  "

		yield atom


def residues_with_atoms(atoms, atom_names, pass_non_atoms):

	atom_names = set(atom_names)

	for uniq_res_num, res_atoms in atoms_by_res(atoms):
		collect_atoms = []

		res_atom_set = set()
		is_non_atom = False
		for atom in res_atoms:
			is_non_atom = (atom["rec_name"] != "ATOM  ")
			if pass_non_atoms:
				has_non_atom = has_non_atom or is_non_atom
			res_atom_set.add(atom["atom_name"])
			collect_atoms.append(atom)

		if atom_names.issubset(res_atom_set) or (pass_non_atoms and has_non_atom):
			for atom in collect_atoms:
				yield atom


def remove_alt_locs(atoms):
	ok_locs = set(" 1A")
	for atom in atoms:
		if atom["alt_loc"] in ok_locs:
			yield atom


MAX_REMARK_WIDTH = 69  # 80 minus "REMARK XXX "


def format_remark(s, number=99):
	# new remarks always start with a blank line
	yield ('REMARK %3d' % number).ljust(80)
	for i in xrange(0, len(s), MAX_REMARK_WIDTH):
		yield 'REMARK %3d %s' % (number, s[i:i + MAX_REMARK_WIDTH])


def aa3_seq(atoms):
	return (atom["res_name"] for atom in ca_atoms(atoms))


def to_fasta(atoms):
	return ''.join(nomenclature.aa3_seq_to_fasta(aa3_seq(atoms)))


# TODO : consider moving this to a separate file, e.g. pdbgeom.py
def apply_transform(atoms, matrix, shift_vec_before, shift_vec_after):
	for atom in atoms:
		atom["pos"] = apply_transform_to_vec(atom["pos"], matrix, shift_vec_before, shift_vec_after)
		yield atom


# matrix is expected to be a 3x3 jagged array
def apply_transform_to_vec(vec, matrix, shift_vec_before, shift_vec_after):
	dim = len(vec)
	ans = [0.0]*dim
	for i in range(dim):
		for j in range(dim):
			ans[i] = ans[i] + matrix[i][j] * (vec[j] + shift_vec_before[j])
		ans[i] = ans[i] + shift_vec_after[i]

	return ans


def clean(atoms,
		convert_mse=True,
		remove_pca_ace=False,
		chain_ranges=None,
		ranges_sequential_numbering=False,
		include_hetatms=False,
		do_renumber_atoms=False,
		do_renumber_residues=False,
		ensure_backbone_exists=True,
		remove_zero_occupancy=False,
		residue_mapping=None):

	if chain_ranges:
		atoms = ranges(atoms, chain_ranges, ranges_sequential_numbering)

	if remove_alt_locs:
		atoms = (atom for atom in atoms if atom["alt_loc"] in set(" 1A"))

	# altering residue/atom names	
	res_map = dict()
	atom_maps = dict()
	if remove_pca_ace:
		for aa in ('PCA', 'ACE'):
			res_map[aa] = None

	if convert_mse:
		res_map['MSE'] = 'MET'
		atom_maps['MSE'] = {'SE  ': (' SD ', 'S')}

	if (res_map or atom_maps):
		atoms = mutate_res(atoms, res_map, atom_maps, map_to_atom=True)

	if (not include_hetatms):
		atoms = (atom for atom in atoms if atom["rec_name"] != "HETATM")

	# remove atoms with zero occupancy
	if remove_zero_occupancy:
		atoms = (atom for atom in atoms if atom["occupancy"]>0)

	# remove residues without any backbone atom (N, Ca, C or O)
	if ensure_backbone_exists:
		atoms = residues_with_atoms(atoms, (" CA ", " N  ", " C  ", " O  "), pass_non_atoms=include_hetatms)

	# atom/residue renumbering
	if do_renumber_atoms:
		atoms = renumber_atoms(atoms)
	if do_renumber_residues:
		atoms = renumber_residues(atoms, residue_mapping)

	return atoms

def get_fasta_subcommand(args):
	import textwrap
	for pdb_file in args.pdb:
		print ">%s" % pdb_file.name
		for line in textwrap.wrap( to_fasta(ca_atoms(pdb_atoms(pdb_file))), 60 ):
			print line

def shift_subcommand(args):
	for atom in pdb_atoms(args.pdb):
		atom.x += args.x
		atom.y += args.y
		atom.z += args.z
		atom.res_num += args.res_num
		atom.chain = args.chain_map.setdefault(atom.chain, atom.chain)
		print str(atom)

def _bool_parse(arg):
	return str(arg).lower() in set(("1", "true", "yes", "t", "y"))

def output_atoms(atoms, add_ter=False):
	# output
	if add_ter:
		for chain, chain_atoms in atoms_by_chain(atoms):
			for atom in chain_atoms:
				yield str(atom)
			yield 'TER'.ljust(80)
	else:
		for atom in atoms:
			yield str(atom)

def clean_subcommand(args):
	range_res = None
	if args.range is not None:
		first_decapped, second = args.range[1:].split('-', 1)
		range_res = (int(args.range[0] + first_decapped), int(second))

	chain_ranges = dict()  # by default, include nothing
	if args.all_chains:
		chain_ranges = None
	elif args.chains is not None:
		# same range (None, or whatever is specified) for all chains
		chain_ranges = {chain: range_res for chain in args.chains}

	atoms = pdb_atoms(args.input, include_hetatms=(args.convert_mse or args.include_hetatms))

	atoms = clean(atoms,
			convert_mse=args.convert_mse,
			remove_pca_ace=args.remove_pca_ace,
			include_hetatms=args.include_hetatms,
			chain_ranges=chain_ranges,
			do_renumber_atoms=args.renumber_atoms,
			do_renumber_residues=args.renumber_residues,
			ranges_sequential_numbering=args.sequential_numbering,
			ensure_backbone_exists=args.ensure_backbone_exists,
			remove_zero_occupancy=args.remove_zero_occupancy)

	for line in output_atoms(atoms, add_ter=args.add_ter):
		print >> args.output, line

def dict_parse(dict_str):
	if not dict_str:
		return dict()
	return dict([ mapping.split(':') for mapping in dict_str.split(',') ])

def main():
	import argparse

	parser = argparse.ArgumentParser()

	subparsers = parser.add_subparsers(help='Sub-command help')

	clean_parser = subparsers.add_parser('clean', help='Clean or renumber a PDB file')
	clean_parser.add_argument('-p', '--input', required=True, type=filetools.opt_compressed_filetype('r'), help='Input PDB')
	clean_parser.add_argument('-o', '--output', required=True, type=filetools.opt_compressed_filetype('w'), help='Output PDB')
	clean_parser.add_argument('-c', '--chains', type=str, help='Extract specific chains to extract')
	clean_parser.add_argument('-r', '--range', type=str, help='Extract a specific range (same range used for all chains)')
	clean_parser.add_argument('-a', '--all-chains', action='store_true',
						help='Extract all chains (use this to just extract all ATOM records from the file).')
	clean_parser.add_argument('-q', '--quiet', action='store_true',
						help='Quiet mode. Does nothing currently, used for compatibility with the Perl verison of this script.')
	clean_parser.add_argument('-n', '--remove-pca-ace', action='store_true', help='Remove PCA and ACE residues.')
	clean_parser.add_argument('-s', '--renumber-atoms', action='store_true', help='Renumber atoms sequentially.')
	clean_parser.add_argument('-e', '--renumber-residues', action='store_true', help='Renumber residues sequentially.')
	clean_parser.add_argument('-m', '--no-convert-mse', dest='convert_mse', action='store_false',
						help='Do not convert MSE HETATM records to MET ATOM records.')
	clean_parser.add_argument('-l', '--include-hetatms', dest='include_hetatms', action='store_true',
						help='Do not remove HETATM records.')
	clean_parser.add_argument('-t', '--add-ter', action='store_true',
						help='Add TER record between chains (and at end of file).')
	clean_parser.add_argument('-x', '--sequential-numbering', action='store_true',
						help='Specify ranges in sequential (as they appear in the order of ATOM records) numbering.')
	clean_parser.add_argument('-b', '--ensure-backbone-exists', type=_bool_parse, default=True, metavar='T/F',
						help='Only print residues with N, CA, C atoms present. True by default.')
	clean_parser.add_argument('-z', '--remove-zero-occupancy', type=_bool_parse, default=True, metavar='T/F',
						help='Only print residues with N, CA, C atoms present. True by default.')

	clean_parser.set_defaults(func=clean_subcommand)

	get_fasta_parser = subparsers.add_parser('get_fasta', help='Clean or renumber a PDB file')
	get_fasta_parser.add_argument('pdb', nargs='+', type=argparse.FileType('r'), help='Input PDBs')
	get_fasta_parser.set_defaults(func=get_fasta_subcommand)

	shift_parser = subparsers.add_parser('shift', help='Shift PDB coordinates')
	shift_parser.add_argument('--x', type=float, help='Shift X coordinates by this amount', default=0)
	shift_parser.add_argument('--y', type=float, help='Shift Y coordinates by this amount', default=0)
	shift_parser.add_argument('--z', type=float, help='Shift Z coordinates by this amount', default=0)
 	shift_parser.add_argument('--res-num', type=float, help='Shift residue numbers by this amount', default=0)
	shift_parser.add_argument('--chain-map', metavar='X:Y,...', type=dict_parse, help='Map chain letter, where X:Y means chain X will be renamed to Y', default='')
	shift_parser.add_argument('pdb', type=argparse.FileType('r'), help='Input PDB')
	shift_parser.set_defaults(func=shift_subcommand)

	args = parser.parse_args()

	args.func(args)


if __name__ == "__main__":
	main()

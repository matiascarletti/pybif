#!/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

import itertools
import collections
import math

import aligntools
import algorithm

def eucledian_distance_sq(vec1, vec2):
	return sum( (x1-x2)**2 for (x1, x2) in zip(vec1, vec2) )


class SpatialHash:

	def __init__(self, distance):
		self.map = collections.defaultdict(list)
		self.distance = float(distance)
		self.dims = None
		self._neighbor_rel_keys = None

	def key(self, pos):
		return tuple( (int(math.floor(float(x)/self.distance)) for x in pos) )

	def add(self, pos, value):
		key = self.key(pos)
		dims = len(key)
		if self.dims is None:
			self.dims = dims
			self._neighbor_rel_keys = list(itertools.product((-1, 0, 1), repeat=dims))
		elif self.dims != dims:
			raise ValueError('keys inserted into SpatialHash must all be of same dimentionality')

		self.map[key].append((pos, value))

	def remove(self, pos):
		key = self.key(pos)
		# note : get() doesn't create an entry
		if key not in self.map:
			raise ValueError('SpatialHash.remove(pos): pos key not populated')
		try:
			self.map[key].remove( next(( pair for pair in self.map[key] if pair[0] == pos )) )
		except StopIteration:
			raise ValueError('SpatialHash.remove(pos): pos not in map')

	def keys_around(self, pos):
		key = self.key(pos)
		# go +1 or -1 in each possible direction from center
		for rel_key in self._neighbor_rel_keys:
			yield tuple(( x+y for x,y in zip(key, rel_key) ))

	def entries(self):
		return itertools.chain(*self.map.values())

	def keys(self):
		# return all populated spots
		return self.map.keys()

	def neighbors(self, pos):
		# note : get() doesn't create an entry
		distance_sq = self.distance ** 2
		for other_pos, value in itertools.chain( *(self.map.get(key, []) for key in self.keys_around(pos)) ):
			other_distance_sq = eucledian_distance_sq(pos, other_pos)
			if other_distance_sq <= distance_sq:
				yield (other_pos, value, math.sqrt(other_distance_sq))


def get_spatial_mapping(vecs1, vecs2, cutoff=4, sequential=True):
	sphash = SpatialHash(cutoff)

	vecs2 = list(vecs2)
	vecs2_len = len(vecs2)

	vecs1_len = 0
	for vec1_index, vec1 in enumerate(vecs1):
		vecs1_len += 1
		sphash.add(vec1, vec1_index)

	forward_map = [None] * vecs1_len
	reverse_map = [None] * vecs2_len
	ambiguous = list()
	# unmapped2 = 0

	# first align those atoms in atoms1 matching a single atom in atoms2
	prev_vec1_index = None
	for vec2_index, vec2 in enumerate(vecs2):
		# get two first neighbors
		# used to do sorted(..., key=lambda (vec1, vec1_index, distance): distance), but we need to search through the entire neighbor list for an ungapped match
		neighbors = list(sphash.neighbors(vec2))

		match = None
		ungapped_match = None
		if len(neighbors) == 1:
			match = neighbors[0]
		elif prev_vec1_index is not None:
			ungapped_match = [ (neighbor_vec1, neighbor_vec1_index, distance) \
					for neighbor_vec1, neighbor_vec1_index, distance in neighbors \
					if neighbor_vec1_index == (prev_vec1_index+1) ]
			if ungapped_match:
				match = ungapped_match[0]

		# print >> sys.stderr, "(%d) has %d neighbors up to %.2f, picked %s: %s" % (vec2_index, len(neighbors), cutoff, str(match), str([ vec1_index for (vec1, vec1_index, distance) in neighbors ]) ) 

		if match:
			# there is only one neighbour or we can make an ungapped alignment
			vec1, vec1_index, distance = match

			existing_match = forward_map[vec1_index]

			# align these two residues if either:
			# - an ungapped alignment is possible
			# - the residue in p1 is not already aligned with any other residue in p2
			#   (False if two residues in p2 are each close to a single residue in p1)
			# - the residue in p2 previously aligned with p1 is farther than the resdiue were evaluating now
			if ungapped_match or (existing_match is None) or (existing_match[1] > distance):
				# print >> sys.stderr, "Match (%d, %d) with distance %.2f" % (vec1_index, vec2_index, distance)
				forward_map[vec1_index] = (vec2_index, distance)
				reverse_map[vec2_index] = vec1_index
				prev_vec1_index = vec1_index

			# else:
			# 	unmapped2 += 1
		else:
			# print >> sys.stderr, "(%d) ambiguous, skipping for now" % (vec2_index) 
			ambiguous.append((vec2_index, neighbors))


	# print >> sys.stderr, "Stage 2 - match ambiguous stuff"
	# then, resolve ambiguities by matching to the nearest neighbor
	# here, there is bias in the order of appearance: early listed atoms will get their closest neighbor first
	for vec2_index, neighbors in ambiguous:
		# mapped = False
		prev_vec1_index = None
		if vec2_index>0:
			prev_vec1_index = reverse_map[vec2_index-1]

		for vec1, vec1_index, distance in neighbors:
			match = forward_map[vec1_index]
			is_ungapped = ((prev_vec1_index is not None) and vec1_index==(prev_vec1_index+1))
			if is_ungapped or (match is None) or (match[1] > distance):
				# print >> sys.stderr, "Match (%d, %d) with distance %.2f" % (vec1_index, vec2_index, distance)

				forward_map[vec1_index] = (vec2_index, distance)
				reverse_map[vec2_index] = vec1_index
				# mapped = True
				# don't continue search if we succeeded in an ungapped alignment
				if is_ungapped:
					# print >> sys.stderr, "Ungapped match found. Moving on to next residue."
					break # out of for .. in neighbors
		# if not mapped:
		#	unmapped2 += 1

	# use 1-based sequential numbering
	mapping = [ (vec1_index+1, pair[0]+1) for vec1_index, pair in enumerate(forward_map) if pair is not None ]

	# the mapping above might produce non-sequential mappings such as
	#     1 <-> 9
    #     . . . 
    #     8 <-> 3
	# is requested, get a (longest) sequential subset of the mapping
	if sequential:
		mapping = algorithm.longest_increasing_subsequence(mapping, lt=lambda this, that: all(( i<j for i,j in zip(this, that) )) )

	return mapping

if __name__ == "__test__":
	import sys
	import random
	import timeit

	for p in xrange(40):
		cutoff = random.random()*0.1
		h = SpatialHash(cutoff)
		(x0, y0) = (random.random(), random.random())
		should_see = list()
		for i in xrange(4000):
			(x,y) = (random.random(), random.random())
			dist = math.sqrt((x-x0)**2 + (y-y0)**2)
			if dist <= cutoff:
				should_see.append(((x,y),i,dist))
			h.add((x,y), i)

	cutoff = 0.01
	for p in xrange(3,8):
		h = SpatialHash(cutoff)
		(x0, y0) = (random.random(), random.random())
		should_see = list()
		for i in xrange(10**p):
			(x,y) = (random.random(), random.random())
			dist = math.sqrt((x-x0)**2 + (y-y0)**2)
			if dist <= cutoff:
				should_see.append(((x,y),i,dist))
			h.add((x,y), i)

		print (p, timeit.timeit("found = list(h.neighbors((x0, y0)))", "from __main__ import h, x0, y0", number=10))
		found = list(h.neighbors((x0, y0)))

		assert(set(should_see)==set(found))

	print should_see
	print found

if __name__ == "__main__":
	import argparse

	import pdb

	parser = argparse.ArgumentParser(description='Gets an alignment of CAs given two superimposed structures')
	parser.add_argument('--cutoff', default=4, help='maximum distance to allow alignment for, in angstroms')
	parser.add_argument("pdb1", type=argparse.FileType('r'))
	parser.add_argument("pdb2", type=argparse.FileType('r'))
	args = parser.parse_args()

	cas1 = list(pdb.ca_atoms(pdb.pdb_atoms(args.pdb1)))
	cas2 = list(pdb.ca_atoms(pdb.pdb_atoms(args.pdb2)))
	vecs1 = [ ca.pos for ca in cas1 ]
	vecs2 = [ ca.pos for ca in cas2 ]
	mapping = get_spatial_mapping(vecs1, vecs2)
	alignment = aligntools.mapping_to_alignment(mapping, len(cas1), len(cas2))
	aligned_fasta1, aligned_fasta2 = aligntools.alignment_fasta(alignment, pdb.to_fasta(cas1), pdb.to_fasta(cas2))
	print ">%s" % args.pdb1.name
	print aligned_fasta1
	print ">%s" % args.pdb2.name
	print aligned_fasta2




def dummy():
	neighbor_map = [ list(sphash.neighbors(vec2)) for vec2 in vecs2 ]

	for vec2_index, vec2 in enumerate(vecs2):
		neighbors = neighbor_map[vec2_index]
		if len(neighbors) == 1:
			vec1, vec1_index, distance = match
			forward_map[vec1_index] = vec2_index
			reverse_map[vec2_index] = vec1_index
			
			prev_vec1_index = vec1_index
			if vec2_index > 0:
				for other_vec2_index, other_vec2 in reversed(vecs2[:vec2_index-1]):
					ungapped_match = [ (other_vec1, other_vec1_index, distance) \
							for other_vec1, other_vec1_index, distance in neighbors \
							if other_vec1_index == (prev_vec1_index-1) ]
					if not ungapped_match:
						break
					forward_map[other_vec1_index] = other_vec2_index
					reverse_map[other_vec2_index] = other_vec1_index
					prev_vec1_index = other_vec1_index

			prev_vec1_index = vec1_index
			if vec2_index < vecs2_len-1:
				for other_vec2_index, other_vec2 in vecs2[vec2_index+1]:
					ungapped_match = [ (other_vec1, other_vec1_index, distance) \
							for other_vec1, other_vec1_index, distance in neighbors \
							if other_vec1_index == (prev_vec1_index+1) ]
					if not ungapped_match:
						break
					forward_map[other_vec1_index] = other_vec2_index
					reverse_map[other_vec2_index] = other_vec1_index
					prev_vec1_index = other_vec1_index




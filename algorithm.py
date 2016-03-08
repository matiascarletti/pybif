#!/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

import math
import operator

def longest_increasing_subsequence(seq, key=None, lt=None):
	"""Returns the longest increasing subsequence in the given sequence"""
	seq_len = len(seq)

	subseq_pred = [0] * (seq_len)
	subseq_last = [0] * (seq_len+1)

	if key is not None:
		seq = [ (key(val), val) for val in seq ]

	subseq_len = 0

	if lt == None:
		lt = operator.lt

	for i in range(seq_len):
		lo = 1
		hi = subseq_len

		while lo <= hi:
			mid = int(math.ceil((lo+hi)/2))

			if (lt(seq[subseq_last[mid]], seq[i])):
				lo = mid+1
			else:
				hi = mid-1

 
		new_subseq_len = lo
		subseq_pred[i] = subseq_last[new_subseq_len-1]
		subseq_last[new_subseq_len] = i
 
		if (new_subseq_len > subseq_len):
			subseq_len = new_subseq_len
 
	result = []
	k = subseq_last[subseq_len]
	for i in range(subseq_len-1, -1, -1):
		result.append(seq[k])
		k = subseq_pred[k]

	if key is not None:
		result = [ v for k, v in result ]

	return result[::-1]
 
if __name__ == '__main__':
	for d in [[3,2,6,4,5,1], [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15]]:
		print('a L.I.S. of %s is %s' % (d, longest_increasing_subsequence(d)))

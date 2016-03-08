#!/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

import os
import gzip
import sys

# adapted from FileType class in argparse module
def open_by_ext_filetype(ext_opener_map, mode):
	def g(filename):
		if filename == '-':
			if 'r' in mode:
				return sys.stdin
			elif 'w' in mode:
				return sys.stdout
			else:
				raise ValueError('Cannot open file "-" with mode %s' % mode)

		ext = os.path.splitext(filename)[1]
		opener = ext_opener_map.get(ext, open)
		return opener(filename, mode)

	return g

def opt_compressed_filetype(mode):
	return open_by_ext_filetype({'.gz': gzip.open}, mode)


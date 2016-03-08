#!/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

import subprocess
import os

_stride = None

def get_stride():
	global _stride
	if ( _stride == None ):
		_stride = os.getenv('STRIDE_EXEC', 'stride')
	return _stride

def stride(pdb_file):
	process = subprocess.Popen([get_stride(), '/dev/stdin', '-o'], stdin=pdb_file, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out_data, err_data = process.communicate()

	seq = list()
	secstr = list()
	file = None
	chain = None

	for line in out_data.splitlines():
		if line.startswith('CHN  '):
			_, file, chain, _ = line.split(' ', 3)
			new_result = True
		elif line.startswith('SEQ  '):
			seq.append(line[10:60])
		elif line.startswith('STR  '):
			secstr.append(line[10:60])
		elif line.startswith('LOC  ') and new_result:
			seq = (''.join(seq)).strip()
			secstr = (''.join(secstr))[:len(seq)]
			yield (file, chain, seq, secstr)
			seq = list()
			secstr = list()
			file = None
			chain = None
			new_result = False

def secstr(pdb_file):
	file, chain, seq, secstr = stride(pdb_file).next()
	return secstr

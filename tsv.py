#!/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

from collections import OrderedDict

def itemmapper(**kwargs):
	if 'kwargs' in kwargs: # workaround for unordering that happens when unpacking an OrderedDict kwargs: https://www.python.org/dev/peps/pep-0468/
		kwargs = kwargs['kwargs']
	def g(obj):
		return OrderedDict(( (key, obj[kwargs[key]]) for key in kwargs))
	return g

def tsv_rows(tsv_file, mapper=None):
	for line in tsv_file:
		if line.startswith('#'):
			if mapper is None:
				mapping = OrderedDict(( (field, index) for index, field in enumerate(line[1:].strip().split('\t')) ))
				mapper = lambda obj: OrderedDict(( (field, obj[mapping[field]]) for field in mapping ))
				# itemmapper(**OrderedDict( ( (field_name, index) for index, field_name in enumerate(line[1:].strip().split('\t')) ) )) 
			continue
		if mapper is None:
			raise ValueError('A header line must appear before data in the file %s' % tsv_file.name)
		yield mapper(line.strip().split('\t'))



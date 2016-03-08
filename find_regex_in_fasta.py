#/usr/bin/env python

# PyBIF library. Copyright (c) 2016 The Furman Lab, HUJI
# Use of this program under the terms of the license specified in the LICENSE
# file, which should be included along with this program.

import re
import argparse

import argparsetools
import fasta

def compile_raw_regex(regex):
	return re.compile('(?=(' + regex + '))')

def parse_elm_motifs(lines):
	motifs = []
	first_line = True
	for line in lines:
		if line.startswith('#'):
			continue
		if first_line:
			# Remove the first elements from the lists, as they are column names
			first_line = False
			continue
		parts = line.replace('"', '').strip('\t').split('\t')
		elm_id      = parts[1]
		motif_regex = compile_raw_regex(parts[3])
		# we don't use this, but if someone else uses this, they might be interested in this:
		# evalue = float(parts[4])
		motifs.append((elm_id, motif_regex))

	return motifs


def parse_regex_arg(regex_arg):
	parts = regex_arg.split(':', 1)
	if len(parts) > 1:
		name, raw_regex = parts
	else:
		name, raw_regex = regex_arg, regex_arg
	return (name, compile_raw_regex(raw_regex))


def parse_acc_filter(seqs, parse):
	for name, seq in seqs:
		if len(name) > 3 and name[3] == '|': # if the header is in the form >xx|<accession>
			_, name, _ = name.split('|', 2)

		yield name, seq


def main():

	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta', required=True, type=argparsetools.optionally_compressed_filetype('r'), help='FASTA file to search (may be gzipped)')
	parser.add_argument('--parse-acc', action='store_true', help='Parse the accession in the sequence identifier, if it is in the format >xx|<accession>|..., such as in Uniprot-provided FASTA files')
	parser.add_argument('--regex', nargs='*', metavar='PATTERN', help='Regular expression pattern(s). PATTERN may be preceded by the name of the pattern followed by a colon, e.g. MY_LIGAND:PP.Y')
	parser.add_argument('--regex-file', help='File containing list of regular expression patterns. See --regex-file-format')
	parser.add_argument('--regex-file-format', choices=['args', 'elm'], help='Format of the input file. Specify "elm" if the input is in ELM elm_classes.tsv format. Specify "args" if the file contains one line for each regular expression, in the format expected by the --regex argument (lines beginning with # will be ignored). Default is "args".', default='args')
	args = parser.parse_args()

	motifs = []
	if args.regex_file is not None:
		if args.regex_file_format == "args":
			for line in args.regex_file:
				if line.startswith('#'):
					continue
				motifs.append( parse_regex_arg(line) )
		else:
			motifs = parse_elm_motifs(args.regex_file)
	
	if args.regex is not None:
		for regex_arg in args.regex:
			motifs.append( parse_regex_arg(regex_arg) )

	if not motifs:
		raise RuntimeError('No motifs were specified. Either specify --regex or --regex-file.')
	
	seqs = fasta.fasta_seqs(args.fasta, strip_leading_bracket=True)
	
	if args.parse_acc:
		seqs = parse_acc_filter(seqs)
	
	for name, seq in seqs:
		for id, motif_regex in motifs:
			for match in motif_regex.finditer(seq):
				if match:
					print '\t'.join(map(str, [name, id, match.start(1)+1, (match.end(1)-match.start(1))]))


if __name__ == '__main__':
	main()

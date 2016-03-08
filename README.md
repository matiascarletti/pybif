# PyBIF

## Overview

This is some initial work on a Python library for Bioinformatics affectionately dubbed PyBIF (Python bioinformatics from the Furman lab).

It was built ad-hoc for construction of the homoogue loop-modeling benchmark, and will hopefully build on from there. As things were built gradually, some duplicacy still exists.

The focus was both on lightness (modules should load as fast as possible) and modularity (keep functions and modules focused on very narrow tasks to minimize interdependece).

Existing libraries, such as @BioPython and @klab, focus on the same problems, and the hope in writing this was to have a fresh take on the same ideas and find new solutions. Ultimately, it'll be good to incorporate any useful ideas presented in this library to the other existing libraries.

### Guidelines for lightness

- modules should not initialize constants upon loading, but upon first use
- use Python slots data model to minimize memory footprint for classes that get instansiated a lot (e.g. PDB atom entries)
- concise naming

### Guidelines for modularity

- keep functions and modules focused on very narrow tasks, trying to stick to very abstract, atomic notions of information processing
- avoid dedicated classes and resort to Python builtins as much as possible, to allow interoperability with other libraries

## Features

- methods for parsing
  - tab-separated tabular data files (`tsv`)
  - PDB files (`pdb`)
  - FASTA files (`fasta`)
  - the ECOD database (`ecod`)
- methods for manipulating and querying PDB files (`pdb`)
- algorithms
  - longest increasing subsequence (`algorithm`)
  - naive spatial alignment (`spatialalign`)
- wrappers for common bioinformatic executables
  - STRIDE
  - MAMMOTH
  - jFATCAT (BioJava implementation of the FATCAT algorithm)
- several command line utilities
  - searching a regex over a fasta file (`find_regex_in_fasta.py`)
  - getting FASTA sequence from a PDB (`pdb.py get_fasta`)
  - cleaning PDB files to include only coordinates (to make them ready for Rosetta, `pdb.py clean`)
  - editing PDB files (`pdb.py shift`)

## Ideas for the future

- separate command-line functionality from library modules
- separate manipulation from parsing in `pdb` library
- provide residue-view in the `pdb` module (`dict` by `atom_name`)
- read other types of records in the `pdb` module
- extend `tsv` to handle more tabular file formats

## License

See `LICENSE` for license information.

## Contact

The code is currently being maintained by @kwikwag
The Furman Lab website is at http://www.cs.huji.ac.il/~fora/

﻿# PDB-Decoder-project-using-Biopython

This program decodes information from different PDB (Protein Bank Data) files and write them down in CSV file.

Information that can be read from a file:

atom serial number,
atom name (spaces stripped, e.g. 'CA'),
residue name, 
chain id (e.g. 'VAL','GLY'),
residue number, 
coordinates x,y,z,
b factor, 
occupancy,
cheking if the atom is a C-alpha atom (if there's return 1, and otherwise 0)



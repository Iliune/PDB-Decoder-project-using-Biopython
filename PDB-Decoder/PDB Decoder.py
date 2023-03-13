from Bio.PDB import PDBParser

#Define the PDB file path and the output CSV file path
pdb_file = "1phm.pdb"
csv_file = "output.csv"

#Create a PDB parser object
parser = PDBParser()

#Parse the PDB file and get the structure object
structure = parser.get_structure("structure", pdb_file)

#Open the CSV file for writing
with open(csv_file, "w") as f:
    
    #Write the header row to the CSV file
    f.write("atom_serial, atom_symbol, CA, residue_name, chain_id, residue_number ,x_coord,y_coord,z_coord, b_factor, occupancy\n")
    
    #Iterate over all atoms in the structure
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    
                    #Extract the atom information
                    atom_serial = atom.get_serial_number()
                    atom_name = atom.get_name()
                    residue_name = residue.get_resname()
                    chain_id = chain.get_id()
                    residue_number = residue.get_id()[1]
                    x_coord, y_coord, z_coord = atom.get_coord()
                    b_factor = atom.get_bfactor()
                    occupancy = atom.get_occupancy()
                    CA = int(atom_name == "CA")
                    
                    # Write the atom information to the CSV file
                    f.write(f"{atom_serial},{atom_name},{CA},{residue_name},{chain_id},{residue_number},{x_coord:.3f},{y_coord:.3f},{z_coord:.3f},{b_factor},{occupancy}\n")

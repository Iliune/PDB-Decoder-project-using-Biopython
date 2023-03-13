'''
It's the same program, but rewrite using classes
'''

from Bio.PDB import PDBParser

class PDBFile:
    def __init__(self, pdb_file_path, csv_file_path):
        self.pdb_file_path = pdb_file_path
        self.csv_file_path = csv_file_path
    
    def parse_and_save_to_csv(self):
        #Create a PDB parser object
        parser = PDBParser()

        #Parse the PDB file and get the structure object
        structure = parser.get_structure("structure", self.pdb_file_path)

        #Open the CSV file for writing
        with open(self.csv_file_path, "w") as f:
            
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
                            
                            #Write the atom information to the CSV file
                            f.write(f"{atom_serial},{atom_name},{CA},{residue_name},{chain_id},{residue_number},{x_coord:.3f},{y_coord:.3f},{z_coord:.3f},{b_factor},{occupancy}\n")

'''
you can upload your own data, here's example of mine
'''

pdb_file = "1phm.pdb"
csv_file = "output.csv"
pdb = PDBFile(pdb_file, csv_file)
pdb.parse_and_save_to_csv()
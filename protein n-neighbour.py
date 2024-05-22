
import numpy as np
from Bio.PDB import PDBParser
from Bio import PDB

def convert_cif_to_pdb(cif_file_path, pdb_file_path):
    parser = PDB.MMCIFParser()
    structure = parser.get_structure('protein', cif_file_path)

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(pdb_file_path)

def extract_pdb_info(pdb_file_path):
    """
    Reads a PDB file and extracts 3D spatial coordinates and amino acids.

    Parameters:
    pdb_file_path (str): Path to the PDB file.

    Returns:
    list: A list of dictionaries, each containing information about an atom.
    """
    # Create a PDBParser object
    parser = PDBParser(PERMISSIVE=1)
    
    # Parse the structure from the PDB file
    structure = parser.get_structure('protein', pdb_file_path)
    
    # List to hold the extracted information
    atom_info_list = []

    # Extract information from the structure
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_info = {
                        'model_id': model.id,
                        'chain_id': chain.id,
                        'residue_name': residue.resname,
                        'residue_id': residue.id[1],
                        'atom_name': atom.name,
                        'atom_coords': atom.coord.tolist()
                    }
                    atom_info_list.append(atom_info)

    return atom_info_list


def count_neighboring_atoms(data):

    def distance_np(coord1, coord2):
        return np.linalg.norm(coord1 - coord2)

    def atom_nearest_neighbours(data):

        coords = np.array([atom['atom_coords'] for atom in data])
        neighbour_counts = np.zeros(len(data), dtype=int)
        
        for i in range(len(coords)):
            for j in range(i+1, len(coords)):
                distances = distance_np(coords[i],coords[j])
                if distances <= 3:
                    neighbour_counts[i] += 1

        for i in range(len(neighbour_counts)):
            data[i]['neighbour_count'] = neighbour_counts[i]
        print(neighbour_counts)
        
        return data

    neighbour_counts = atom_nearest_neighbours(data)

    '''
    #delete irrelevant ids in data
    for entry in range(len(neighbour_counts)):
        del neighbour_counts[entry]['model_id']
        del neighbour_counts[entry]['chain_id']
        # del neighbour_counts[entry]['residue_name']
        del neighbour_counts[entry]['residue_id']
    '''
        
    return neighbour_counts

'''
def count_neighbouring_amino_acids(data):

    def split_into_acids(data):
        amino_acid = []
        amino_acid_list = []

        for i in range(len(data)):
            acid_name = data[i]['residue_name']

            if i+1 < len(data) and data[i+1]['residue_name'] == acid_name:
                amino_acid.append(data[i])
            else:
                amino_acid.append(data[i])
                amino_acid_list.append(amino_acid)
                amino_acid = []

        return amino_acid_list
    
    def find_highest_count(amino_acid_list):

        for i in range(len(amino_acid_list)):
            highest_value = 0
            filtered_list = []

            for j in range(len(amino_acid_list[i])):
                

                counts = [for amino_acid_list[i][j]['neighbour_count'] in amino_acid_list[i]]
                highest_value = max(counts)
                print(highest_value)

                if amino_acid_list[i][j]['neighbour_count'] == highest_value:
                    filtered_list.append(amino_acid_list[i][j])

        return filtered_list

    split = split_into_acids(data)
    filtered_acids = find_highest_count(split)
    print(filtered_acids)

    return(filtered_acids)
'''

def count_neighbouring_amino_acids(data):
    # Dictionary to store the highest 'neighbour_count' for each 'residue_id'
    max_neighbour_counts = {}

    # Iterate through the data to find the highest 'neighbour_count' for each 'residue_id'
    for atom in data:
        residue_id = atom['residue_id']
        neighbour_count = atom.get('neighbour_count', 0)
        if residue_id not in max_neighbour_counts or neighbour_count > max_neighbour_counts[residue_id]:
            max_neighbour_counts[residue_id] = neighbour_count

    # Filter the original data to retain only one entry for each 'residue_id' with the highest 'neighbour_count'
    filtered_data = []
    for atom in data:
        residue_id = atom['residue_id']
        neighbour_count = atom.get('neighbour_count', 0)
        if neighbour_count == max_neighbour_counts[residue_id]:
            # Remove other entries for the same residue_id
            filtered_data = [entry for entry in filtered_data if entry['residue_id'] != residue_id]
            filtered_data.append(atom)
            atom.pop('model_id')
            atom.pop('atom_name')
            atom.pop('atom_coords')

    return filtered_data

# ============== MAIN =======================

cif_file_path = '/Users/reishuenng/Downloads/1vyb-assembly1.cif'
pdb_file_path = '/Users/reishuenng/Downloads/1vyb-assembly1.pdb'

'''
cif_file_path = '/Users/reishuenng/Downloads/7sfc_updated.cif'
pdb_file_path = '/Users/reishuenng/Downloads/pdb7sfc.pdb'
'''

convert_cif_to_pdb(cif_file_path, pdb_file_path)
atom_info_list = extract_pdb_info(pdb_file_path)
print(len(atom_info_list))
atom_info_count = count_neighboring_atoms(atom_info_list[:])

#Check loading
for atom_info in atom_info_count[:20]:  
    print(atom_info)

amino_acid_info_count = count_neighbouring_amino_acids(atom_info_count)

#Check amino_acid_neighbour, should have 1389
for amino_acid in amino_acid_info_count:  
    print(f"{amino_acid}")

print(f"Number of Amino Acids: {len(amino_acid_info_count)}")






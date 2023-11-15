import numpy as np
from Bio import PDB
import matplotlib.pyplot as plt
from Bio.PDB import Polypeptide
from Bio.PDB.DSSP import DSSP

def calculate_phi_psi(structure, pdb_filename):
    phi_angles = []
    psi_angles = []
    secondary_structure = []

    for model in structure:
        for chain in model:
            polypeptides = Polypeptide.Polypeptide(chain)
            phi_psi_angles = polypeptides.get_phi_psi_list()
            dssp = PDB.DSSP(model, pdb_filename)

            residues = list(chain.get_residues())
            for phi, psi in phi_psi_angles:
                if phi is not None and psi is not None:
                    residue_index = polypeptides.get_phi_psi_list().index((phi, psi))
                    residue_id = (chain.id, residues[residue_index].get_id()[1])
                    phi_angles.append(np.degrees(phi))
                    psi_angles.append(np.degrees(psi))
                    secondary_structure.append(dssp[residue_id][2])

    return phi_angles, psi_angles, secondary_structure

def plot_ramachandran(phi_angles, psi_angles, secondary_structure):
    plt.figure(figsize=(8, 6))
    
    # Map secondary structure types to colors
    color_map = {'H': 'red', 'E': 'blue', '-': 'green', 'G': 'orange', 'T': 'purple', 'S': 'yellow', 'I': 'cyan', 'B': 'pink'}

    # Scatter points with color based on secondary structure
    for phi, psi, ss in zip(phi_angles, psi_angles, secondary_structure):
        plt.scatter(phi, psi, s=10, alpha=0.5, edgecolors='none', c=color_map[ss])

    plt.title('Ramachandran Plot with Secondary Structure Coloring')
    plt.xlabel('Phi Angle (degrees)')
    plt.ylabel('Psi Angle (degrees)')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

    # Create a legend
    legend_labels = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label)
                     for label, color in color_map.items()]
    plt.legend(handles=legend_labels, title='Secondary Structure', loc='upper right')
    plt.savefig('ramachandran_plot.pdf')
    plt.show()

if __name__ == "__main__":
    pdb_id = "1HBB"
    pdb_filename = f"{pdb_id.lower()}.pdb"

    parser = PDB.PDBParser(QUIET=True)

    structure = parser.get_structure(pdb_id, pdb_filename)

    phi_angles, psi_angles, secondary_structure = calculate_phi_psi(structure, pdb_filename)

    plot_ramachandran(phi_angles, psi_angles, secondary_structure)

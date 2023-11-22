import sys
import numpy as np
from Bio import PDB     # do obliczania kątów phi i psi
import matplotlib.pyplot as plt  # do rysowania Ramachandran Plot
from Bio.PDB import Polypeptide # stąd  informacje o kątach torsyjnych (z obiektu Polypeptide)
from Bio.PDB.DSSP import DSSP  # do uzyskania informacji o strukturze wtórnej.
import os

def calculate_phi_psi(structure):
    # Inicjalizacja list do przechowywania kątów phi, psi i struktury wtórnej
    phi_angles = []
    psi_angles = []
    secondary_structure = []

    for model in structure:      # Iteracja po modelach w strukturze białka
        for chain in model:      # Iteracja po łańcuchach w każdym modelu
            polypeptides = Polypeptide.Polypeptide(chain)    # Utworzenie obiektu Polypeptide dla danego łańcucha -  klasy, która umożliwia analizę sekwencji aminokwasów w ramach łańcucha białkowego.
            phi_psi_angles = polypeptides.get_phi_psi_list() # Pobranie listy kątów phi i psi dla danego łańcucha
            dssp = PDB.DSSP(model, pdb_filename)  #Utworzenie obiektu DSSP dla danego modelu do analizy struktury wtórnej (jak helisy alfa, arkusze beta, zgięcia itp.)

            residues = list(chain.get_residues()) # uzyskujemy pełną listę reszt aminokwasowych w danym łańcuchu
            for phi, psi in phi_psi_angles:       #  Iteracja po parach kątów phi i psi
                if phi is not None and psi is not None: # Sprawdzenie, czy oba kąty są zdefiniowane
                    residue_index = polypeptides.get_phi_psi_list().index((phi, psi)) # Indeks reszty aminokwasowej w Polypeptide
                    residue_id = (chain.id, residues[residue_index].get_id()[1]) #tworzenia identyfikatora reszty aminokwasowej, który składa się z identyfikatora łańcucha oraz numeru reszty
                    phi_angles.append(np.degrees(phi))              # Dodanie kątów phi i psi w stopniach do list
                    psi_angles.append(np.degrees(psi))
                    secondary_structure.append(dssp[residue_id][2])  # Dodanie informacji o strukturze wtórnej do listy

    return phi_angles, psi_angles, secondary_structure  # Zwrócenie list kątów phi, psi i struktury wtórnej

def plot_ramachandran(phi_angles, psi_angles, secondary_structure, pdb_filename):
    plt.figure(figsize=(8, 6))  # Utworzenie nowego obiektu figury o rozmiarze 8x6 cali
    
    # Mapowanie typów struktury wtórnej na kolory
    color_map = {'H': 'red', 'E': 'blue', '-': 'green', 'G': 'orange', 'T': 'purple', 'S': 'yellow', 'I': 'cyan', 'B': 'pink'}

    # Rozrzucanie punktów na wykresie, kolorując je w zależności od struktury wtórnej
    for phi, psi, ss in zip(phi_angles, psi_angles, secondary_structure):
        plt.scatter(phi, psi, s=10, alpha=0.5, edgecolors='none', c=color_map[ss]) #Dodaje punkt na wykresie.

    plt.title('Ramachandran Plot with Secondary Structure Coloring')  # Dodanie tytułu i opisów osi
    plt.xlabel('Phi Angle (degrees)')
    plt.ylabel('Psi Angle (degrees)')
    plt.xlim(-180, 180)                                               # Ustawienie zakresów osi
    plt.ylim(-180, 180)
    plt.axhline(0, color='black', linewidth=0.5)                       # Dodanie linii odniesienia w poziomie i pionie
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.5)   # Dodanie siatki pomocniczej


    legend_labels = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label) #Jest to lista obiektów Line2D, reprezentujących punkty w legendzie.
                     for label, color in color_map.items()]
    plt.legend(handles=legend_labels, title='Secondary Structure', loc='upper right') #Tworzy legendę na podstawie utworzonych wcześniej obiektów Line2D.
    plt.savefig(pdb_filename+'_ramachandran_plot.pdf')  #zapisanie
    plt.show()                            #wyświetlenie

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <pdb_file_path>")
        sys.exit(1)

    pdb_filename = sys.argv[1]

    parser = PDB.PDBParser(QUIET=True)  # utworzenie instancji parsera PDB z biblioteki BioPython

    try:
        structure = parser.get_structure("structure", pdb_filename) # uzyskania obiektu reprezentującego strukturę białka z pliku PDB 
    except FileNotFoundError:
        print(f"File {pdb_filename} not found.")
        sys.exit(1)

    phi_angles, psi_angles, secondary_structure = calculate_phi_psi(structure)
    #phi_angles: Lista zawierająca kąty phi (ϕ) w stopniach dla poszczególnych reszt aminokwasowych w strukturze białka.

    #psi_angles: Lista zawierająca kąty psi (ψ) w stopniach dla poszczególnych reszt aminokwasowych w strukturze białka.

    #secondary_structure: Lista zawierająca informacje o strukturze wtórnej dla poszczególnych reszt aminokwasowych w formie oznaczeń typów struktury wtórnej (np. helisa, arkusz beta, zgięcie itp.).


    filename, file_extension = os.path.splitext(os.path.basename(pdb_filename))
    plot_ramachandran(phi_angles, psi_angles, secondary_structure, os.path.basename(filename+"_"+file_extension[1:]))

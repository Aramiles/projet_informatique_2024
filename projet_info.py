import pandas as pd
import os
import streamlit as st
os.chdir(os.path.dirname(os.path.abspath(__file__)))
organic_atoms = ['H', 'B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']

#to reduce the dimension of a list

def flat_list(list):
    flat = []
    [flat.extend(mini_list) for mini_list in list]
    return flat

#class for each molecule

class molecule:
    def __init__(self, fichier_name):
        self.fichier = open(os.path.join(os.path.dirname(os.path.abspath(__file__))+"\MOL_files", fichier_name),"r").readlines()         #read the file 
        self.fichier_name = fichier_name

        #overall informations found in the head of the .mol file
        self.get_infos = {'brut formula':self.fichier[0][:-1], 'nb of atoms':int(self.fichier[3][:-1].split()[0]), 'nb of bonds':int(self.fichier[3][:-1].split()[1])}
        list_fichier = [k[:-1].split() for k in self.fichier[4:]]
        
        #atom bloc where we can find informations about atoms in a dataframe 
        self.get_atoms = pd.DataFrame([k for k in list_fichier if len(k)==len(list_fichier[0])], #atom bloc begins at line 4. We take all the line of the same lenght in the file
                                      index = [k+1 for k in range (len([k for k in list_fichier if len(k)==len(list_fichier[0])]))], #index corresponds to atom number, it begins at 1
                                      columns = ['x', 'y', 'z','atom', 'mass difference', 'charge', 'atom stereo parity']+
                                      [k for k in range (len([k[:-1].split() for k in self.fichier[4:] if len(k)==len(self.fichier[4])][0])-7)])       #we don't know the title of the other columns

        # bond bloc where we can find informations about bonds in a dataframe 
        self.get_bonds = pd.DataFrame([k[:-1].split() for k in self.fichier[4:-2] if len(k[:-1].split())==7],        #As there is 7 features to describe a bond, we take the lines of lenght 7 of the file
                                      columns = ['first atom nb','second atom nb', 'order', 'stereo','useless','topology', 'reacting center status'])

        

    def load_verify(self):
        if all(atom in organic_atoms for atom in self.get_atoms['atom']) is False:           #we verify that all the atoms in the atom_bloc are in the list of the organic atoms
            print(f"File {self.fichier_name} could not be loaded. Contains non-organic atom elements.")
            return False
        elif len(self.get_bonds) != int(self.get_infos['nb of bonds']):                #we  verify that lenght of the bond bloc corresponds to the number given in the head bloc
            print(f'File {self.fichier_name} could not be loaded. Bonds block is incomplete.')
            return False
        elif len(set(list(self.get_bonds['first atom nb'])+list(self.get_bonds['second atom nb']))) != len(self.get_atoms): #we  verify that all the atoms are named somewhere in the bond bloc
            print(f'File {self.fichier_name} could not be loaded. Bonds block does not cover all atoms.')
            return False
        else:
            print(f"File {self.fichier_name} validated and loaded!")
            return True



    def load_verify_st(self):
        if all(atom in organic_atoms for atom in self.get_atoms['atom']) is False:           #we verify that all the atoms in the atom_bloc are in the list of the organic atoms
            st.write(f"File {self.fichier_name} could not be loaded. Contains non-organic atom elements.")
            return False
        elif len(self.get_bonds) != int(self.get_infos['nb of bonds']):                #we  verify that lenght of the bond bloc corresponds to the number given in the head bloc
            st.write(f'File {self.fichier_name} could not be loaded. Bonds block is incomplete.')
            return False
        elif len(set(list(self.get_bonds['first atom nb'])+list(self.get_bonds['second atom nb']))) != len(self.get_atoms): #we  verify that all the atoms are named somewhere in the bond bloc
            st.write(f'File {self.fichier_name} could not be loaded. Bonds block does not cover all atoms.')
            return False
        else:
            st.write(f"File {self.fichier_name} validated and loaded!")
            return True

    def count_elements(self, element):
        return len([k for k in list(self.get_atoms['atom']) if k == element])   #we add a string in the list each time the letter of the element desired is named. the lenght give us the number of atoms of the given element. 


    def distance_3D(self, atom1, atom2):
        return round(((float(self.get_atoms['x'][atom1])-float(self.get_atoms['x'][atom2]))**2 + (float(self.get_atoms['y'][atom1])-float(self.get_atoms['y'][atom2]))**2 + (float(self.get_atoms['z'][atom1])-float(self.get_atoms['z'][atom2]))**2)**(1/2),2)                #calcul of the distance
        

    def atom_neighbours(self, atom_identifier):
        return [int(self.get_bonds['second atom nb'][i]) for i, atoms in enumerate(self.get_bonds['first atom nb']) if int(atoms) == atom_identifier]+[int(self.get_bonds['first atom nb'][i]) for i, atoms in enumerate(self.get_bonds['second atom nb']) if int(atoms) == atom_identifier]      #we put in a list the second atom's identifier each time the targeted atom is seen in the 'atom' column of the atom bloc. We'll return this list because it can be useful for the topological distance between the atoms 


    def distance_topo_2D(self, atom1, atom2):
        neighbours_path = self.atom_neighbours(atom1)   #As mentioned, we use the method above to find the first circle of atoms around the atom1
        k = 1
        while atom2 not in neighbours_path:   #We increase the diameter of the circle until the circle meet the atom2. Meaning that, for each atom in the neighbours_path's list, we put there neighbours in the list, we increase the k counter and we check that atom2 is not in the list. We repeat this until he is.
            neighbours_path += flat_list([self.atom_neighbours(neighbour) for neighbour in neighbours_path])
            k += 1
        return k

    def ring_finding(self, ring_path):
        graph = {}
        for row in self.get_atoms.index:
            graph[row] = self.atom_neighbours(row)

        # Conversion of ring_path into integers
        ring_atoms = list(map(int, ring_path.split()))
        if len(ring_atoms) == 0:
            print("no path provided")
        if len(ring_atoms) < 3:
            print("Invalid input: not enough atoms to form a cycle.")
        if len(ring_atoms) == 1:
            print("Invalid input: only one atom provided.")

        for i in range(len(ring_atoms)):
            if ring_atoms[(i + 1) % len(ring_atoms)] not in graph[ring_atoms[i]]:
                print(f"Invalid cycle: no bond between {current_atom} and {next_atom}.")
        print("Valid cycle.")

    def ring_finding_st(self, ring_path):
        graph = {}
        for row in self.get_atoms.index:
            graph[row] = self.atom_neighbours(row)

        # Conversion of ring_path into integers
        ring_atoms = list(map(int, ring_path.split()))
        if len(ring_atoms) == 0:
            return st.warning("no path provided")
        if len(ring_atoms) < 3:
            return st.warning("Invalid input: not enough atoms to form a cycle.")
        if len(ring_atoms) == 1:
            return st.warning("Invalid input: only one atom provided.")

        for i in range(len(ring_atoms)):
            if ring_atoms[(i + 1) % len(ring_atoms)] not in graph[ring_atoms[i]]:
                return st.warning(f"Invalid cycle: no bond between {current_atom} and {next_atom}.")

        return st.success("Valid cycle !")

#This class will be useful for subsequents methods comparing molecules
class all_molecule:
    def __init__(self):
        self.molecules = {}

    def add_molecule(self, molecule, name):
        self.molecules[name] = {'overall_info':molecule.get_infos, 'atom_bloc' :molecule.get_atoms, 'bond_bloc' : molecule.get_bonds}


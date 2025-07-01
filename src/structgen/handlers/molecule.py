import sys, re
import numpy as np
import random
import copy
from collections import deque

from .atomic import ATOMIC

from structgen import constants

class MOLECULE(ATOMIC):
    
    def __init__(self, args):

        super().__init__(args)

        self.ATOM_STYLE = "MOLECULE"
        self.ATOM_ATTRS = ["atom-id", "molecule-id", "type", "x", "y", "z"]
        self.BOND_ATTRS = ["id", "type", "atom1-id", "atom2-id"]

        self.INT_ATOM_ATTRS = [0, 1, 2]

        molecule_help_msg = "\n\nPlease enter the molecules you wish to generate in the structure\n"
        molecule_help_msg += "in the following format: type(num atoms in molecule)-type(num atoms in molecule).\n"
        molecule_help_msg += "One the next line please enter the number of molecules of that type.\n"
        molecule_help_msg += "There can be more than two atom types in a molecule.\n"
        molecule_help_msg += "Example 1: 1(2)-2(5) | Example 2: 2(10)-3(1)-1(4)\n"
        molecule_help_msg += "Types you have provided:\n"
        molecule_help_msg += "\n".join(f"type {i + 1}: MASS({self.MASSES[i]}) COUNT({self.COUNTS[i]}) " for i in range(self.NUM_TYPES))
        molecule_help_msg += "\nEach new molecule goes on the next line. When done please input empty line by pressing Enter."

        molecule_help_msg += "\n\n"
        print(molecule_help_msg)

        self.MOLECULES = list()
        self.MOLECULE_STOICHEOMETRIES = list()
        got_all_molecules = False
        atoms_in_molecules = self.COUNTS.copy()
        i = 1
        while not got_all_molecules:
            molecule_str = input(f"Molecule {i}: ")
            if molecule_str:
                molecule = list()
                try:
                    for atom_str in re.sub(r'\s+', '', molecule_str).split("-"):
                        atom_data = list(map(int, atom_str[:-1].split("(")))
                        molecule.append(atom_data)
                except:
                    sys.exit(f"ERROR! Incorrect molecule format: {molecule_str}")

                self.MOLECULE_STOICHEOMETRIES.append(copy.deepcopy(molecule))

                mol_count = input(f"Number of molecule {i}: ")

                if not mol_count.isdigit():
                    sys.exit("ERROR! Number of molecules has to be a non-negative integer")

                for type_count_pair in molecule:
                    type_count_pair[1] *= int(mol_count)
                    type_ind = type_count_pair[0] - 1
                    atoms_in_molecules[type_ind] -= type_count_pair[1]

                    if atoms_in_molecules[type_ind] < 0:
                        sys.exit(f"ERROR! Not enough atoms of type {type_count_pair[0]}")

                self.MOLECULES.append(molecule)

                i += 1
            else:
                got_all_molecules = True

        self.NUM_MOLECULES = i - 1
        self.BONDS = None

    def molecule_precondition(self):
        self.__to_distribute = copy.deepcopy(self.MOLECULES)
        self.__mols_left = self.NUM_MOLECULES
        self.__current_mol_id = 1
        self.__current_bond_id = 1
        self.__incomplete_molecules = list()
        self.__num_incomplete = 0
        self.BONDS = list()
            
    def get_atom_data(self):
        self.molecule_precondition()
        return super().get_atom_data()

    def get_extra_atom_props(self, id, type, position):
        for mol_id, mol in enumerate(self.__to_distribute):
            for type_count_pair in mol:
                if type == type_count_pair[0] and type_count_pair[1]:
                    i = 0
                    found_matching_incomplete = False
                    while i < self.__num_incomplete and not found_matching_incomplete:
                        if type in self.__incomplete_molecules[i]:
                            found_matching_incomplete = True
                            self.__incomplete_molecules[i].remove(type)

                        else:
                            i += 1

                    if i == self.__num_incomplete:
                        rest_of_mol_types = list()
                        found_one_of_target_type = False
                        for sp_id, species in enumerate(mol):
                            sp_type, _ = species
                            if sp_type == type and not found_one_of_target_type:
                                rest_of_mol_types += [sp_type for _ in range(self.MOLECULE_STOICHEOMETRIES[mol_id][sp_id][1] - 1)]
                                found_one_of_target_type = True
                            else:
                                rest_of_mol_types += [sp_type for _ in range(self.MOLECULE_STOICHEOMETRIES[mol_id][sp_id][1])]


                        self.__incomplete_molecules.append(rest_of_mol_types)
                        self.__num_incomplete += 1

                    type_count_pair[1] -= 1

                    return [id, i + 1, type]


        return [id, 0, type]

    def get_nbonds(self):
        return 0

    def get_nbond_types(self):
        return 0

    def get_bond_data(self):
        return []

    def get_bond_attrs(self):
        return self.BOND_ATTRS

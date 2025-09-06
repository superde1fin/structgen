import sys, re
import numpy as np
import random
import copy
from collections import deque

from .atomic import ATOMIC

from structgen import constants

class MOLECULAR(ATOMIC):

    def __init__(self, args):

        super().__init__(args)

        self.ATOM_STYLE = "MOLECULE"
        self.ATOM_ATTRS = ["atom-id", "molecule-id", "type", "x", "y", "z"]
        self.BOND_ATTRS = ["id", "type", "atom1-id", "atom2-id"]

        self.INT_ATOM_ATTRS = [0, 1, 2]

        self.BOND_TYPES = 0
        self.BOND_COUNTS = list()
        self.BONDS = list()

        #Scan bond numbers
        for bond_data in args.bond:
            num_attrs = len(bond_data)

            if num_attrs > 2:
                sys.exit(f"ERROR! Too many per-atom arguments ({num_attrs}).")

            pair_str = bond_data[0]
            try:
                b_types = list(map(int, re.sub(r'\s+', '', pair_str).split("-")))
                if b_types[0] < 1 or b_types[1] < 1:
                    raise ValueError
            except:
                sys.exit("ERROR! Atom types making up a bond have to be positive integers")

            self.BONDS.append(b_types)

            if num_attrs >= 2:
                count_str = bond_data[1]

                if count_str.isdigit():
                    self.BOND_COUNTS.append(int(count_str))
                else:
                    sys.exit(f"ERROR!: Non-numeric value for bond count {count_str}")

            else:
                self.BOND_COUNTS.append(constants.DEFAULT_BOND_COUNT)

            self.BOND_TYPES += 1


        self.INTRAMOL = [ct for cr in self.BOND_COUNTS]

        if len(args.intermol) > len(self.INTRAMOL):
            sys.exit("ERROR! Please make sure that no extra --intermol parameters are provided.")

        for i, inter in args.intermol:
            if inter > self.INTRAMOL[i]:
                sys.exit(f"ERROR! Cannot make {inter} bonds intermolecular. Only {self.INTRAMOL[i]} bonds of this type requested.")

            self.INTRAMOL[i] -= inter


        #Request molecule information
        molecule_help_msg = "\n\nPlease enter the molecules you wish to generate in the structure\n"
        molecule_help_msg += "in the following format: type(num atoms in molecule)-type(num atoms in molecule).\n"
        molecule_help_msg += "One the next line please enter the number of molecules of that type.\n"
        molecule_help_msg += "There can be more than two atom types in a molecule.\n"
        molecule_help_msg += "Example 1: 1(2)-2(5) | Example 2: 2(10)-3(1)-1(4)\n"
        molecule_help_msg += "Types you have provided:\n"
        molecule_help_msg += "\n".join(f"type {i + 1}: MASS({self.MASSES[i]}) COUNT({self.COUNTS[i]}) " for i in range(self.NUM_TYPES))
        molecule_help_msg += "\nEach new molecule goes on the next line. When done please input empty line by pressing Enter.\n"
        molecule_help_msg += "If the total number of atoms is less than the number of atoms in provided molecules each of the leftover atoms will be assigned a separate molecule id.\n"

        molecule_help_msg += "\n\n"
        print(molecule_help_msg)

        self.NUM_MOLECULES = 0
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

                mol_count = int(mol_count)

                self.NUM_MOLECULES += mol_count

                for type_count_pair in molecule:
                    type_count_pair[1] *= mol_count
                    type_ind = type_count_pair[0] - 1
                    atoms_in_molecules[type_ind] -= type_count_pair[1]

                    if atoms_in_molecules[type_ind] < 0:
                        sys.exit(f"ERROR! Not enough atoms of type {type_count_pair[0]}")

                self.MOLECULES.append(molecule)

                i += 1
            else:
                got_all_molecules = True

        self.EXTRA_MOL_TYPES = dict()

        self.NUM_MOLECULE_TYPES = i - 1


    def molecule_precondition(self):
        self.__extra_mol_ctr = self.NUM_MOLECULES
        self.__to_distribute = copy.deepcopy(self.MOLECULES)
        self.__intramol_counter = copy.deepcopy(self.INTRAMOL)
        self.__intermol = [count - self.INTRAMOL[i] for i, count in enumerate(self.BOND_COUNTS)]
        self.__current_bond_id = 1
        self.__molecules_leftover = list()
        self.__incomplete_molecules = list()
        self.__num_incomplete = 0
        self.__available_for_non_mo2mo_bonding = list()
        self.__bond_data = list()
            
    #Have to modify so that atoms are produced in an a molecular fashion
    def get_atom_data(self):
        self.molecule_precondition()
        return super().get_atom_data()

    def get_extra_atom_props(self, id, type, position):
        for mol_id, mol in enumerate(self.__to_distribute):
            for type_count_pair in mol:
                if type == type_count_pair[0] and type_count_pair[1]:
                    self.__available_for_non_mo2mo_bonding.append((id, type))
                    i = 0
                    found_matching_incomplete = False
                    while i < self.__num_incomplete and not found_matching_incomplete:

                        available_bonds = [[] for _ in self.BOND_COUNTS]
                        for prev_id, prev_type in self.__incomplete_molecules[i]:
                            for btype_ctr, bpair in self.BOND_COUNTS:
                                if [prev_type, type] == bpair or [type, prev_type] == bpair:
                                    available_bonds[btype_ctr].append((prev_id, id))

                        if type in self.__molecules_leftover[i]:
                            found_matching_incomplete = True
                            self.__molecules_leftover[i].remove(type)
                            self.__incomplete_molecules[i].append((id, type))

                            for av_bond_ctr, av_bond_pairs in enumerate(available_bonds):
                                while av_bond_pairs and self.__intramol_counter[av_bond_ctr]:
                                    bond_pair = bond_pairs.pop()
                                    self.__bond_data.append([self.__current_bond_id, av_bond_ctr + 1,  *bond_pair])
                                    self.__current_bond_id += 1
                                    self.__intramol_counter[av_bond_ctr] -= 1

                        else:
                            for av_bond_ctr, av_bond_pairs in enumerate(available_bonds):
                                while av_bond_pairs and self.__intermol[av_bond_ctr]:
                                    bond_pair = bond_pairs.pop()
                                    self.__bond_data.append([self.__current_bond_id, av_bond_ctr + 1,  *bond_pair])
                                    self.__current_bond_id += 1
                                    self.__intermol[av_bond_ctr] -= 1

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


                        self.__molecules_leftover.append(rest_of_mol_types)
                        self.__incomplete_molecules.append([(id, type)])
                        self.__num_incomplete += 1

                    type_count_pair[1] -= 1

                    return [id, i + 1, type]

        #Have to handle cases where one (ore more) of the atoms in a bond are not in a molecule

        self.__extra_mol_ctr += 1
        return [id, self.__extra_mol_ctr, type]

    def get_nbonds(self):
        return 0

    def get_nbond_types(self):
        return 0

    def get_bond_data(self):
        return []

    def get_bond_attrs(self):
        return self.BOND_ATTRS

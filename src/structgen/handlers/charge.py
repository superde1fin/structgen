import sys
import numpy as np
import random

from structgen import constants

class CHARGE:
    def __init__(self, args):
        self.ATOM_STYLE = "CHARGE"
        self.ATOM_ATTRS = ["id", "type", "q", "x", "y", "z"]

        self.MASSES = list()
        self.COUNTS = list()
        self.CHARGES = list()
        self.NATOMS = 0
        self.TOTAL_MASS = 0
        self.NUM_TYPES = 0
        self.ATOM_DATA = None

        if args.factor < 1:
            sys.exit("ERROR! Factor must be a positive integer.")

        for atom_data in args.atom:
            num_attrs = len(atom_data)

            if num_attrs > 3:
                sys.exit(f"ERROR! Too many per-atom arguments ({num_attrs}).")

            self.NUM_TYPES += 1

            mass_str = atom_data[0]
            try:
                mass_val = float(mass_str)

                if mass_val <= 0:
                    raise ValueError

                self.MASSES.append(mass_val)
            except ValueError:
                sys.exit(f"ERROR! Non‑numeric or non-positive value for atom mass: '{mass_str}'.")


            if num_attrs >= 2:
                count_str = atom_data[1]

                if count_str.isdigit():
                    add_atoms = args.factor*int(count_str)
                else:
                    sys.exit(f"ERROR! Non‑numeric value for atom mass: '{mass_str}'.")

            else:
                add_atoms = constants.DEFAULT_ATOM_COUNT*args.factor


            if num_attrs == 3:
                charge_str = atom_data[2]

                try:
                    self.CHARGES.append(float(charge_str))
                except ValueError:
                    sys.exit(f"ERROR! Non‑numeric or non-positive value for atom mass: '{mass_str}'.")
            else:
                self.CHARGES.append(constants.DEFAULT_ATOM_CHARGE)

            self.COUNTS.append(add_atoms)
            self.NATOMS += add_atoms
            self.TOTAL_MASS += (mass_val * add_atoms)

        volume = self.TOTAL_MASS/(args.density * 0.6022136) #Conversion from g to amu
        self.REGION_SIDE = volume**(1/3)

        self.BUFFER = self.REGION_SIDE*args.buffer/100

        self.OUTFILE = args.output

    def get_atom_style(self):
        return self.ATOM_STYLE.lower()

    def get_atom_attrs(self):
        return self.ATOM_ATTRS

    def get_natoms(self):
        return self.NATOMS

    def get_natom_types(self):
        return self.NUM_TYPES

    def get_sim_region_sides(self):
        return [self.REGION_SIDE]*3

    def get_outfile_name(self):
        return self.OUTFILE

    def get_atom_data(self):
        
        atoms_per_side = int(self.NATOMS**(1/3)) + 1 #Round up
        pos_options = np.linspace(self.BUFFER, self.REGION_SIDE - self.BUFFER, atoms_per_side)

        pos_visited = list()

        atoms_counter = 0
        atoms_added = 0 #Sum of self.COUNTS up to the current type

        self.ATOM_DATA = np.empty((self.NATOMS, len(self.ATOM_ATTRS)))

        for i in range(self.NUM_TYPES):
            while atoms_counter < self.COUNTS[i] + atoms_added:
                found_unused_pos = False
                while not found_unused_pos:
                    x_rand = random.randint(0, atoms_per_side - 1)
                    y_rand = random.randint(0, atoms_per_side - 1)
                    z_rand = random.randint(0, atoms_per_side - 1)
                    cur_pos = (pos_options[x_rand], pos_options[y_rand], pos_options[z_rand])
                    if cur_pos not in pos_visited:
                        found_unused_pos = True

                self.ATOM_DATA[atoms_counter] = np.array([atoms_counter + 1, i + 1, self.CHARGES[i], *cur_pos])
                pos_visited.append(cur_pos)
                atoms_counter += 1

            atoms_added += self.COUNTS[i]

        return self.ATOM_DATA

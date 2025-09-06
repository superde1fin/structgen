import sys
import numpy as np
import random

from .atomic import ATOMIC

from structgen import constants

class CHARGE(ATOMIC):

    extra_atom_args = 1 #Charge added

    def __init__(self, args):

        super().__init__(args)

        self.ATOM_STYLE = "CHARGE"
        self.ATOM_ATTRS = ["id", "type", "q", "x", "y", "z"]

        self.CHARGES = list()

        atom_args_limit = self.get_max_args()

        for atom_data in args.atom:
            num_attrs = len(atom_data)

            if num_attrs > atom_args_limit:
                sys.exit(f"ERROR! Too many per-atom arguments ({num_attrs}).")

            if num_attrs == atom_args_limit:
                charge_str = atom_data[2]

                try:
                    self.CHARGES.append(float(charge_str))
                except ValueError:
                    sys.exit(f"ERROR! Non‑numeric  value for atom charge: '{charge_str}'.")
            else:
                self.CHARGES.append(constants.DEFAULT_ATOM_CHARGE)

    def get_max_args(self):
        return super().get_max_args() + self.extra_atom_args

    def get_extra_atom_props(self, id, type, position):
        return [id, type, self.CHARGES[type - 1]]

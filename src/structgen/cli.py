import sys
import numpy as np
import random
import argparse

#from .handlers import CHARGE, ATOMIC, FULL, MOLECULAR
from .handlers import CHARGE, ATOMIC
from . import constants
from . import utils

def parse_args():
    parser = utils.ErrorHandlingParser(
        prog="structgen",
        description = (
            "Generate initial structures in LAMMPS format for glassy systems "
            "given density, atom masses and counts."
        ),
        formatter_class = argparse.RawTextHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest = "command", required = True, metavar = "atom_style", action = utils.StrictSubParsersAction)

    atomic = subparsers.add_parser("atomic", help = "Argument parser for generating initial glass structures in lammps format compliant with atom_style ATOMIC.", formatter_class = utils.NoMetavarHelpFormatter)

    atomic.add_argument(
        "-a", "--atom",
        metavar = "",
        nargs = "+",
        action = "append",
        required = True,
        help = (
            "Provide <mass> [<count>] for an atom type.\n"
            "Defaults: count=1\n"
            "Repeat for multiple atom types: -a 28.085 -a 16.00 2"
        ),
    )

    atomic.add_argument(
        "-d", "--density",
        type = float,
        required = True,
        metavar = "",
        help = "Desired density of the output structure in g/cm^3.",
    )

    atomic.add_argument(
        "-b", "--buffer",
        type=float,
        default = constants.DEFAULT_BUFFER_PERCENT,
        metavar = "",
        help = "Spacing between simulation region walls in %% of box size (default: 5).",
    )

    atomic.add_argument(
        "-f", "--factor",
        type = int,
        default = constants.DEFAULT_ATOM_FACTOR,
        metavar = "",
        help = "Scale factor applied to all atom counts (default: 1).",
    )

    atomic.add_argument(
        "-o", "--output",
        default = constants.DEFAULT_FILENAME,
        metavar = "",
        help = "Name of the resulting structure file (default: initial.structure).",
    )

    atomic.set_defaults(handler_class = ATOMIC)

    charge = subparsers.add_parser("charge", help = "Argument parser for generating initial glass structures in lammps format compliant with atom_style CHARGE.", formatter_class = utils.NoMetavarHelpFormatter)

    charge.add_argument(
        "-a", "--atom",
        metavar = "",
        nargs = "+",
        action = "append",
        required = True,
        help = (
            "Provide <mass> [<count> <charge>] for an atom type.\n"
            "Defaults: count=1, charge=0.0.\n"
            "Repeat for multiple atom types: -a 28.085 -a 16.00 2 -1.2"
        ),
    )


    charge.add_argument(
        "-d", "--density",
        type = float,
        required = True,
        metavar = "",
        help = "Desired density of the output structure in g/cm^3.",
    )

    charge.add_argument(
        "-b", "--buffer",
        type=float,
        default = constants.DEFAULT_BUFFER_PERCENT,
        metavar = "",
        help = "Spacing between simulation region walls in %% of box size (default: 5).",
    )

    charge.add_argument(
        "-f", "--factor",
        type = int,
        default = constants.DEFAULT_ATOM_FACTOR,
        metavar = "",
        help = "Scale factor applied to all atom counts (default: 1).",
    )

    charge.add_argument(
        "-o", "--output",
        default = constants.DEFAULT_FILENAME,
        metavar = "",
        help = "Name of the resulting structure file (default: initial.structure).",
    )

    charge.set_defaults(handler_class = CHARGE)


    """
    full = subparsers.add_parser("full", help = "Argument parser for generating initial glass structures in lammps format compliant with atom_style FULL.")

    full.add_argument(
        "-a", "--atom",
        metavar = "",
        nargs = "+",
        action = "append",
        required = True,
        help = (
            "Provide <mass> [<count> <charge>] for an atom type.\n"
            "Defaults: count=1, charge=0.0.\n"
            "Repeat for multiple atom types: -a 28.085 -a 16.00 2 -1.2"
        ),
    )

    full.add_argument(
        "-d", "--density",
        type = float,
        required = True,
        metavar = "",
        help = "Desired density of the output structure in g/cm^3.",
    )

    full.add_argument(
        "-b", "--buffer",
        type=float,
        default = constants.DEFAULT_BUFFER_PERCENT,
        metavar = "",
        help = "Spacing between simulation region walls in %% of box size (default: 5).",
    )

    full.add_argument(
        "-f", "--factor",
        type = int,
        default = constants.DEFAULT_ATOM_FACTOR,
        metavar = "",
        help = "Scale factor applied to all atom counts (default: 1).",
    )

    full.add_argument(
        "-o", "--output",
        default = constants.DEFAULT_FILENAME,
        metavar = "",
        help = "Name of the resulting structure file (default: initial.structure).",
    )

    full.set_defaults(handler_class = FULL)

    molecular = subparsers.add_parser("molecular", help = "Argument parser for generating initial glass structures in lammps format compliant with atom_style MOLECULAR.")

    molecular.add_argument(
        "-a", "--atom",
        metavar = "",
        nargs = "+",
        action = "append",
        required = True,
        help = (
            "Provide <mass> [<count>] for an atom type.\n"
            "Defaults: count=1.\n"
            "Repeat for multiple atom types: -a 28.085 -a 16.00 2"
        ),
    )

    molecular.add_argument(
        "-b", "--bond",
        metavar = "",
        nargs = "+",
        action = "append",
        help = (
            "Provide <type pair> [<count>] specifying a number of bonds of a particular type.\n"
            "Types are given to atoms according to the order specified (starting at 1)\n"
            "Defaults: count=1.\n"
            "Repeat for multiple bond types: -b 1-2 10 -b 3-4"
        ),
    )

    molecular.add_argument(
        "-d", "--density",
        type = float,
        required = True,
        metavar = "",
        help = "Desired density of the output structure in g/cm^3.",
    )

    molecular.add_argument(
        "-bf", "--buffer",
        type=float,
        default = constants.DEFAULT_BUFFER_PERCENT,
        metavar = "",
        help = "Spacing between simulation region walls in %% of box size (default: 5).",
    )

    molecular.add_argument(
        "-f", "--factor",
        type = int,
        default = constants.DEFAULT_ATOM_FACTOR,
        metavar = "",
        help = "Scale factor applied to all atom counts (default: 1).",
    )

    molecular.add_argument(
        "--intermol",
        type = int,
        nargs = "+",
        metavar = "",
        help = (
            "If specified this argument determines the number of bonds of each type that will be assigned to the intermolecular pairs. "
            "The number of values provided should be less than or equal to the number of bonds requested. "
            "When omitted this is NOT the default for all bond types meaning that structgen will attempt to assign all bonds to intramolecular pairs."
        ),
    )

    molecular.add_argument(
        "-o", "--output",
        default = constants.DEFAULT_FILENAME,
        metavar = "",
        help = "Name of the resulting structure file (default: initial.structure).",
    )

    molecular.set_defaults(handler_class = MOLECULAR)
    """

    return parser.parse_args()

def gen_header(handler):
    handler_style = handler.get_atom_style()

    header = str()
    header += "#Generated by structgen utility git@github.com:superde1fin/structgen.git\n\n"

    header += f"{handler.get_natoms()} atoms\n"
    if handler_style in constants.NEED_BONDS:
        header += f"{handler.get_nbonds()} bonds\n"

    header += f"{handler.get_natom_types()} atom types\n"
    if handler_style in constants.NEED_BONDS:
        header += f"{handler.get_nbond_types()} bond types\n"

    sides = handler.get_sim_region_sides()
    header += f"\n0 {sides[0]} xlo xhi"
    header += f"\n0 {sides[1]} ylo yhi"
    header += f"\n0 {sides[2]} zlo zhi"

    return header


def gen_atoms(handler):
    atoms = str()
    atoms += "\n\n"
    column_headers = ' '.join(handler.get_atom_attrs())
    atoms += f"Atoms # {handler.get_atom_style()} | {column_headers}\n\n"

    atom_data = handler.get_atom_data()

    #Add atom lines and convert necesary elements to integers
    indices_to_int = handler.get_int_atom_indices()

    atoms += "\n".join(" ".join(str(int(x)) if i in indices_to_int else str(x) for i, x in enumerate(row)) for row in atom_data)

    return atoms

def gen_bonds(handler):
    bonds = str()
    bonds += "\n\n"
    column_headers = " ".join(handler.get_bond_attrs())
    bonds += f"Bonds # {column_headers}"

    bond_data = handler.get_bond_data()

    bonds += "\n".join(" ".join(str(x) for x in row) for row in bond_data)

    return bonds


def main():
    args = parse_args()

    handler = args.handler_class(args)

    handler_style = handler.get_atom_style()

    file_str = str()

    file_str += gen_header(handler)
    file_str += gen_atoms(handler)

    if handler_style in constants.NEED_BONDS:
        file_str += gen_bonds(handler)

    with open(handler.get_outfile_name(), "w") as file:
        file.write(file_str)

        

if __name__ == "__main__":
    main()

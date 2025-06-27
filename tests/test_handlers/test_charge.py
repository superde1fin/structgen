from types import SimpleNamespace
from structgen.handlers.charge import CHARGE

def test_charge_handler_atom_counts():
    args = SimpleNamespace(
        atom=[["28.085"], ["16.00", "2", "-1.2"]],
        density=2.2,
        buffer=5.0,
        factor=1,
        output="test.structure"
    )
    handler = CHARGE(args)
    assert handler.get_natoms() == 3
    assert handler.get_natom_types() == 2


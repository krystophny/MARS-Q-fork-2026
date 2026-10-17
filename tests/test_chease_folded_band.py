"""Folded active equations solve independent constrained native systems."""
from pathlib import Path
from support import NativeTest, NATIVE, routine

PRECISION = "module prec_const\ninteger,parameter::rkind=kind(1.d0)\nend module\n"


class CheaseFoldedBand(NativeTest):
    def run_oracle(self, name):
        files = {"precision.f90": PRECISION,
                 "DECLAR.inc": (NATIVE / "DECLAR.inc").read_text(),
                 "folded_band.f90": (NATIVE / "folded_band.f90").read_text()}
        files["native.f"] = "\n".join(routine(name) for name in
                                     ("ALDLT", "LYV", "DWY", "LTXW", "SAXPY", "SDOT"))
        files["oracle.f90"] = Path(__file__).with_name(name).read_text()
        self.execute(self.compile(files))

    def test_physical_cells_seam_center_boundary_and_guards(self):
        self.run_oracle("chease_folded_cells.f90")

    def test_scaled_physical_pivots_with_inactive_unit_constraints(self):
        self.run_oracle("chease_scaled_boundary.f90")

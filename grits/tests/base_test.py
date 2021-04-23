from os import path
import pytest

import mbuild as mb

from grits import CG_Compound


test_dir = path.dirname(__file__)


class BaseTest:
    @pytest.fixture
    def p3ht(self):
        p3ht = mb.load(path.join(test_dir, "assets/P3HT_16.mol2"))
        return p3ht

    @pytest.fixture
    def methane(self):
        methane =  mb.load("C", smiles=True)
        return methane

    @pytest.fixture
    def cg_methane(self, methane):
        cg_beads = [("_A", "C")]

        cg_methane = CG_Compound(methane, cg_beads)
        return cg_methane

    @pytest.fixture
    def cg_p3ht(self, p3ht):
        cg_beads = [("_B", "c1sccc1"), ("_S", "CCC")]

        cg_p3ht = CG_Compound(p3ht, cg_beads)
        return cg_p3ht

import pytest
from base_test import BaseTest
from mbuild import Compound
import gsd.hoomd

from grits import backmap, CG_System
from grits.utils import amber_dict
from os import path

asset_dir = path.join(path.dirname(__file__), "assets")

class Test_Backmap(BaseTest):
    def test_backmapnobonds(self, methane, cg_methane):
        fg_methane = backmap(cg_methane)

        assert isinstance(fg_methane, Compound)
        assert fg_methane.n_particles == methane.n_particles

    def test_backmapbonds(self, p3ht, cg_p3ht):
        fg_p3ht = backmap(cg_p3ht)

        assert isinstance(fg_p3ht, Compound)
        assert fg_p3ht.n_particles == p3ht.n_particles
        assert fg_p3ht.n_bonds == p3ht.n_bonds

    def test_alkane(self, alkane, cg_alkane):
        fg_alkane = backmap(cg_alkane)

        assert fg_alkane.n_bonds == alkane.n_bonds
        assert fg_alkane.n_particles == alkane.n_particles

    def test_backmapsystem(self, tmp_path):
        gsdfile = path.join(asset_dir, "benzene-aa.gsd")
        cg_system = CG_System(
            gsdfile,
            beads={"_B": "c1ccccc1"},
            conversion_dict=amber_dict,
            mass_scale=12.011,
        )
        cg_filename = 'benzene-cg.gsd'#path.join(tmp_path, "benzene-cg.gsd")
        cg_system.save(cg_filename)
        backmapped_system = backmap(cg_system, cg_filename)
        fg_filename = cg_filename.split('.gsd')[0] + "-finegrained.gsd"
        #path.join(tmp_path, "benzene-fg.gsd")
        #backmapped_system.save(fg_filename)
        with gsd.hoomd.open(gsdfile, "r") as orig, gsd.hoomd.open(
                fg_filename, "r") as backmapped:
            # should have same number of frames
            assert len(orig) == len(backmapped)
            # should have same number of atoms, bonds, etc
            for frame1, frame2 in zip(orig, backmapped):
                assert frame1.particles.N == frame2.particles.N
                

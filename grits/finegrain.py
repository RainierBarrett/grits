"""GRiTS: Fine-graining tools."""
__all__ = ["backmap"]

import itertools as it
import gsd.hoomd
import freud
import numpy as np
from collections import defaultdict

from mbuild import Compound, Particle, load, formats, conversion
import mbuild.formats.hoomd_snapshot

from grits.utils import align, get_hydrogen, get_index
from grits import CG_System


def backmap(cg_compound, cg_gsd_filename=None):
    """Backmap a fine-grained representation onto a coarse one.

    Creates a fine-grained compound from a coarse one using the attributes
    created during CG_Compound initialization.

    Parameters
    ----------
    cg_compound : CG_Compound
        Coarse-grained compound
    cg_gsd_filename : string
        location of GSD file of CG trajectory to backmap
        NOTE: this is *NOT* the gsd file corresponding to
        the CG_Compound, which is the original (FG) trajectory

    Returns
    -------
    :py:class:`mbuild.Compound`
        The atomistic structure mapped onto the coarse-grained one.
    """

    if type(cg_compound) is CG_System:
        assert cg_gsd_filename is not None
               
        # TODO: change logic to save a single gsd file and/or return one single mbuild Compound
        backmapped_system = Compound()
        # get all gsd frames
        cg_gsd_name = cg_gsd_filename.split('.')[0]
        fg_gsd_path = cg_gsd_name + '-finegrained' + '.gsd'
        with gsd.hoomd.open(cg_gsd_filename, 'r') as coarse, gsd.hoomd.open(fg_gsd_path,'w') as fine:
            for frame_i, cg_frame in enumerate(coarse):
                new_frame = gsd.hoomd.Frame()
                position = []
                mass = []
                #TODO: does this need to exist yet?
                orientation = [] if cg_compound.aniso_beads else None
                fine_grained = Compound()
                f_box = freud.Box.from_box(cg_frame.configuration.box)
                unwrap_pos = f_box.unwrap(
                    cg_frame.particles.position, cg_frame.particles.image
                )
                typeid = []
                types = [i.split("...")[0] for i in cg_compound.mapping]
                smarts_strings = [i.split("...")[1] for i in cg_compound.mapping]
                for i, inds in enumerate(cg_compound.mapping.values()):
                    typeid.append(np.ones(len(inds), dtype=np.int64) * i)
                typeid = np.hstack(typeid)
                np.savetxt("typeid.txt",typeid)
                print(f"DEBUG. UNWRAP_POS IS {unwrap_pos} ({unwrap_pos.shape})")
                for i, pos in enumerate(unwrap_pos):
                    print(f"DEBUG. MADE IT TO {i} ({pos})")
                    smarts = smarts_strings[typeid[i]]
                    b = load(smarts, smiles=True)
                    b.translate_to(pos)
                    # TODO: translate orientations if included
                    # TODO: figure out bonds
                    fine_grained.add(b, str(i))
                parmed_compound = conversion.to_parmed(fine_grained,)
                new_frame.configuration.box = cg_frame.configuration.box
                new_frame.configuration.step = frame_i
                new_frame.particles.N = fine_grained.n_particles
                particles = [particle for particle in fine_grained.particles()]
                new_frame.particles.position = [particle.pos for particle in particles]
                types = list(set([particle.name for particle in particles]))
                print(f"DEBUG: types is {types}")
                typeid_vals = [i for i, _ in enumerate(types)]
                # TODO: typeids?, masses? (infer these?)
                new_frame.particles.types = types

                print(f"DEBUG: dir(new_frame): {dir(new_frame)}")
                fine.append(new_frame)
                
        #for compound in cg_compound._compounds:
        #    backmapped_system.add(backmap(compound))
        # save to a gsd here?
        return

    def fg_particles():
        """Set the particles of the fine-grained structure."""
        fine_grained = Compound()

        anchors = dict()
        for i, bead in enumerate(cg_compound):
            smiles = bead.smarts
            b = load(smiles, smiles=True)
            b.translate_to(bead.pos)
            anchors[i] = dict()
            if cg_compound.anchors is not None:
                for index in cg_compound.anchors[bead.name]:
                    anchors[i][index] = b[index]
            fine_grained.add(b, str(i))
        return fine_grained, anchors

    def fg_bonds():
        """Set the bonds for the fine-grained structure."""
        bonded_atoms = []
        remove_hs = []
        rotated = {k: False for k in anchors.keys()}
        for name, inds in cg_compound.bond_map:
            for ibead, jbead in cg_compound.bonds():
                names = [ibead.name, jbead.name]
                if "-".join(names) == name:
                    fi, fj = inds
                elif "-".join(names[::-1]) == name:
                    fj, fi = inds
                else:
                    continue

                i = get_index(cg_compound, ibead)
                j = get_index(cg_compound, jbead)
                try:
                    iatom = anchors[i].pop(fi)
                except KeyError:
                    fi = [x for x in inds if x in anchors[i]][0]
                    iatom = anchors[i].pop(fi)
                try:
                    jatom = anchors[j].pop(fj)
                except KeyError:
                    fj = [x for x in inds if x in anchors[j]][0]
                    jatom = anchors[j].pop(fj)

                hi = get_hydrogen(fine_grained, iatom)
                hj = get_hydrogen(fine_grained, jatom)
                # each part can be rotated
                if not rotated[i]:
                    # rotate
                    align(fine_grained[str(i)], hi, jbead)
                    rotated[i] = True
                if not rotated[j]:
                    # rotate
                    align(fine_grained[str(j)], hj, ibead)
                    rotated[j] = True

                fine_grained.add_bond((iatom, jatom))

                bonded_atoms += (iatom, jatom)
                remove_hs += (hi, hj)

        for atom in remove_hs:
            fine_grained.remove(atom)
        return fine_grained

    fine_grained, anchors = fg_particles()

    if cg_compound.bond_map is None:
        return fine_grained

    fine_grained = fg_bonds()
    return fine_grained

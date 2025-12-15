"""
`test_misc` contains unit tests of auxiliary functions.

Test data is stored in `./ToyData/Simulations.2`

-------------------------------------------------------------------------------
NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import os

import pytest
import pytest_check as check

# run only on sim2 mocking data
pytestmark = [pytest.mark.nodata, pytest.mark.min]


def test_uname2element():
    """Test uname2element function."""
    from fairmd.lipids.auxiliary.elements import uname2element

    check.equal(uname2element("M_C1_M"), "C")
    check.equal(uname2element("M_G1_M"), "C")
    check.equal(uname2element("M_C1N3_M"), "N")
    check.equal(uname2element("M_X_M"), "Dummy")
    check.equal(uname2element("M_D_M"), "Dummy")

    with pytest.raises(KeyError):
        uname2element("UnknownElement")


@pytest.fixture
def tip4p_trajectory(tmpdir):
    import shutil

    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.coordinates.memory import MemoryReader

    folder = str(tmpdir)

    pdb_content = """TITLE     Single 4-site water
CRYST1   63.646   63.646   86.762  90.00  90.00  90.00 P 1           1
ATOM      1  OW  wate    1       4.370  31.010  61.030  1.00  0.00
ATOM      2  HW1 wate    1       3.800  30.820  61.780  1.00  0.00
ATOM      3  HW2 wate    1       5.250  30.760  61.320  1.00  0.00
ATOM      4  MW  wate    1       4.410  30.950  61.170  1.00  0.00
END
"""
    tmpfile = os.path.join(folder, "test_water_maicos.pdb")
    with open(tmpfile, "w") as f:
        f.write(pdb_content)
    u = mda.Universe(
        tmpfile,
        format="PDB",
    )
    nul_dim = u.dimensions
    # Store coordinates for each frame
    n_frames = 10
    trajectory = np.zeros((n_frames, u.atoms.n_atoms, 3), dtype=np.float32)

    for i in range(n_frames):
        coords = u.atoms.positions.copy()
        # move all atoms a little bit in x direction
        coords = np.add(coords, [0.1 * i, 0.0, 0.0])
        trajectory[i] = coords

    # Assign trajectory to Universe
    u.load_new(trajectory, order="fac", format=MemoryReader)
    for ts in u.trajectory:
        ts.dimensions = nul_dim

    yield u, folder
    # tear down
    shutil.rmtree(folder)


def test_maicos_interface(tip4p_trajectory):
    from fairmd.lipids.maicos import DensityPlanar

    u, folder = tip4p_trajectory

    # Now we are done!
    u.add_TopologyAttr("elements")
    u.atoms.elements = ["O", "H", "H", "Dummy"]  # Assign elements manually

    # Skip unwrap/pack for speed - trajectories are already centered and whole
    base_options = {"unwrap": False, "bin_width": 1, "pack": False}
    zlim = {"zmin": 0, "zmax": 8.67623}
    dens_options = {**zlim, **base_options}
    ofname = os.path.join(folder, "DiporderWater.json")
    # Simulate MAICoS calls
    cos_water = DensityPlanar(
        u.atoms,
        dens="electron",
        **dens_options,
        output=ofname,
    )
    cos_water.run()
    cos_water.save()

    fexists = os.path.isfile(ofname)
    assert fexists


def test_maicos_what_to_compute(caplog, logger):
    import logging

    from fairmd.lipids import RCODE_ERROR
    from fairmd.lipids.analyze import computeMAICOS
    from fairmd.lipids.core import System

    s = System(
        data={
            "DOI": "00.0000/abcd",
            "path": "does-not-matter",
            "TYPEOFSYSTEM": "lipid bilayer",
            "SOFTWARE": "GROMACS",
            "TPR": [["file.tpr"]],
            "TRJ": [["file.xtc"]],
            "COMPOSITION": {
                "SOL": {
                    "NAME": "SPC",
                    "MAPPING": "mappingSPCwater.yaml",
                    "COUNT": 1000,
                },
                "DPPC": {
                    "NAME": "DPPC",
                    "MAPPING": "mappingDPPCberger.yaml",
                    "COUNT": 100,
                },
            },
        },
    )
    caplog.clear()
    logger.info("Testing MAICOS interface against GROMACS-like setup")
    with caplog.at_level(logging.INFO):
        rcode = computeMAICOS(s, logging.getLogger("test_logger"), ffonly=False)
    for line in caplog.text.splitlines():
        if "Files to be computed:" in line:
            break
    check.is_in("TotalDensity", line)
    check.is_in("DiporderWater", line)
    check.equal(rcode, RCODE_ERROR)
    # Testing maicos-default ffonly=True
    caplog.clear()
    logger.info("Testing MAICOS default-mode interface against GROMACS-like setup")
    with caplog.at_level(logging.INFO):
        rcode = computeMAICOS(s, logging.getLogger("test_logger"))
    for line in caplog.text.splitlines():
        if "Files to be computed:" in line:
            break
    check.is_in("TotalDensity", line)
    check.is_not_in("DiporderWater", line)
    check.equal(rcode, RCODE_ERROR)
    # Now it will be NAMD-like setup
    logger.info("Testing MAICOS interface against NAMD-like setup")
    del s["TPR"]
    s["PDB"] = [["file.pdb"]]
    s["SOFTWARE"] = "NAMD"
    # and we will not compute DiporderWater, Dielectric and ChargeDensity

    caplog.clear()
    with caplog.at_level(logging.INFO):
        rcode = computeMAICOS(s, logging.getLogger("test_logger"), ffonly=False)
    for line in caplog.text.splitlines():
        if "Files to be computed:" in line:
            break
    check.is_in("TotalDensity", line)
    check.is_not_in("DiporderWater", line)
    check.is_not_in("Dielectric", line)
    check.is_not_in("ChargeDensity", line)
    check.equal(rcode, RCODE_ERROR)

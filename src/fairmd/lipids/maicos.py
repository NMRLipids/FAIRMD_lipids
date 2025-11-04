"""
Wrappers for MAICoS calculations adapted to Databank needs.

- Checks if a system is suitable for maicos calculations
- Custom JSON encoder for numpy arrays
- Custom maicos analysis classes with save methods adapted to Databank needs
"""

import contextlib
import json
import os
import re
import subprocess
from collections import deque
from logging import Logger

import maicos
import MDAnalysis as mda
import numpy as np
from maicos.core import ProfilePlanarBase
from maicos.lib.weights import density_weights

from fairmd.lipids.core import System
from fairmd.lipids.jsonEncoders import CompactJSONEncoder
from fairmd.lipids.settings.molecules import lipids_set


def is_system_suitable_4_maicos(system: System) -> bool:
    """
    Check if the system is suitable for maicos calculations."

    :param system: Databank System object (System)
    :return: False if system should be skipped
    """
    if system["TYPEOFSYSTEM"] == "miscellaneous":
        return False
    try:
        if system["WARNINGS"]["ORIENTATION"]:
            print("Skipping due to ORIENTATION warning:", system["WARNINGS"]["ORIENTATION"])
            return False
    except (KeyError, TypeError):
        pass
    try:
        if system["WARNINGS"]["PBC"] == "hexagonal-box":
            print("Skipping due to PBC warning:", system["WARNINGS"]["PBC"])
            return False
    except (KeyError, TypeError):
        pass
    try:
        if system["WARNINGS"]["NOWATER"]:
            print("Skipping because there is not water in the trajectory.")
            return False
    except (KeyError, TypeError):
        pass
    return True


def first_last_carbon(system: System, logger: Logger) -> tuple[str, str]:
    """Find last carbon of sn-1 tail and g3 carbon."""
    g3_atom = ""
    last_atom = ""
    for molecule in system["COMPOSITION"]:
        if molecule in lipids_set:
            mapping = system.content[molecule].mapping_dict

            # TODO: rewrite via lipid dictionary!
            for nm in ["M_G3_M", "M_G13_M", "M_C32_M"]:
                _ga = mapping.get(nm, {}).get("ATOMNAME")
                g3_atom = _ga if _ga else g3_atom

            # TODO: rewrite via lipid dictionary
            for c_idx in range(4, 30):
                if "M_G1C4_M" in mapping:
                    atom = "M_G1C" + str(c_idx) + "_M"
                elif "M_G11C4_M" in mapping:
                    atom = "M_G11C" + str(c_idx) + "_M"
                elif "M_CA4_M" in mapping:
                    atom = "M_CA" + str(c_idx) + "_M"
                else:
                    # cannot be determined for this particular lipid. Maybe another ..
                    break
                _la = mapping.get(atom, {}).get("ATOMNAME")
                last_atom = _la if _la else last_atom
    logger.info(f"Found last atom {last_atom} and g3 atom {g3_atom} for system {system['ID']}")

    return (last_atom, g3_atom)


def traj_centering_for_maicos(
    system_path: str,
    trj_name: str,
    tpr_name: str,
    last_atom: str,
    g3_atom: str,
    eq_time: int = 0,
    *,
    recompute: bool = False,
) -> str:
    """Center trajectory around the center of mass of all methyl carbons."""
    xtccentered = os.path.join(system_path, "centered.xtc")
    if os.path.isfile(xtccentered) and not recompute:
        return xtccentered  # already done
    if recompute:
        with contextlib.suppress(FileNotFoundError):
            os.remove(xtccentered)

    # make index
    # TODO refactor to MDAnalysis
    ndxpath = os.path.join(system_path, "foo.ndx")
    try:
        echo_input = f"a {last_atom}\nq\n".encode()
        subprocess.run(["gmx", "make_ndx", "-f", tpr_name, "-o", ndxpath], input=echo_input, check=True)
    except subprocess.CalledProcessError as e:
        msg = f"Subprocess failed during ndx file creation: {e}"
        raise RuntimeError(msg) from e
    try:
        with open(ndxpath) as f:
            last_lines = deque(f, 1)
        last_atom_id = int(re.split(r"\s+", last_lines[0].strip())[-1])
        with open(ndxpath, "a") as f:
            f.write("[ centralAtom ]\n")
            f.write(f"{last_atom_id}\n")
    except Exception as e:
        msg = f"Some error occurred while reading the foo.ndx {ndxpath}"
        raise RuntimeError(msg) from e

    # start preparing centered trajectory
    xtcwhole = os.path.join(system_path, "whole.xtc")
    print("Make molecules whole in the trajectory")
    with contextlib.suppress(FileNotFoundError):
        os.remove(xtcwhole)
    try:
        echo_proc = b"System\n"
        subprocess.run(
            [
                "gmx",
                "trjconv",
                "-f",
                trj_name,
                "-s",
                tpr_name,
                "-o",
                xtcwhole,
                "-pbc",
                "mol",
                "-b",
                str(eq_time),
            ],
            input=echo_proc,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        msg = "trjconv for whole.xtc failed"
        raise RuntimeError(msg) from e

    # centering irt methyl-groups
    xtcfoo = os.path.join(system_path, "foo2.xtc")
    with contextlib.suppress(FileNotFoundError):
        os.remove(xtcfoo)
    try:
        echo_input = b"centralAtom\nSystem"
        subprocess.run(
            [
                "gmx",
                "trjconv",
                "-center",
                "-pbc",
                "mol",
                "-n",
                ndxpath,
                "-f",
                xtcwhole,
                "-s",
                tpr_name,
                "-o",
                xtcfoo,
            ],
            input=echo_input,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        msg = f"trjconv for center failed: {e}"
        raise RuntimeError(msg) from e

    try:
        os.remove(ndxpath)
        os.remove(xtcwhole)
    except OSError as e:
        msg = f"Failed to remove temporary files: {e}"
        raise RuntimeError(msg) from e

    # Center around the center of mass of all the g_3 carbons
    try:
        echo_input = f"a {g3_atom}\nq\n".encode()
        subprocess.run(["gmx", "make_ndx", "-f", tpr_name, "-o", ndxpath], input=echo_input, check=True)
        echo_input = f"{g3_atom}\nSystem".encode()
        subprocess.run(
            [
                "gmx",
                "trjconv",
                "-center",
                "-pbc",
                "mol",
                "-n",
                ndxpath,
                "-f",
                xtcfoo,
                "-s",
                tpr_name,
                "-o",
                xtccentered,
            ],
            input=echo_input,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        msg = "Failed during centering on g3 carbons."
        raise RuntimeError(msg) from e

    try:
        os.remove(xtcfoo)
        os.remove(ndxpath)
    except OSError as e:
        msg = f"A error occurred during removing temporary files {ndxpath} & {xtcfoo}."
        raise RuntimeError(msg) from e

    return xtccentered


class NumpyArrayEncoder(CompactJSONEncoder):
    """Encoder for 2xN numpy arrays to be used with json.dump."""

    def encode(self, o) -> str:
        """Encode numpy arrays as lists."""
        if isinstance(o, np.ndarray):
            return CompactJSONEncoder.encode(self, o.tolist())
        return CompactJSONEncoder.encode(self, o)


class FormFactorPlanar(ProfilePlanarBase):
    """Form factor of a planar system based on the linear electron density profile."""

    def __init__(
        self,
        atomgroup: mda.AtomGroup,
        unwrap: bool = True,
        dim: int = 2,
        zmin: float | None = None,
        zmax: float | None = None,
        bin_width: float = 1,
        refgroup: mda.AtomGroup | None = None,
        pack: bool = True,
        output: str = "form_factor.dat",
        concfreq: int = 0,
        jitter: float = 0.0,
    ) -> None:
        self._locals = locals()
        super().__init__(
            atomgroup=atomgroup,
            unwrap=unwrap,
            pack=pack,
            jitter=jitter,
            concfreq=concfreq,
            dim=dim,
            zmin=zmin,
            zmax=zmax,
            bin_width=bin_width,
            refgroup=refgroup,
            sym=False,
            sym_odd=False,
            grouping="atoms",
            bin_method="com",
            output=output,
            weighting_function=density_weights,
            weighting_function_kwargs={"dens": "electron"},
            normalization="volume",
        )

        self.results.scattering_vectors = np.linspace(0, 1, 1000)

    def _single_frame(self) -> float:
        super()._single_frame()

        bin_pos = self._obs.bin_pos - self._obs.box_center[self.dim]

        # Define bulk region as the first 3.3 Å (two water layers) from the box edges
        # `bin_pos[-1]` is the bin corresponding to the box edge
        bulk_mask = np.abs(bin_pos) > bin_pos[-1] - 3.3
        bulk = self._obs.profile[bulk_mask].mean()

        delta = (self._obs.profile - bulk)[:, np.newaxis]
        angles = self.results.scattering_vectors * bin_pos[:, np.newaxis]

        # Here delta is e/A^3 and _bin_width in A,
        # so the result is 100 times lower than if we calculate in nm
        self._obs.ff_real = np.sum(delta * np.cos(angles) * self._bin_width, axis=0)
        self._obs.ff_imag = np.sum(delta * np.sin(angles) * self._bin_width, axis=0)

        # This value at q=0 will be used for a correlation analysis and error estimate
        return np.sqrt(self._obs.ff_real[0] ** 2 + self._obs.ff_imag[0] ** 2)

    def _conclude(self) -> None:
        super()._conclude()

        self.results.form_factor = np.sqrt(
            self.means.ff_real**2 + self.means.ff_imag**2,
        )

        # error from error propagation of the form factor
        self.results.dform_factor = np.sqrt(
            (self.sems.ff_real * self.means.ff_real / self.results.form_factor) ** 2
            + (self.sems.ff_imag * self.means.ff_imag / self.results.form_factor) ** 2,
        )

    def save(self) -> None:
        """
        Save performing unit conversion from Å to nm.

        NOTE: see comments in _single_frame
        """
        output = np.vstack(
            [
                self.results.scattering_vectors,
                self.results.form_factor * 1e2,
                self.results.dform_factor * 1e2,
            ],
        ).T
        with open(self.output, "w") as f:
            json.dump(output, f, cls=NumpyArrayEncoder)


class DensityPlanar(maicos.DensityPlanar):
    """Density profiler for planar system."""

    def save(self) -> None:
        """Save performing unit conversion from Å to nm and e/Å^3 to e/nm^3"""
        outdata = np.vstack(
            [
                self.results.bin_pos / 10,
                self.results.profile * 1e3,
                self.results.dprofile * 1e3,
            ],
        ).T
        with open(self.output, "w") as f:
            json.dump(outdata, f, cls=NumpyArrayEncoder)


class DielectricPlanar(maicos.DielectricPlanar):
    """Dielectric profile for planar system."""

    def save(self) -> None:
        """Save performing unit conversion from Å to nm for the distance."""
        outdata_perp = np.vstack(
            [
                self.results.bin_pos / 10,  # Convert from Å to nm
                self.results.eps_perp,
                self.results.deps_perp,
            ],
        ).T

        with open(f"{self.output_prefix}_perp.json", "w") as f:
            json.dump(outdata_perp, f, cls=NumpyArrayEncoder)

        outdata_par = np.vstack(
            [
                self.results.bin_pos / 10,  # Convert from Å to nm
                self.results.eps_par,
                self.results.deps_par,
            ],
        ).T

        with open(f"{self.output_prefix}_par.json", "w") as f:
            json.dump(outdata_par, f, cls=NumpyArrayEncoder)


class DiporderPlanar(maicos.DiporderPlanar):
    """Dipole order parameter profile for planar system."""

    def save(self) -> None:
        """Save performing unit conversion from Å to nm for the distance."""
        outdata = np.vstack(
            [
                self.results.bin_pos / 10,  # Convert from Å to nm
                self.results.profile,
                self.results.dprofile,
            ],
        ).T
        with open(self.output, "w") as f:
            json.dump(outdata, f, cls=NumpyArrayEncoder)

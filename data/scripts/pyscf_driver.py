#!/usr/bin/env python3
"""PySCF driver for Schrödinger's Sandbox."""

import json
import os
import sys
import time
import traceback

import numpy as np


def _solvent_dielectric(name):
    """Return dielectric constant for common solvents."""
    solvents = {
        "water": 78.36,
        "dmso": 46.7,
        "thf": 7.58,
        "dcm": 8.93,
        "acetonitrile": 35.7,
        "methanol": 32.7,
        "ethanol": 24.3,
        "acetone": 20.7,
        "toluene": 2.38,
        "hexane": 1.88,
        "chloroform": 4.81,
        "dmf": 36.7,
        "benzene": 2.27,
    }
    return solvents.get(name.lower(), 78.36)


def _element_symbol(z):
    from pyscf.data import elements

    if isinstance(elements.ELEMENTS, dict):
        symbol = elements.ELEMENTS.get(z)
    else:
        symbol = elements.ELEMENTS[z]
    if not symbol:
        raise ValueError(f"Unsupported atomic number: {z}")
    return symbol


def _sqrtm_symmetric(matrix):
    vals, vecs = np.linalg.eigh(matrix)
    vals = np.clip(vals, 1.0e-14, None)
    return vecs @ np.diag(np.sqrt(vals)) @ vecs.T


def _mulliken_charges(mol, dm, overlap):
    if dm.ndim == 3:
        dm_total = dm[0] + dm[1]
    else:
        dm_total = dm
    gross = np.einsum("ij,ji->i", dm_total, overlap)
    charges = []
    aoslices = mol.aoslice_by_atom()
    atom_charges = mol.atom_charges()
    for atom_idx, (_, _, p0, p1) in enumerate(aoslices):
        charges.append(float(atom_charges[atom_idx] - np.sum(gross[p0:p1])))
    return charges


def _lowdin_charges(mol, dm, overlap):
    if dm.ndim == 3:
        dm_total = dm[0] + dm[1]
    else:
        dm_total = dm
    s_half = _sqrtm_symmetric(overlap)
    dm_lowdin = s_half @ dm_total @ s_half
    populations = np.diag(dm_lowdin)
    charges = []
    aoslices = mol.aoslice_by_atom()
    atom_charges = mol.atom_charges()
    for atom_idx, (_, _, p0, p1) in enumerate(aoslices):
        charges.append(float(atom_charges[atom_idx] - np.sum(populations[p0:p1])))
    return charges


def _mayer_bond_orders(mol, dm, overlap):
    if dm.ndim == 3:
        dm_total = dm[0] + dm[1]
    else:
        dm_total = dm

    ps = dm_total @ overlap
    natm = mol.natm
    bond_orders = np.zeros((natm, natm))
    aoslices = mol.aoslice_by_atom()

    for i in range(natm):
        _, _, i0, i1 = aoslices[i]
        for j in range(i + 1, natm):
            _, _, j0, j1 = aoslices[j]
            block_ij = ps[i0:i1, j0:j1]
            block_ji = ps[j0:j1, i0:i1]
            bo = float(np.sum(block_ij * block_ji.T))
            bond_orders[i, j] = bo
            bond_orders[j, i] = bo

    return bond_orders.tolist()


def _write_progress(progress_path, stage, iteration=0, energy=0.0, message=""):
    with open(progress_path, "w", encoding="utf-8") as pf:
        json.dump(
            {
                "stage": stage,
                "iteration": iteration,
                "energy": energy,
                "message": message,
                "timestamp": time.time(),
            },
            pf,
        )


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 pyscf_driver.py job.json", file=sys.stderr)
        sys.exit(1)

    with open(sys.argv[1], "r", encoding="utf-8") as f:
        job = json.load(f)

    output_dir = job["output_dir"]
    os.makedirs(output_dir, exist_ok=True)
    progress_path = os.path.join(output_dir, "progress.json")

    result = {
        "success": False,
        "error": "",
        "total_energy": 0.0,
        "dipole_moment": [0.0, 0.0, 0.0],
        "mulliken_charges": [],
        "lowdin_charges": [],
        "mayer_bond_orders": [],
        "orbital_energies": [],
        "orbital_occupations": [],
        "scf_history": [],
        "wall_time": 0.0,
    }

    start_time = time.time()

    try:
        from pyscf import cc, dft, gto, mp, scf
        from pyscf.tools import cubegen
        from pyscf.tools import molden as molden_tools

        _write_progress(progress_path, "building_molecule")

        atom_list = []
        for atomic_number, coords in job["geometry"]:
            symbol = _element_symbol(int(atomic_number))
            x_ang, y_ang, z_ang = [float(c) * 0.529177 for c in coords]
            atom_list.append([symbol, (x_ang, y_ang, z_ang)])

        mol = gto.Mole()
        mol.atom = atom_list
        mol.basis = job["basis"]
        mol.charge = int(job.get("charge", 0))
        mol.spin = int(job.get("multiplicity", 1)) - 1
        mol.unit = "Angstrom"
        mol.verbose = 4
        mol.output = os.path.join(output_dir, "pyscf.log")
        mol.build()

        _write_progress(progress_path, "starting_scf")

        method = job["method"].lower()
        is_open_shell = int(job.get("multiplicity", 1)) > 1
        if method in ("hf", "rhf"):
            mf = scf.UHF(mol) if is_open_shell else scf.RHF(mol)
        elif method == "uhf":
            mf = scf.UHF(mol)
        elif method in ("b3lyp", "pbe", "pbe0", "tpss", "m06-2x", "m062x"):
            xc_map = {
                "b3lyp": "b3lyp",
                "pbe": "pbe",
                "pbe0": "pbe0",
                "tpss": "tpss",
                "m06-2x": "m06-2x",
                "m062x": "m06-2x",
            }
            mf = dft.UKS(mol) if is_open_shell else dft.RKS(mol)
            mf.xc = xc_map[method]
        elif method == "mp2":
            mf = scf.UHF(mol) if is_open_shell else scf.RHF(mol)
        elif method == "ccsd":
            if is_open_shell:
                raise ValueError("CCSD driver currently supports closed-shell systems only")
            mf = scf.RHF(mol)
        else:
            raise ValueError(f"Unknown method: {method}")

        mf.max_cycle = int(job.get("max_scf_cycles", 200))
        mf.conv_tol = float(job.get("scf_convergence", 1.0e-8))

        solvent = job.get("solvent", "")
        if solvent:
            mf = mf.ddCOSMO()
            mf.with_solvent.eps = _solvent_dielectric(solvent)

        scf_iter = [0]

        def scf_callback(envs):
            scf_iter[0] += 1
            energy = float(envs.get("e_tot", 0.0))
            result["scf_history"].append(
                {
                    "iteration": scf_iter[0],
                    "energy": energy,
                }
            )
            _write_progress(
                progress_path,
                "scf",
                scf_iter[0],
                energy,
                f"SCF iteration {scf_iter[0]}",
            )

        mf.callback = scf_callback
        scf_energy = mf.kernel()
        result["total_energy"] = float(scf_energy)

        if not getattr(mf, "converged", False):
            result["error"] = "SCF did not converge"

        post_hf_energy = None
        if method == "mp2":
            _write_progress(progress_path, "post_hf", energy=result["total_energy"], message="MP2")
            post_hf = mp.UMP2(mf) if is_open_shell else mp.MP2(mf)
            corr_energy = post_hf.kernel()[0]
            post_hf_energy = float(mf.e_tot + corr_energy)
            result["total_energy"] = post_hf_energy
        elif method == "ccsd":
            _write_progress(progress_path, "post_hf", energy=result["total_energy"], message="CCSD")
            post_hf = cc.CCSD(mf)
            corr_energy = post_hf.kernel()[0]
            post_hf_energy = float(mf.e_tot + corr_energy)
            result["total_energy"] = post_hf_energy

        _write_progress(progress_path, "computing_properties", energy=result["total_energy"])

        properties = set(job.get("properties", []))
        dm = mf.make_rdm1()
        overlap = mf.get_ovlp()

        if "dipole" in properties:
            dip = mf.dip_moment(verbose=0)
            result["dipole_moment"] = [float(x) for x in dip]

        if "mulliken" in properties:
            result["mulliken_charges"] = _mulliken_charges(mol, dm, overlap)

        if "lowdin" in properties:
            try:
                result["lowdin_charges"] = _lowdin_charges(mol, dm, overlap)
            except Exception:
                result["lowdin_charges"] = []

        try:
            result["mayer_bond_orders"] = _mayer_bond_orders(mol, dm, overlap)
        except Exception:
            result["mayer_bond_orders"] = []

        mo_energy = np.asarray(mf.mo_energy).reshape(-1)
        mo_occ = np.asarray(mf.mo_occ).reshape(-1)
        result["orbital_energies"] = [float(x) for x in mo_energy]
        result["orbital_occupations"] = [float(x) for x in mo_occ]

        if "molden" in properties:
            _write_progress(progress_path, "writing_molden", energy=result["total_energy"])
            molden_path = os.path.join(output_dir, "result.molden")
            molden_tools.from_scf(mf, molden_path)

        cube_res = int(job.get("cube_resolution", 80))
        nx = ny = nz = cube_res

        if "cube_density" in properties:
            _write_progress(progress_path, "writing_cube", energy=result["total_energy"], message="Electron density")
            cubegen.density(
                mol,
                os.path.join(output_dir, "density.cube"),
                dm,
                nx=nx,
                ny=ny,
                nz=nz,
            )

        occupied = np.where(mo_occ > 0.5)[0]
        homo_idx = int(occupied[-1]) if occupied.size else -1

        if "cube_homo" in properties and homo_idx >= 0:
            _write_progress(progress_path, "writing_cube", energy=result["total_energy"], message="HOMO")
            cubegen.orbital(
                mol,
                os.path.join(output_dir, "homo.cube"),
                mf.mo_coeff[:, homo_idx],
                nx=nx,
                ny=ny,
                nz=nz,
            )

        lumo_idx = homo_idx + 1 if homo_idx >= 0 else -1
        if "cube_lumo" in properties and 0 <= lumo_idx < mo_energy.size:
            _write_progress(progress_path, "writing_cube", energy=result["total_energy"], message="LUMO")
            cubegen.orbital(
                mol,
                os.path.join(output_dir, "lumo.cube"),
                mf.mo_coeff[:, lumo_idx],
                nx=nx,
                ny=ny,
                nz=nz,
            )

        if "cube_esp" in properties:
            _write_progress(progress_path, "writing_cube", energy=result["total_energy"], message="Electrostatic potential")
            cubegen.mep(
                mol,
                os.path.join(output_dir, "esp.cube"),
                dm,
                nx=nx,
                ny=ny,
                nz=nz,
            )

        if job.get("optimize", False):
            _write_progress(progress_path, "optimizing", energy=result["total_energy"])
            try:
                from pyscf.geomopt.geometric_solver import optimize as geom_optimize

                mol_eq = geom_optimize(mf, maxsteps=int(job.get("max_opt_steps", 100)))
                opt_coords = mol_eq.atom_coords(unit="Bohr").tolist()
                opt_atoms = [[int(mol_eq.atom_charge(i)), [float(v) for v in opt_coords[i]]] for i in range(mol_eq.natm)]
                result["optimized_geometry"] = opt_atoms
                result["total_energy"] = float(getattr(mf, "e_tot", result["total_energy"]))
            except ImportError:
                result["error"] = "geomeTRIC not installed, optimization unavailable"
            except Exception as exc:
                result["error"] = f"Optimization failed: {exc}"

        if result["error"]:
            result["success"] = False
        else:
            result["success"] = True

        if post_hf_energy is not None:
            result["total_energy"] = post_hf_energy

    except ImportError as exc:
        result["success"] = False
        result["error"] = f"Missing Python package: {exc}"
    except Exception as exc:
        result["success"] = False
        result["error"] = f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"

    result["wall_time"] = time.time() - start_time

    _write_progress(progress_path, "done", energy=result["total_energy"], message=result["error"])
    result_path = os.path.join(output_dir, "result.json")
    with open(result_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, default=str)


if __name__ == "__main__":
    main()

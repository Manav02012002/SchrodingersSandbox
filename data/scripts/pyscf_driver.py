#!/usr/bin/env python3
"""PySCF driver for Schrödinger's Sandbox."""

import json
import os
import sys
import time
import traceback
import importlib

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


def _is_dft_method(method):
    return method in ("b3lyp", "pbe", "pbe0", "tpss", "m06-2x", "m062x")


def _build_mean_field(mol, method, multiplicity, max_scf_cycles, scf_convergence, solvent):
    from pyscf import dft, scf

    is_open_shell = int(multiplicity) > 1
    if method in ("hf", "rhf"):
        mf = scf.UHF(mol) if is_open_shell else scf.RHF(mol)
    elif method == "uhf":
        mf = scf.UHF(mol)
    elif _is_dft_method(method):
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

    mf.max_cycle = int(max_scf_cycles)
    mf.conv_tol = float(scf_convergence)

    if solvent:
        mf = mf.ddCOSMO()
        mf.with_solvent.eps = _solvent_dielectric(solvent)

    return mf


def _is_linear_molecule(mol):
    coords = np.asarray(mol.atom_coords(unit="Bohr"), dtype=float)
    natom = coords.shape[0]
    if natom <= 2:
        return True
    ref = coords[-1] - coords[0]
    ref_norm = np.linalg.norm(ref)
    if ref_norm < 1.0e-10:
        return False
    ref /= ref_norm
    for i in range(1, natom - 1):
        vec = coords[i] - coords[0]
        if np.linalg.norm(np.cross(ref, vec)) > 1.0e-4:
            return False
    return True


def _trim_vibrational_data(mol, frequencies_cm1, normal_modes=None, ir_intensities=None):
    natom = mol.natm
    total_modes = 3 * natom
    if len(frequencies_cm1) != total_modes:
        return frequencies_cm1, normal_modes, ir_intensities

    n_remove = 5 if _is_linear_molecule(mol) else 6
    trimmed_freqs = list(frequencies_cm1[n_remove:])
    trimmed_modes = None if normal_modes is None else list(normal_modes[n_remove:])
    trimmed_intensities = None if ir_intensities is None else list(ir_intensities[n_remove:])
    return trimmed_freqs, trimmed_modes, trimmed_intensities


def _manual_frequency_analysis(mol, hessian_matrix):
    natom = mol.natm
    masses = np.asarray(mol.atom_mass_list(), dtype=float)
    mass_vec = np.repeat(masses, 3)

    h_flat = np.asarray(hessian_matrix, dtype=float).reshape(3 * natom, 3 * natom)
    mass_mat = np.sqrt(np.outer(mass_vec, mass_vec))
    h_mw = h_flat / mass_mat

    eigenvalues, eigenvectors = np.linalg.eigh(h_mw)

    hartree_to_j = 4.3597447222071e-18
    amu_to_kg = 1.66053906660e-27
    bohr_to_m = 5.29177210903e-11
    c_si = 299792458.0
    freq_factor = np.sqrt(hartree_to_j / (amu_to_kg * bohr_to_m * bohr_to_m)) / (2.0 * np.pi * c_si * 100.0)

    frequencies_cm1 = []
    normal_modes = []
    for mode_idx, eigenvalue in enumerate(eigenvalues):
        freq = np.sqrt(abs(eigenvalue)) * freq_factor
        frequencies_cm1.append(float(freq if eigenvalue >= 0.0 else -freq))

        raw_mode = eigenvectors[:, mode_idx] / np.sqrt(np.clip(mass_vec, 1.0e-12, None))
        norm = np.linalg.norm(raw_mode)
        if norm > 1.0e-12:
            raw_mode = raw_mode / norm
        normal_modes.append(raw_mode.astype(float).tolist())

    return frequencies_cm1, normal_modes


def _extract_ir_intensities(mf):
    module_names = []
    class_name = mf.__class__.__name__.lower()
    if "uks" in class_name:
        module_names.extend(["pyscf.prop.infrared.uks", "pyscf.prop.infrared.rks"])
    elif "rks" in class_name:
        module_names.extend(["pyscf.prop.infrared.rks"])
    elif "uhf" in class_name:
        module_names.extend(["pyscf.prop.infrared.uhf", "pyscf.prop.infrared.rhf"])
    else:
        module_names.extend(["pyscf.prop.infrared.rhf"])
    module_names.extend(["pyscf.prop.infrared"])

    for module_name in module_names:
        try:
            module = importlib.import_module(module_name)
        except Exception:
            continue

        infrared_cls = getattr(module, "Infrared", None)
        if infrared_cls is None:
            continue
        try:
            ir_obj = infrared_cls(mf)
            kernel_result = ir_obj.kernel()
            for attr_name in ("ir_intensity", "ir_intensities", "intensities"):
                intensities = getattr(ir_obj, attr_name, None)
                if intensities is not None:
                    return np.asarray(intensities, dtype=float).reshape(-1).tolist()
            if isinstance(kernel_result, dict):
                for key in ("ir_intensity", "ir_intensities", "intensities"):
                    if key in kernel_result:
                        return np.asarray(kernel_result[key], dtype=float).reshape(-1).tolist()
        except Exception:
            continue

    return None


def _run_frequency_analysis(mol, mf, method, progress_path, energy):
    if method not in ("hf", "rhf") and not _is_dft_method(method):
        raise ValueError("Frequency calculations are currently supported for HF and DFT methods only")

    _write_progress(progress_path, "hessian", energy=energy, message="Computing Hessian")
    hessian_matrix = mf.Hessian().kernel()

    frequencies_cm1 = None
    normal_modes = None

    try:
        from pyscf.hessian import thermo

        analysis = thermo.harmonic_analysis(mol, hessian_matrix)
        if isinstance(analysis, dict):
            for key in ("freq_wavenumber", "freq_wavenumber_cm1", "freq_cm1"):
                if key in analysis:
                    frequencies_cm1 = np.asarray(analysis[key], dtype=float).reshape(-1).tolist()
                    break
            raw_modes = None
            for key in ("norm_mode", "normal_mode", "modes"):
                if key in analysis:
                    raw_modes = np.asarray(analysis[key], dtype=float)
                    break
            if raw_modes is not None:
                if raw_modes.ndim == 3:
                    normal_modes = [raw_modes[i].reshape(-1).astype(float).tolist() for i in range(raw_modes.shape[0])]
                elif raw_modes.ndim == 2:
                    normal_modes = [raw_modes[:, i].reshape(-1).astype(float).tolist() for i in range(raw_modes.shape[1])]
    except Exception:
        frequencies_cm1 = None
        normal_modes = None

    if frequencies_cm1 is None:
        frequencies_cm1, normal_modes = _manual_frequency_analysis(mol, hessian_matrix)

    ir_intensities = _extract_ir_intensities(mf)
    frequencies_cm1, normal_modes, ir_intensities = _trim_vibrational_data(
        mol, frequencies_cm1, normal_modes, ir_intensities
    )
    if ir_intensities is None or len(ir_intensities) != len(frequencies_cm1):
        ir_intensities = [1.0] * len(frequencies_cm1)

    return frequencies_cm1, ir_intensities, normal_modes


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
        "frequencies_cm1": [],
        "ir_intensities": [],
        "normal_modes": [],
        "scf_history": [],
        "wall_time": 0.0,
    }

    start_time = time.time()

    try:
        from pyscf import cc, gto, mp
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
        solvent = job.get("solvent", "")
        mf = _build_mean_field(
            mol,
            method,
            int(job.get("multiplicity", 1)),
            int(job.get("max_scf_cycles", 200)),
            float(job.get("scf_convergence", 1.0e-8)),
            solvent,
        )

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

        optimized_mol = None
        if job.get("optimize", False):
            _write_progress(progress_path, "optimizing", energy=result["total_energy"])
            try:
                from pyscf.geomopt.geometric_solver import optimize as geom_optimize

                mol_eq = geom_optimize(mf, maxsteps=int(job.get("max_opt_steps", 100)))
                optimized_mol = mol_eq
                opt_coords = mol_eq.atom_coords(unit="Bohr").tolist()
                opt_atoms = [[int(mol_eq.atom_charge(i)), [float(v) for v in opt_coords[i]]] for i in range(mol_eq.natm)]
                result["optimized_geometry"] = opt_atoms
                result["total_energy"] = float(getattr(mf, "e_tot", result["total_energy"]))
            except ImportError:
                result["error"] = "geomeTRIC not installed, optimization unavailable"
            except Exception as exc:
                result["error"] = f"Optimization failed: {exc}"

        if job.get("frequencies", False):
            freq_mol = optimized_mol if optimized_mol is not None else mol
            freq_mf = mf
            if optimized_mol is not None:
                _write_progress(progress_path, "rebuilding_scf", energy=result["total_energy"], message="SCF on optimized geometry")
                freq_mf = _build_mean_field(
                    freq_mol,
                    method,
                    int(job.get("multiplicity", 1)),
                    int(job.get("max_scf_cycles", 200)),
                    float(job.get("scf_convergence", 1.0e-8)),
                    solvent,
                )
                freq_mf.kernel()
                if not getattr(freq_mf, "converged", False):
                    raise RuntimeError("SCF on optimized geometry did not converge for frequency analysis")

            frequencies_cm1, ir_intensities, normal_modes = _run_frequency_analysis(
                freq_mol,
                freq_mf,
                method,
                progress_path,
                result["total_energy"],
            )
            result["frequencies_cm1"] = [float(v) for v in frequencies_cm1]
            result["ir_intensities"] = [float(v) for v in ir_intensities]
            result["normal_modes"] = [[float(v) for v in mode] for mode in normal_modes]

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

#!/usr/bin/env python3
"""PES scan driver for Schrödinger's Sandbox."""

import copy
import json
import os
import sys
import time
import traceback

import numpy as np


def _write_progress(progress_path, stage, step=0, total=0, energy=0.0, message=""):
    with open(progress_path, "w", encoding="utf-8") as pf:
        json.dump(
            {
                "stage": stage,
                "step": step,
                "total": total,
                "energy": energy,
                "message": message,
                "timestamp": time.time(),
            },
            pf,
        )


def _element_symbol(z):
    from pyscf.data import elements

    if isinstance(elements.ELEMENTS, dict):
        symbol = elements.ELEMENTS.get(z)
    else:
        symbol = elements.ELEMENTS[z]
    if not symbol:
        raise ValueError(f"Unsupported atomic number: {z}")
    return symbol


def _build_method(mol, job):
    from pyscf import dft, scf

    method = str(job["method"]).lower()
    multiplicity = int(job.get("multiplicity", 1))
    open_shell = multiplicity > 1
    if method in ("hf", "rhf"):
        mf = scf.UHF(mol) if open_shell else scf.RHF(mol)
    elif method == "uhf":
        mf = scf.UHF(mol)
    elif method in ("b3lyp", "pbe", "pbe0"):
        mf = dft.UKS(mol) if open_shell else dft.RKS(mol)
        mf.xc = method
    else:
        mf = scf.UHF(mol) if open_shell else scf.RHF(mol)

    mf.verbose = 0
    mf.max_cycle = int(job.get("max_scf_cycles", 200))
    mf.conv_tol = float(job.get("scf_convergence", 1.0e-8))
    return mf


def _constraint_line(coord, value):
    coord_type = coord["type"]
    atoms = coord["atoms"]
    if coord_type == "distance":
        i, j = atoms
        return f"$set\ndistance {i + 1} {j + 1} {float(value):.10f}\n$end"
    if coord_type == "angle":
        i, j, k = atoms
        return f"$set\nangle {i + 1} {j + 1} {k + 1} {float(value):.10f}\n$end"
    if coord_type == "dihedral":
        i, j, k, l = atoms
        return f"$set\ndihedral {i + 1} {j + 1} {k + 1} {l + 1} {float(value):.10f}\n$end"
    raise ValueError(f"Unsupported scan coordinate type: {coord_type}")


def _geometry_to_result(mol):
    coords = np.asarray(mol.atom_coords(unit="Bohr"), dtype=float)
    geometry = []
    for i in range(mol.natm):
        geometry.append([int(mol.atom_charge(i)), [float(v) for v in coords[i]]])
    return geometry


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 pes_scan_driver.py job.json", file=sys.stderr)
        sys.exit(1)

    with open(sys.argv[1], "r", encoding="utf-8") as f:
        job = json.load(f)

    output_dir = job["output_dir"]
    os.makedirs(output_dir, exist_ok=True)
    progress_path = os.path.join(output_dir, "progress.json")

    result = {
        "success": False,
        "error": "",
        "scan_type": "",
        "energies": [],
        "geometries": [],
        "coordinate_values_1": [],
        "coordinate_values_2": [],
        "wall_time": 0.0,
    }

    start_time = time.time()

    def write_result_snapshot():
        with open(os.path.join(output_dir, "result.json"), "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2, default=str)

    try:
        from pyscf import gto
        from pyscf.geomopt.geometric_solver import optimize as geom_optimize

        scan = job["scan"]
        scan_type = str(scan["type"]).lower()
        result["scan_type"] = scan_type

        coord1 = scan["coordinate_1"]
        values_1 = np.linspace(float(coord1["start"]), float(coord1["end"]), int(coord1["steps"]))

        if scan_type == "2d":
            coord2 = scan["coordinate_2"]
            values_2 = np.linspace(float(coord2["start"]), float(coord2["end"]), int(coord2["steps"]))
        else:
            coord2 = None
            values_2 = [None]

        total_points = len(values_1) * len(values_2)
        point_idx = 0

        bohr_to_ang = 0.529177
        atom_list = []
        for atomic_number, coords in job["geometry"]:
            symbol = _element_symbol(int(atomic_number))
            x_ang, y_ang, z_ang = [float(c) * bohr_to_ang for c in coords]
            atom_list.append([symbol, (x_ang, y_ang, z_ang)])

        prev_geometry = copy.deepcopy(atom_list)

        for value_2 in values_2:
            for value_1 in values_1:
                point_idx += 1
                _write_progress(
                    progress_path,
                    "scanning",
                    point_idx,
                    total_points,
                    message=f"Point {point_idx}/{total_points}",
                )

                mol = gto.Mole()
                mol.atom = copy.deepcopy(prev_geometry)
                mol.basis = job["basis"]
                mol.charge = int(job.get("charge", 0))
                mol.spin = int(job.get("multiplicity", 1)) - 1
                mol.unit = "Angstrom"
                mol.verbose = 0
                mol.build()

                constraints_lines = [_constraint_line(coord1, value_1)]
                if coord2 is not None and value_2 is not None:
                    constraints_lines.append(_constraint_line(coord2, value_2))

                constraints_path = os.path.join(output_dir, f"constraints_{point_idx}.txt")
                with open(constraints_path, "w", encoding="utf-8") as f:
                    f.write("\n".join(constraints_lines) + "\n")

                mf = _build_method(mol, job)

                try:
                    mol_opt = geom_optimize(
                        mf,
                        maxsteps=int(job.get("max_opt_steps", 50)),
                        constraints=constraints_path,
                    )
                    optimized_mf = _build_method(mol_opt, job)
                    optimized_mf.kernel()
                    if not getattr(optimized_mf, "converged", False):
                        raise RuntimeError("SCF on optimized scan geometry did not converge")
                    energy = float(optimized_mf.e_tot)
                    final_mol = mol_opt

                    opt_coords = np.asarray(final_mol.atom_coords(unit="Bohr"), dtype=float)
                    prev_geometry = []
                    for atom_idx in range(final_mol.natm):
                        x, y, z = opt_coords[atom_idx] * bohr_to_ang
                        prev_geometry.append([final_mol.atom_symbol(atom_idx), (float(x), float(y), float(z))])
                except Exception:
                    mf.kernel()
                    energy = float(mf.e_tot)
                    final_mol = mol

                result["energies"].append(energy)
                result["coordinate_values_1"].append(float(value_1))
                if value_2 is not None:
                    result["coordinate_values_2"].append(float(value_2))
                result["geometries"].append(_geometry_to_result(final_mol))
                write_result_snapshot()

        result["success"] = True
    except Exception as exc:
        result["success"] = False
        result["error"] = f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"

    result["wall_time"] = time.time() - start_time
    _write_progress(
        progress_path,
        "done",
        total=len(result["energies"]),
        energy=result["energies"][-1] if result["energies"] else 0.0,
        message=result["error"],
    )

    write_result_snapshot()


if __name__ == "__main__":
    main()

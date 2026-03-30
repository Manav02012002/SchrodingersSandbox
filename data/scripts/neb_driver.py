#!/usr/bin/env python3
"""NEB driver for Schrödinger's Sandbox."""

import json
import os
import sys
import time
import traceback

import numpy as np


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 neb_driver.py job.json", file=sys.stderr)
        sys.exit(1)

    with open(sys.argv[1], "r", encoding="utf-8") as f:
        job = json.load(f)

    output_dir = job["output_dir"]
    os.makedirs(output_dir, exist_ok=True)

    progress_path = os.path.join(output_dir, "progress.json")

    def write_progress(stage, step=0, message=""):
        with open(progress_path, "w", encoding="utf-8") as pf:
            json.dump(
                {"stage": stage, "step": step, "message": message, "timestamp": time.time()},
                pf,
            )

    result = {
        "success": False,
        "error": "",
        "path_energies": [],
        "path_geometries": [],
        "ts_index": -1,
        "ts_energy": 0.0,
        "forward_barrier": 0.0,
        "reverse_barrier": 0.0,
        "converged": False,
        "wall_time": 0.0,
    }

    start_time = time.time()

    try:
        from pyscf import dft, gto, scf
        from pyscf.data import elements

        bohr_to_ang = 0.529177

        def geom_to_pyscf(geom):
            atom_list = []
            for atomic_number, coords in geom:
                symbol = elements.ELEMENTS[atomic_number]
                x, y, z = [float(c) * bohr_to_ang for c in coords]
                atom_list.append([symbol, (x, y, z)])
            return atom_list

        reactant_atoms = geom_to_pyscf(job["reactant"])
        product_atoms = geom_to_pyscf(job["product"])
        num_images = int(job.get("num_images", 9))

        write_progress("interpolating", message="Generating initial path")

        reactant_coords = np.array([coords for _, coords in job["reactant"]], dtype=float)
        product_coords = np.array([coords for _, coords in job["product"]], dtype=float)

        images_coords = []
        for i in range(num_images):
            t = i / max(1, num_images - 1)
            images_coords.append(reactant_coords * (1.0 - t) + product_coords * t)

        try:
            from ase import Atoms
            from ase.calculators.calculator import Calculator, all_changes
            from ase.mep import NEB
            from ase.optimize import BFGS

            class PySCFCalculator(Calculator):
                implemented_properties = ["energy", "forces"]

                def __init__(self, method, basis, charge, spin, **kwargs):
                    super().__init__(**kwargs)
                    self.method = method
                    self.basis = basis
                    self.charge = charge
                    self.spin = spin

                def calculate(self, atoms=None, properties=None, system_changes=all_changes):
                    super().calculate(atoms, properties or ["energy"], system_changes)

                    mol = gto.Mole()
                    atom_list = []
                    for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
                        atom_list.append([symbol, tuple(pos)])
                    mol.atom = atom_list
                    mol.basis = self.basis
                    mol.charge = self.charge
                    mol.spin = self.spin
                    mol.unit = "Angstrom"
                    mol.verbose = 0
                    mol.build()

                    method = str(self.method).lower()
                    if method in ("hf", "rhf"):
                        mf = scf.RHF(mol)
                    elif method in ("b3lyp", "pbe", "pbe0"):
                        mf = dft.RKS(mol)
                        mf.xc = method
                    else:
                        mf = scf.RHF(mol)

                    mf.verbose = 0
                    mf.kernel()
                    self.results["energy"] = float(mf.e_tot) * 27.2114
                    grad = mf.nuc_grad_method().kernel()
                    self.results["forces"] = -grad * 27.2114 / 0.529177

            atomic_numbers = [atomic_number for atomic_number, _ in job["reactant"]]
            ase_images = []
            for i in range(num_images):
                coords_ang = images_coords[i] * bohr_to_ang
                atoms = Atoms(numbers=atomic_numbers, positions=coords_ang)
                atoms.calc = PySCFCalculator(
                    method=job["method"],
                    basis=job["basis"],
                    charge=int(job.get("charge", 0)),
                    spin=int(job.get("multiplicity", 1)) - 1,
                )
                ase_images.append(atoms)

            write_progress("neb", message="Running NEB optimization")
            neb = NEB(ase_images, climb=True)
            optimizer = BFGS(neb, trajectory=os.path.join(output_dir, "neb.traj"))
            optimizer.run(fmax=0.05, steps=int(job.get("max_neb_steps", 50)))

            for atoms in ase_images:
                energy_hartree = atoms.get_potential_energy() / 27.2114
                result["path_energies"].append(float(energy_hartree))
                coords_bohr = (atoms.get_positions() / bohr_to_ang).tolist()
                geom = [[int(z), [float(v) for v in coords]] for z, coords in zip(atomic_numbers, coords_bohr)]
                result["path_geometries"].append(geom)

            result["converged"] = bool(getattr(optimizer, "converged", lambda: False)())
        except ImportError:
            write_progress("single_points", message="ASE not available, computing single-point energies along linear path")

            for i in range(num_images):
                write_progress("single_points", i, message=f"Image {i + 1}/{num_images}")

                mol = gto.Mole()
                coords_ang = images_coords[i] * bohr_to_ang
                atom_list = []
                for atom_index, (atomic_number, _) in enumerate(job["reactant"]):
                    symbol = elements.ELEMENTS[atomic_number]
                    atom_list.append([symbol, tuple(coords_ang[atom_index])])
                mol.atom = atom_list
                mol.basis = job["basis"]
                mol.charge = int(job.get("charge", 0))
                mol.spin = int(job.get("multiplicity", 1)) - 1
                mol.unit = "Angstrom"
                mol.verbose = 0
                mol.build()

                method = str(job["method"]).lower()
                if method in ("hf", "rhf"):
                    mf = scf.RHF(mol)
                elif method in ("b3lyp", "pbe", "pbe0"):
                    mf = dft.RKS(mol)
                    mf.xc = method
                else:
                    mf = scf.RHF(mol)
                mf.verbose = 0
                mf.kernel()

                result["path_energies"].append(float(mf.e_tot))
                result["path_geometries"].append(
                    [[int(atomic_number), [float(v) for v in images_coords[i][atom_index]]]
                     for atom_index, (atomic_number, _) in enumerate(job["reactant"])]
                )

            result["converged"] = False

        energies = result["path_energies"]
        if len(energies) > 2:
            interior = energies[1:-1]
            ts_idx = int(np.argmax(interior)) + 1
            result["ts_index"] = ts_idx
            result["ts_energy"] = float(energies[ts_idx])
            result["forward_barrier"] = float(energies[ts_idx] - energies[0])
            result["reverse_barrier"] = float(energies[ts_idx] - energies[-1])

        result["success"] = True
    except Exception as exc:
        result["success"] = False
        result["error"] = f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"

    result["wall_time"] = time.time() - start_time
    write_progress("done")

    with open(os.path.join(output_dir, "result.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, default=str)


if __name__ == "__main__":
    main()

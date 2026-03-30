#!/usr/bin/env python3
"""xTB driver for Schrödinger's Sandbox."""

import json
import os
import shutil
import subprocess
import sys
import time
import traceback

import numpy as np

SYMBOLS = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am",
    "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc",
    "Lv", "Ts", "Og",
]


def write_progress(progress_path, stage, message=""):
    with open(progress_path, "w", encoding="utf-8") as pf:
        json.dump({"stage": stage, "message": message, "timestamp": time.time()}, pf)


def _method_name(method):
    if method == "gfn2-xtb":
        return "GFN2-xTB"
    if method == "gfn1-xtb":
        return "GFN1-xTB"
    return "GFN-FF"


def _gfn_level(method):
    if method == "gfn1-xtb":
        return 1
    if method == "gfn-ff":
        return 0
    return 2


def _find_xtb():
    path = shutil.which("xtb")
    if path:
        return path
    for candidate in [
        os.path.expanduser("~/miniconda3/bin/xtb"),
        os.path.expanduser("~/anaconda3/bin/xtb"),
        "/opt/homebrew/bin/xtb",
        "/usr/local/bin/xtb",
    ]:
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None


def _write_xyz(path, geometry):
    bohr_to_ang = 0.529177
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{len(geometry)}\n")
        f.write("sbox job\n")
        for atomic_number, coords in geometry:
            z = int(atomic_number)
            sym = SYMBOLS[z] if 0 < z < len(SYMBOLS) else "X"
            x_ang, y_ang, z_ang = [float(c) * bohr_to_ang for c in coords]
            f.write(f"{sym} {x_ang:.10f} {y_ang:.10f} {z_ang:.10f}\n")


def _read_xyz_frames(path):
    frames = []
    with open(path, "r", encoding="utf-8") as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue
            nat = int(line)
            comment = f.readline()
            atoms = []
            for _ in range(nat):
                atom_line = f.readline()
                if not atom_line:
                    raise RuntimeError(f"Unexpected end of XYZ file: {path}")
                parts = atom_line.split()
                if len(parts) < 4:
                    raise RuntimeError(f"Invalid XYZ line in {path}: {atom_line}")
                sym = parts[0]
                x, y, z = map(float, parts[1:4])
                atoms.append((sym, x, y, z))
            frames.append((comment.strip(), atoms))
    return frames


def _xyz_frames_to_bohr_geometry(frames):
    ang_to_bohr = 1.0 / 0.529177
    converted = []
    for _, atoms in frames:
        geometry = []
        for sym, x, y, z in atoms:
            geometry.append([SYMBOLS.index(sym), [x * ang_to_bohr, y * ang_to_bohr, z * ang_to_bohr]])
        converted.append(geometry)
    return converted


def _write_trajectory_xyz(output_dir, source_path):
    if not os.path.exists(source_path):
        return []
    frames = _read_xyz_frames(source_path)
    traj_path = os.path.join(output_dir, "trajectory.xyz")
    shutil.copyfile(source_path, traj_path)
    return _xyz_frames_to_bohr_geometry(frames)


def _is_xyz_atom_line(line):
    parts = line.split()
    if len(parts) < 4 or parts[0] not in SYMBOLS:
        return False
    try:
        float(parts[1])
        float(parts[2])
        float(parts[3])
    except ValueError:
        return False
    return True


def _extract_embedded_xyz_frames(text):
    lines = text.splitlines()
    frames = []
    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        try:
            nat = int(stripped)
        except ValueError:
            i += 1
            continue

        if nat <= 0 or i + 1 + nat >= len(lines):
            i += 1
            continue

        comment = lines[i + 1].rstrip("\n")
        atom_lines = lines[i + 2 : i + 2 + nat]
        if all(_is_xyz_atom_line(line) for line in atom_lines):
            atoms = []
            for line in atom_lines:
                parts = line.split()
                atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
            frames.append((comment.strip(), atoms))
            i += nat + 2
            continue

        i += 1
    return frames


def _write_trajectory_from_text(output_dir, text):
    frames = _extract_embedded_xyz_frames(text)
    if not frames:
        return []
    traj_path = os.path.join(output_dir, "trajectory.xyz")
    with open(traj_path, "w", encoding="utf-8") as f:
        for comment, atoms in frames:
            f.write(f"{len(atoms)}\n")
            f.write(f"{comment}\n")
            for sym, x, y, z in atoms:
                f.write(f"{sym} {x:.10f} {y:.10f} {z:.10f}\n")
    return _xyz_frames_to_bohr_geometry(frames)


def _run_tblite(job, result):
    from tblite.interface import Calculator

    geometry = job["geometry"]
    method = job["method"].lower()
    numbers = np.array([int(z) for z, _ in geometry], dtype=np.int32)
    positions = np.array([[float(v) for v in coords] for _, coords in geometry], dtype=np.float64)

    calc = Calculator(_method_name(method), numbers, positions)
    calc.set("verbosity", 0)
    if job.get("charge", 0) != 0:
        calc.set("charge", float(job["charge"]))
    if job.get("multiplicity", 1) != 1:
        calc.set("spin", int(job["multiplicity"]) - 1)
    if job.get("solvent", ""):
        calc.set("solvent", str(job["solvent"]))

    res = calc.singlepoint()
    result["total_energy"] = float(res.get("energy"))

    charges = res.get("charges")
    if charges is not None:
        result["mulliken_charges"] = [float(x) for x in charges]

    dipole = res.get("dipole")
    if dipole is not None:
        result["dipole_moment"] = [float(x) for x in dipole]

    orbital_energies = res.get("orbital-energies")
    if orbital_energies is not None:
        result["orbital_energies"] = [float(x) for x in np.asarray(orbital_energies).reshape(-1)]

    occupations = res.get("orbital-occupations")
    if occupations is not None:
        result["orbital_occupations"] = [float(x) for x in np.asarray(occupations).reshape(-1)]

    result["success"] = True


def _parse_xtb_stdout_energy(stdout_text):
    for line in stdout_text.splitlines():
        if "TOTAL ENERGY" not in line:
            continue
        parts = line.replace(":", " ").split()
        for token in reversed(parts):
            try:
                return float(token)
            except ValueError:
                continue
    return None


def _run_xtb_cli(job, output_dir, xyz_path, result):
    xtb_path = _find_xtb()
    if not xtb_path:
        raise RuntimeError(
            "Neither tblite Python package nor xtb command-line tool found. "
            "Install tblite or xtb."
        )

    method = job["method"].lower()
    gfn_level = _gfn_level(method)

    cmd = [xtb_path, xyz_path]
    if gfn_level > 0:
        cmd.extend(["--gfn", str(gfn_level)])
    else:
        cmd.append("--gfnff")
    cmd.extend(["--chrg", str(job.get("charge", 0))])
    cmd.extend(["--uhf", str(int(job.get("multiplicity", 1)) - 1)])
    cmd.extend(["--json"])
    properties = set(job.get("properties", []))
    if "molden" in properties:
        cmd.append("--molden")

    if job.get("optimize", False):
        cmd.append("--opt")

    if job.get("solvent", ""):
        cmd.extend(["--alpb", str(job["solvent"])])

    proc = subprocess.run(
        cmd,
        cwd=output_dir,
        capture_output=True,
        text=True,
        timeout=300,
    )

    if proc.returncode != 0:
        stderr = proc.stderr.strip() or proc.stdout.strip()
        raise RuntimeError(f"xtb failed: {stderr[:500]}")

    xtb_json_path = os.path.join(output_dir, "xtbout.json")
    if os.path.exists(xtb_json_path):
        with open(xtb_json_path, "r", encoding="utf-8") as f:
            xtb_result = json.load(f)

        result["total_energy"] = float(
            xtb_result.get("total energy", xtb_result.get("total_energy", 0.0))
        )

        charges = xtb_result.get("partial charges", xtb_result.get("partial_charges"))
        if charges is not None:
            result["mulliken_charges"] = [float(x) for x in charges]

        dipole = xtb_result.get("dipole")
        if dipole is not None:
            result["dipole_moment"] = [float(x) for x in dipole]

        orbital_block = xtb_result.get("orbital energies", xtb_result.get("orbital_energies"))
        if isinstance(orbital_block, dict):
            alpha = orbital_block.get("alphaorbitals", orbital_block.get("alpha"))
            occ = orbital_block.get("occupations")
            if alpha is not None:
                result["orbital_energies"] = [float(x) for x in alpha]
            if occ is not None:
                result["orbital_occupations"] = [float(x) for x in occ]
        elif orbital_block is not None:
            result["orbital_energies"] = [float(x) for x in orbital_block]
    else:
        energy = _parse_xtb_stdout_energy(proc.stdout)
        if energy is None:
            raise RuntimeError("xtb completed but no JSON or total energy was produced")
        result["total_energy"] = float(energy)

    molden_path = os.path.join(output_dir, "molden.input")
    if os.path.exists(molden_path):
        shutil.copy(molden_path, os.path.join(output_dir, "result.molden"))

    if job.get("optimize", False):
        result["optimization_converged"] = True
        traj_frames = []
        for candidate in ("xtbopt.trj", "xtb.trj"):
            candidate_path = os.path.join(output_dir, candidate)
            if os.path.exists(candidate_path):
                try:
                    traj_frames = _write_trajectory_xyz(output_dir, candidate_path)
                    if traj_frames:
                        break
                except Exception:
                    traj_frames = []
        if not traj_frames:
            for candidate in ("xtbopt.log", "xtbopt.xyzlog"):
                candidate_path = os.path.join(output_dir, candidate)
                if os.path.exists(candidate_path):
                    try:
                        with open(candidate_path, "r", encoding="utf-8", errors="ignore") as f:
                            traj_frames = _write_trajectory_from_text(output_dir, f.read())
                        if traj_frames:
                            break
                    except Exception:
                        traj_frames = []
        if not traj_frames and proc.stdout:
            try:
                traj_frames = _write_trajectory_from_text(output_dir, proc.stdout)
            except Exception:
                traj_frames = []
        opt_xyz = os.path.join(output_dir, "xtbopt.xyz")
        if os.path.exists(opt_xyz):
            final_frames = _xyz_frames_to_bohr_geometry(_read_xyz_frames(opt_xyz))
            if final_frames:
                result["optimized_geometry"] = final_frames[-1]
                if not traj_frames:
                    shutil.copyfile(opt_xyz, os.path.join(output_dir, "trajectory.xyz"))
                    traj_frames = final_frames
        result["opt_steps"] = len(traj_frames)

    result["success"] = True


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 xtb_driver.py job.json", file=sys.stderr)
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
        "optimization_converged": False,
        "dipole_moment": [0.0, 0.0, 0.0],
        "mulliken_charges": [],
        "orbital_energies": [],
        "orbital_occupations": [],
        "wall_time": 0.0,
    }

    start_time = time.time()

    try:
        geometry = job["geometry"]
        xyz_path = os.path.join(output_dir, "input.xyz")
        _write_xyz(xyz_path, geometry)

        write_progress(progress_path, "running_xtb", "Preparing xTB calculation")

        properties = set(job.get("properties", []))
        need_molden = "molden" in properties
        used_tblite = False
        if (not need_molden and not job.get("optimize", False)) or _find_xtb() is None:
            try:
                _run_tblite(job, result)
                used_tblite = True
            except ImportError:
                used_tblite = False
            except Exception:
                used_tblite = False

        if not used_tblite:
            _run_xtb_cli(job, output_dir, xyz_path, result)

    except Exception as exc:
        result["success"] = False
        result["error"] = f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"

    result["wall_time"] = time.time() - start_time
    write_progress(progress_path, "done", message=f"Energy: {result['total_energy']:.6f} Hartree")

    result_path = os.path.join(output_dir, "result.json")
    with open(result_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, default=str)


if __name__ == "__main__":
    main()

#include "ui/molecule_info.h"

#include "core/elements.h"
#include "core/symmetry.h"

#include <imgui.h>

#include <map>
#include <sstream>
#include <string>

namespace sbox::ui {
namespace {

std::string hill_formula(const sbox::chem::MolecularSystem& mol) {
    std::map<std::string, int> counts;
    for (const sbox::chem::Atom& atom : mol.atoms()) {
        counts[sbox::elements::get_element(atom.Z).symbol] += 1;
    }

    std::ostringstream out;
    auto emit = [&](const std::string& symbol) {
        const auto it = counts.find(symbol);
        if (it == counts.end()) {
            return;
        }
        out << it->first;
        if (it->second != 1) {
            out << it->second;
        }
        counts.erase(it);
    };

    emit("C");
    emit("H");
    for (const auto& entry : counts) {
        out << entry.first;
        if (entry.second != 1) {
            out << entry.second;
        }
    }
    return out.str();
}

std::string signature_for_molecule(const sbox::chem::MolecularSystem& mol) {
    std::ostringstream out;
    out << mol.name() << '|' << mol.charge() << '|' << mol.multiplicity() << '|' << mol.num_atoms();
    for (const sbox::chem::Atom& atom : mol.atoms()) {
        out << '|' << atom.Z << ':'
            << atom.position.x() << ','
            << atom.position.y() << ','
            << atom.position.z();
    }
    return out.str();
}

}  // namespace

void draw_molecule_info(const AppState& state, const sbox::chem::MolecularSystem& mol) {
    if (!ImGui::Begin("Molecule Info")) {
        ImGui::End();
        return;
    }

    if (mol.num_atoms() == 0) {
        ImGui::TextUnformatted("No molecule loaded.");
        ImGui::End();
        return;
    }

    static std::string cached_signature;
    static sbox::chem::PointGroup cached_pg = sbox::chem::PointGroup::Unknown;

    const std::string signature = signature_for_molecule(mol);
    if (signature != cached_signature) {
        cached_signature = signature;
        cached_pg = sbox::chem::detect_point_group(mol, 0.3);
    }

    const Eigen::Vector3d com = mol.center_of_mass();

    ImGui::SetWindowFontScale(1.2f);
    ImGui::Text("%s", mol.name().empty() ? "Unnamed Molecule" : mol.name().c_str());
    ImGui::SetWindowFontScale(1.0f);

    const std::string formula = hill_formula(mol);
    ImGui::Text("Formula: %s", formula.c_str());
    ImGui::Text("Atoms: %d", mol.num_atoms());
    ImGui::Text("Bonds: %d", mol.num_bonds());
    ImGui::Text("Charge: %d", mol.charge());
    ImGui::Text("Multiplicity: %d", mol.multiplicity());
    ImGui::Text("Point Group: %s", sbox::chem::point_group_name(cached_pg).c_str());
    ImGui::Text("Center of mass: (%.3f, %.3f, %.3f)", com.x(), com.y(), com.z());

    if (state.mol_has_mo_summary) {
        ImGui::Separator();
        ImGui::Text("Basis functions: %d", state.mol_num_basis);
        ImGui::Text("MOs: %d", state.num_mo);
        ImGui::Text("Total energy: %.6f Ha", state.mol_total_energy_h);
        ImGui::Text("Total energy: %.3f eV", state.mol_total_energy_h * 27.2114);
        if (state.mol_homo_lumo_gap_ev > 0.0) {
            ImGui::Text("HOMO-LUMO gap: %.3f eV", state.mol_homo_lumo_gap_ev);
        }
    }

    ImGui::End();
}

}  // namespace sbox::ui

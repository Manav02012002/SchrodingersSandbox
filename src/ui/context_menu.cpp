#include "ui/context_menu.h"

#include "core/elements.h"

#include <imgui.h>

#include <cstdio>
#include <queue>
#include <string>
#include <utility>
#include <vector>

namespace sbox::ui {

namespace {

std::string atom_label(const sbox::chem::MolecularSystem& mol, int index) {
    return std::string(sbox::elements::get_element(mol.atom(index).Z).symbol) + std::to_string(index + 1);
}

void select_connected_atoms(const sbox::chem::MolecularSystem& mol,
                            int seed_atom,
                            sbox::editor::Selection& selection) {
    if (seed_atom < 0 || seed_atom >= mol.num_atoms()) {
        return;
    }
    std::vector<bool> visited(static_cast<std::size_t>(mol.num_atoms()), false);
    std::queue<int> pending;
    selection.clear();
    pending.push(seed_atom);
    visited[static_cast<std::size_t>(seed_atom)] = true;

    while (!pending.empty()) {
        const int atom = pending.front();
        pending.pop();
        selection.atoms.push_back(atom);
        for (int neighbor : mol.neighbors(atom)) {
            if (!visited[static_cast<std::size_t>(neighbor)]) {
                visited[static_cast<std::size_t>(neighbor)] = true;
                pending.push(neighbor);
            }
        }
    }
}

}  // namespace

void draw_context_menu(ContextMenuState& ctx,
                       sbox::chem::MolecularSystem& mol,
                       sbox::editor::Selection& selection,
                       sbox::editor::CommandStack& commands) {
    if (ctx.show) {
        ImGui::SetNextWindowPos(ctx.position);
        ImGui::OpenPopup("##editor_context");
        ctx.show = false;
    }

    if (!ImGui::BeginPopup("##editor_context")) {
        return;
    }

    if (ctx.clicked_atom >= 0 && ctx.clicked_atom < mol.num_atoms()) {
        ImGui::Text("Atom: %s", atom_label(mol, ctx.clicked_atom).c_str());
        ImGui::Separator();
        if (ImGui::MenuItem("Delete Atom")) {
            commands.execute(std::make_unique<sbox::editor::RemoveAtomCommand>(ctx.clicked_atom), mol);
            selection.clear();
        }
        if (ImGui::MenuItem("Add Hydrogen")) {
            commands.execute(std::make_unique<sbox::editor::AddHydrogensCommand>(ctx.clicked_atom), mol);
        }
        if (ImGui::BeginMenu("Set Element")) {
            const int common[] = {1, 6, 7, 8, 16, 15, 9, 17, 35};
            for (int Z : common) {
                if (ImGui::MenuItem(sbox::elements::get_element(Z).symbol)) {
                    commands.execute(std::make_unique<sbox::editor::SetElementCommand>(ctx.clicked_atom, Z), mol);
                }
            }
            if (ImGui::BeginMenu("Other...")) {
                for (int Z = 1; Z <= 118; ++Z) {
                    if (ImGui::MenuItem(sbox::elements::get_element(Z).symbol)) {
                        commands.execute(std::make_unique<sbox::editor::SetElementCommand>(ctx.clicked_atom, Z), mol);
                    }
                }
                ImGui::EndMenu();
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Set Formal Charge")) {
            for (int charge = -2; charge <= 2; ++charge) {
                char label[8];
                std::snprintf(label, sizeof(label), "%+d", charge);
                if (ImGui::MenuItem(label)) {
                    commands.execute(std::make_unique<sbox::editor::SetChargeCommand>(ctx.clicked_atom, charge), mol);
                }
            }
            ImGui::EndMenu();
        }
        ImGui::Separator();
        if (ImGui::MenuItem("Select Connected")) {
            select_connected_atoms(mol, ctx.clicked_atom, selection);
        }
        if (ImGui::MenuItem("Select Fragment")) {
            select_connected_atoms(mol, ctx.clicked_atom, selection);
        }
    } else if (ctx.clicked_bond >= 0 && ctx.clicked_bond < mol.num_bonds()) {
        const auto& bond = mol.bond(ctx.clicked_bond);
        ImGui::Text("Bond: %s-%s",
                    atom_label(mol, bond.atom_i).c_str(),
                    atom_label(mol, bond.atom_j).c_str());
        ImGui::Separator();
        if (ImGui::MenuItem("Delete Bond")) {
            commands.execute(std::make_unique<sbox::editor::RemoveBondCommand>(ctx.clicked_bond), mol);
            selection.clear();
        }
        if (ImGui::BeginMenu("Set Bond Order")) {
            const std::pair<const char*, sbox::chem::BondOrder> entries[] = {
                {"Single", sbox::chem::BondOrder::Single},
                {"Double", sbox::chem::BondOrder::Double},
                {"Triple", sbox::chem::BondOrder::Triple},
                {"Aromatic", sbox::chem::BondOrder::Aromatic},
            };
            for (const auto& [label, order] : entries) {
                if (ImGui::MenuItem(label)) {
                    commands.execute(std::make_unique<sbox::editor::ChangeBondOrderCommand>(ctx.clicked_bond, order), mol);
                }
            }
            ImGui::EndMenu();
        }
        ImGui::MenuItem("Rotate Around Bond", nullptr, false, false);
    } else {
        ImGui::MenuItem("Paste", nullptr, false, false);
        if (ImGui::MenuItem("Center View")) {
            ctx.center_view_requested = true;
        }
        if (ImGui::MenuItem("Fit to View")) {
            ctx.fit_view_requested = true;
        }
    }

    ImGui::EndPopup();
}

}  // namespace sbox::ui

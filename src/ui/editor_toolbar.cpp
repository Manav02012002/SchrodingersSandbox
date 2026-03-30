#include "ui/editor_toolbar.h"

#include "core/elements.h"
#include "core/valence.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace sbox::ui {

namespace {

constexpr double kBohrToAngstrom = 1.0 / 1.8897259886;

const char* mode_name(EditorState::Mode mode) {
    switch (mode) {
    case EditorState::Mode::Select: return "Select";
    case EditorState::Mode::Draw: return "Draw";
    case EditorState::Mode::Erase: return "Erase";
    case EditorState::Mode::Measure: return "Measure";
    case EditorState::Mode::Fragment: return "Fragment";
    }
    return "Select";
}

const char* bond_order_name(sbox::chem::BondOrder order) {
    switch (order) {
    case sbox::chem::BondOrder::Single: return "Single";
    case sbox::chem::BondOrder::Double: return "Double";
    case sbox::chem::BondOrder::Triple: return "Triple";
    case sbox::chem::BondOrder::Aromatic: return "Aromatic";
    case sbox::chem::BondOrder::Unknown: return "Unknown";
    }
    return "Unknown";
}

bool toolbar_button(const char* label, bool active) {
    if (active) {
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.16f, 0.55f, 0.52f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.20f, 0.65f, 0.60f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.14f, 0.45f, 0.42f, 1.0f));
    }
    const bool pressed = ImGui::Button(label);
    if (active) {
        ImGui::PopStyleColor(3);
    }
    return pressed;
}

std::vector<int> all_atom_indices(const sbox::chem::MolecularSystem& mol) {
    std::vector<int> indices;
    indices.reserve(static_cast<std::size_t>(mol.num_atoms()));
    for (int i = 0; i < mol.num_atoms(); ++i) {
        indices.push_back(i);
    }
    return indices;
}

}  // namespace

EditorState::EditorState()
    : select_mode(std::make_unique<sbox::editor::SelectMode>()),
      draw_mode(std::make_unique<sbox::editor::DrawMode>()),
      erase_mode(std::make_unique<sbox::editor::EraseMode>()),
      measure_mode(std::make_unique<sbox::editor::MeasureMode>()),
      fragment_mode(std::make_unique<sbox::editor::FragmentMode>(&fragment_library)) {
    select_mode->set_context_menu_state(&context_menu);
}

sbox::editor::EditorMode* EditorState::active_mode() {
    switch (current_mode) {
    case Mode::Select: return select_mode.get();
    case Mode::Draw: return draw_mode.get();
    case Mode::Erase: return erase_mode.get();
    case Mode::Measure: return measure_mode.get();
    case Mode::Fragment: return fragment_mode.get();
    }
    return select_mode.get();
}

const sbox::editor::EditorMode* EditorState::active_mode() const {
    switch (current_mode) {
    case Mode::Select: return select_mode.get();
    case Mode::Draw: return draw_mode.get();
    case Mode::Erase: return erase_mode.get();
    case Mode::Measure: return measure_mode.get();
    case Mode::Fragment: return fragment_mode.get();
    }
    return select_mode.get();
}

void draw_editor_toolbar(EditorState& editor, sbox::chem::MolecularSystem& mol) {
    if (!ImGui::Begin("Editor")) {
        ImGui::End();
        return;
    }

    if (toolbar_button("Select", editor.current_mode == EditorState::Mode::Select)) {
        editor.current_mode = EditorState::Mode::Select;
    }
    ImGui::SameLine();
    if (toolbar_button("Draw", editor.current_mode == EditorState::Mode::Draw)) {
        editor.current_mode = EditorState::Mode::Draw;
    }
    ImGui::SameLine();
    if (toolbar_button("Erase", editor.current_mode == EditorState::Mode::Erase)) {
        editor.current_mode = EditorState::Mode::Erase;
    }
    ImGui::SameLine();
    if (toolbar_button("Measure", editor.current_mode == EditorState::Mode::Measure)) {
        editor.current_mode = EditorState::Mode::Measure;
    }
    ImGui::SameLine();
    if (toolbar_button("Fragment", editor.current_mode == EditorState::Mode::Fragment)) {
        editor.current_mode = EditorState::Mode::Fragment;
    }

    ImGui::Separator();

    if (editor.current_mode == EditorState::Mode::Draw && editor.draw_mode) {
        ImGui::TextUnformatted("Elements");
        const int common_elements[] = {1, 6, 7, 8, 16, 15, 9, 17, 35};
        for (int idx = 0; idx < static_cast<int>(std::size(common_elements)); ++idx) {
            const int Z = common_elements[idx];
            if (idx > 0) {
                ImGui::SameLine();
            }
            const bool active = editor.draw_mode->element() == Z;
            if (toolbar_button(sbox::elements::get_element(Z).symbol, active)) {
                editor.draw_mode->set_element(Z);
            }
        }
        if (ImGui::Button("Other...")) {
            ImGui::OpenPopup("editor_periodic_popup");
        }
        if (ImGui::BeginPopup("editor_periodic_popup")) {
            for (int Z = 1; Z <= 118; ++Z) {
                if (ImGui::Selectable(sbox::elements::get_element(Z).symbol, editor.draw_mode->element() == Z)) {
                    editor.draw_mode->set_element(Z);
                    ImGui::CloseCurrentPopup();
                }
                if (Z % 10 != 0) {
                    ImGui::SameLine();
                }
            }
            ImGui::EndPopup();
        }

        ImGui::Spacing();
        ImGui::TextUnformatted("Bond Order");
        const sbox::chem::BondOrder orders[] = {
            sbox::chem::BondOrder::Single,
            sbox::chem::BondOrder::Double,
            sbox::chem::BondOrder::Triple,
        };
        for (int i = 0; i < 3; ++i) {
            if (i > 0) {
                ImGui::SameLine();
            }
            const bool active = editor.draw_mode->bond_order() == orders[i];
            if (toolbar_button(bond_order_name(orders[i]), active)) {
                editor.draw_mode->set_bond_order(orders[i]);
            }
        }
    } else if (editor.current_mode == EditorState::Mode::Fragment && editor.fragment_mode) {
        static int current_category = 0;
        const std::vector<std::string> categories = editor.fragment_library.categories();
        if (current_category >= static_cast<int>(categories.size())) {
            current_category = 0;
        }
        std::vector<const char*> labels;
        labels.reserve(categories.size());
        for (const auto& category : categories) {
            labels.push_back(category.c_str());
        }
        if (!labels.empty()) {
            ImGui::Combo("Category", &current_category, labels.data(), static_cast<int>(labels.size()));
            const std::vector<const sbox::editor::Fragment*> fragments =
                editor.fragment_library.by_category(categories[static_cast<std::size_t>(current_category)]);
            if (ImGui::BeginListBox("Fragments", ImVec2(-FLT_MIN, 180.0f))) {
                for (const sbox::editor::Fragment* fragment : fragments) {
                    const bool selected = editor.fragment_mode->selected_fragment() == fragment;
                    if (ImGui::Selectable(fragment->name.c_str(), selected)) {
                        editor.fragment_mode->set_fragment(fragment);
                    }
                }
                ImGui::EndListBox();
            }
        }
        if (const sbox::editor::Fragment* fragment = editor.fragment_mode->selected_fragment()) {
            ImGui::Text("Selected: %s", fragment->name.c_str());
            ImGui::Text("Atoms: %d", fragment->molecule.num_atoms());
        }
    } else if (editor.current_mode == EditorState::Mode::Measure && editor.measure_mode) {
        ImGui::TextUnformatted("Measurements");
        const auto& measurements = editor.measure_mode->measurements();
        if (measurements.empty()) {
            ImGui::TextDisabled("No measurements");
        } else {
            for (const auto& measurement : measurements) {
                ImGui::BulletText("%s", measurement.label.c_str());
            }
        }
        if (ImGui::Button("Clear All Measurements")) {
            editor.measure_mode->clear_measurements();
        }
    } else {
        ImGui::TextDisabled("%s mode", mode_name(editor.current_mode));
    }

    ImGui::Separator();
    if (ImGui::Button("Add Hydrogens")) {
        if (!editor.selection.atoms.empty()) {
            for (int atom_index : editor.selection.atoms) {
                editor.commands.execute(std::make_unique<sbox::editor::AddHydrogensCommand>(atom_index), mol);
            }
        } else {
            editor.commands.execute(std::make_unique<sbox::editor::AddHydrogensCommand>(), mol);
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Remove Hydrogens")) {
        editor.commands.execute(std::make_unique<sbox::editor::RemoveHydrogensCommand>(), mol);
    }
    if (ImGui::Button("Center Molecule")) {
        std::vector<int> indices = all_atom_indices(mol);
        std::vector<Eigen::Vector3d> positions;
        positions.reserve(indices.size());
        const Eigen::Vector3d center = mol.center_of_mass();
        for (int idx : indices) {
            positions.push_back(mol.atom(idx).position - center);
        }
        if (!indices.empty()) {
            editor.commands.execute(std::make_unique<sbox::editor::MoveAtomsCommand>(indices, positions), mol);
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Perceive Bonds")) {
        mol.perceive_bonds();
    }

    ImGui::Separator();
    const std::string undo_label = editor.commands.can_undo()
        ? "Undo: " + editor.commands.undo_description()
        : "Undo";
    const std::string redo_label = editor.commands.can_redo()
        ? "Redo: " + editor.commands.redo_description()
        : "Redo";

    if (!editor.commands.can_undo()) {
        ImGui::BeginDisabled();
    }
    if (ImGui::Button(undo_label.c_str())) {
        editor.commands.undo(mol);
    }
    if (!editor.commands.can_undo()) {
        ImGui::EndDisabled();
    }
    ImGui::SameLine();
    if (!editor.commands.can_redo()) {
        ImGui::BeginDisabled();
    }
    if (ImGui::Button(redo_label.c_str())) {
        editor.commands.redo(mol);
    }
    if (!editor.commands.can_redo()) {
        ImGui::EndDisabled();
    }

    ImGui::Separator();
    if (editor.selection.num_atoms() > 0) {
        ImGui::Text("%d atoms selected", editor.selection.num_atoms());
    }
    if (editor.selection.num_bonds() > 0) {
        ImGui::Text("%d bonds selected", editor.selection.num_bonds());
    }
    if (editor.selection.num_atoms() == 1) {
        const int atom_index = editor.selection.atoms.front();
        if (atom_index >= 0 && atom_index < mol.num_atoms()) {
            const auto& atom = mol.atom(atom_index);
            ImGui::Text("Atom: %s%d", sbox::elements::get_element(atom.Z).symbol, atom_index + 1);
            ImGui::Text("Position: (%.3f, %.3f, %.3f)", atom.position.x(), atom.position.y(), atom.position.z());
            ImGui::Text("Formal charge: %d", atom.formal_charge);
            ImGui::Text("Coordination: %d", mol.coordination_number(atom_index));
        }
    }
    if (editor.selection.num_bonds() == 1) {
        const int bond_index = editor.selection.bonds.front();
        if (bond_index >= 0 && bond_index < mol.num_bonds()) {
            const auto& bond = mol.bond(bond_index);
            ImGui::Text("Bond: %d-%d", bond.atom_i + 1, bond.atom_j + 1);
            ImGui::Text("Order: %s", bond_order_name(bond.order));
            ImGui::Text("Length: %.3f A", mol.distance(bond.atom_i, bond.atom_j) * kBohrToAngstrom);
        }
    }

    ImGui::End();
}

}  // namespace sbox::ui

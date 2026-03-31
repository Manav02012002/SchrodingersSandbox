#include "ui/annotation_editor.h"

#include "core/elements.h"

#include <imgui.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <string>

namespace sbox::ui {

namespace {

constexpr double kAngstromToBohr = 1.0 / 0.529177210903;

const char* type_icon(Annotation::Type type) {
    switch (type) {
    case Annotation::Type::Text: return "T";
    case Annotation::Type::Arrow: return "->";
    case Annotation::Type::DimensionLine: return "[]";
    case Annotation::Type::AngleArc: return "<)";
    case Annotation::Type::Circle: return "O";
    }
    return "?";
}

std::string atom_label(const sbox::chem::MolecularSystem& mol, int atom_index) {
    return std::string(sbox::elements::get_element(mol.atom(atom_index).Z).symbol) + std::to_string(atom_index + 1);
}

std::string preview_text(const Annotation& annotation) {
    if (!annotation.text.empty()) {
        return annotation.text;
    }
    switch (annotation.type) {
    case Annotation::Type::Text: return "Text";
    case Annotation::Type::Arrow: return "Arrow";
    case Annotation::Type::DimensionLine: return annotation.auto_value ? "Auto distance" : "Dimension";
    case Annotation::Type::AngleArc: return annotation.auto_value ? "Auto angle" : "Angle";
    case Annotation::Type::Circle: return "Circle";
    }
    return "Annotation";
}

void apply_publication_style(AnnotationManager& manager) {
    manager.set_default_style(ImVec4(1.0f, 1.0f, 1.0f, 1.0f), 1.0f, 1.5f);
}

void apply_presentation_style(AnnotationManager& manager) {
    manager.set_default_style(ImVec4(1.0f, 0.92f, 0.25f, 1.0f), 1.6f, 2.5f);
}

void apply_teaching_style(AnnotationManager& manager) {
    manager.set_default_style(ImVec4(0.35f, 0.92f, 1.0f, 1.0f), 1.3f, 3.0f);
}

}  // namespace

void draw_annotation_editor(AnnotationManager& manager,
                            AppState& state,
                            const sbox::chem::MolecularSystem& mol,
                            const sbox::editor::Selection& selection) {
    (void)state;
    if (!ImGui::Begin("Annotations")) {
        ImGui::End();
        return;
    }

    static int editing_id = -1;

    ImGui::TextUnformatted("Add from Selection");
    const auto& atoms = selection.atoms;
    if (atoms.size() == 1 && atoms[0] < mol.num_atoms()) {
        if (ImGui::Button("Add Label")) {
            const int atom_index = atoms[0];
            manager.add_text(mol.atom(atom_index).position, atom_label(mol, atom_index));
        }
    } else if (atoms.size() == 2 && atoms[0] < mol.num_atoms() && atoms[1] < mol.num_atoms()) {
        if (ImGui::Button("Add Dimension Line")) {
            manager.add_dimension_line(mol.atom(atoms[0]).position, mol.atom(atoms[1]).position, true);
        }
        ImGui::SameLine();
        if (ImGui::Button("Add Arrow")) {
            manager.add_arrow(mol.atom(atoms[0]).position, mol.atom(atoms[1]).position);
        }
    } else if (atoms.size() == 3 &&
               atoms[0] < mol.num_atoms() &&
               atoms[1] < mol.num_atoms() &&
               atoms[2] < mol.num_atoms()) {
        if (ImGui::Button("Add Angle Arc")) {
            manager.add_angle_arc(mol.atom(atoms[1]).position,
                                  mol.atom(atoms[0]).position,
                                  mol.atom(atoms[2]).position,
                                  true);
        }
    } else {
        ImGui::TextDisabled("Select 1, 2, or 3 atoms to create annotation geometry.");
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Quick Labels");
    if (ImGui::Button("Label All Atoms")) {
        for (int i = 0; i < mol.num_atoms(); ++i) {
            if (mol.atom(i).Z != 1) {
                manager.add_text(mol.atom(i).position, sbox::elements::get_element(mol.atom(i).Z).symbol);
            }
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Label Selected")) {
        for (int atom_index : atoms) {
            if (atom_index >= 0 && atom_index < mol.num_atoms()) {
                manager.add_text(mol.atom(atom_index).position, atom_label(mol, atom_index));
            }
        }
    }
    if (ImGui::Button("Add Scale Bar") && mol.num_atoms() > 0) {
        Eigen::Vector3d min_corner = mol.atom(0).position;
        Eigen::Vector3d max_corner = mol.atom(0).position;
        for (int i = 1; i < mol.num_atoms(); ++i) {
            min_corner = min_corner.cwiseMin(mol.atom(i).position);
            max_corner = max_corner.cwiseMax(mol.atom(i).position);
        }
        const Eigen::Vector3d margin(0.4, 0.6, 0.0);
        const Eigen::Vector3d start(min_corner.x(), min_corner.y() - margin.y(), min_corner.z());
        const Eigen::Vector3d end = start + Eigen::Vector3d(1.0 * kAngstromToBohr, 0.0, 0.0);
        const int id = manager.add_dimension_line(start, end, false);
        Annotation& scale_bar = manager.get(id);
        scale_bar.text = "1.00 \xC3\x85";
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Annotation List");
    int remove_id = -1;
    for (const Annotation& annotation_view : manager.all()) {
        Annotation& annotation = manager.get(annotation_view.id);
        ImGui::PushID(annotation.id);
        ImGui::Checkbox("##visible", &annotation.visible);
        ImGui::SameLine();
        ImGui::Text("%s %s", type_icon(annotation.type), preview_text(annotation).c_str());
        ImGui::SameLine();
        if (ImGui::Button("Edit")) {
            editing_id = (editing_id == annotation.id) ? -1 : annotation.id;
        }
        ImGui::SameLine();
        if (ImGui::Button("x")) {
            remove_id = annotation.id;
        }
        if (editing_id == annotation.id) {
            char buffer[256];
            std::snprintf(buffer, sizeof(buffer), "%s", annotation.text.c_str());
            if (ImGui::InputText("Text", buffer, sizeof(buffer))) {
                annotation.text = buffer;
                annotation.auto_value = false;
            }
            ImGui::ColorEdit4("Color", &annotation.color.x);
            ImGui::SliderFloat("Font Scale", &annotation.font_scale, 0.5f, 3.0f, "%.2f");
            ImGui::SliderFloat("Thickness", &annotation.thickness, 0.5f, 6.0f, "%.1f");
        }
        ImGui::Separator();
        ImGui::PopID();
    }
    if (remove_id >= 0) {
        if (editing_id == remove_id) {
            editing_id = -1;
        }
        manager.remove(remove_id);
    }
    if (ImGui::Button("Clear All Annotations")) {
        editing_id = -1;
        manager.clear();
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Style Presets");
    if (ImGui::Button("Publication Style")) {
        apply_publication_style(manager);
    }
    ImGui::SameLine();
    if (ImGui::Button("Presentation Style")) {
        apply_presentation_style(manager);
    }
    ImGui::SameLine();
    if (ImGui::Button("Teaching Style")) {
        apply_teaching_style(manager);
    }
    ImGui::Text("Defaults: color(%.2f, %.2f, %.2f), font %.2f, line %.1f",
                manager.default_color().x,
                manager.default_color().y,
                manager.default_color().z,
                manager.default_font_scale(),
                manager.default_thickness());

    ImGui::End();
}

}  // namespace sbox::ui

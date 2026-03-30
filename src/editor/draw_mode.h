#pragma once

#include "editor/editor_mode.h"

namespace sbox::editor {

class DrawMode : public EditorMode {
public:
    void on_mouse_down(const Ray& ray, int button, bool shift,
                       sbox::chem::MolecularSystem& mol,
                       Selection& selection,
                       CommandStack& commands) override;

    void on_mouse_up(const Ray& ray, int button,
                     sbox::chem::MolecularSystem& mol,
                     Selection& selection,
                     CommandStack& commands) override;

    void on_mouse_move(const Ray& ray, float dx, float dy, bool dragging,
                       sbox::chem::MolecularSystem& mol,
                       Selection& selection,
                       CommandStack& commands) override;

    void on_key(int key, bool ctrl, bool shift,
                sbox::chem::MolecularSystem& mol,
                Selection& selection,
                CommandStack& commands) override;

    void draw_overlay(ImDrawList* draw_list,
                      const sbox::chem::MolecularSystem& mol,
                      const Selection& selection,
                      const Eigen::Matrix4f& vp_matrix,
                      const ImVec2& viewport_pos,
                      const ImVec2& viewport_size) override;

    const char* name() const override { return "Draw"; }
    const char* cursor() const override { return "crosshair"; }

    void set_element(int Z) { current_element_ = Z; }
    int element() const { return current_element_; }
    void set_bond_order(sbox::chem::BondOrder order) { current_bond_order_ = order; }
    sbox::chem::BondOrder bond_order() const { return current_bond_order_; }

private:
    int current_element_ = 6;
    sbox::chem::BondOrder current_bond_order_ = sbox::chem::BondOrder::Single;

    bool drawing_bond_ = false;
    int bond_start_atom_ = -1;
    int hovered_atom_ = -1;
    Eigen::Vector3f bond_preview_end_ = Eigen::Vector3f::Zero();
    Eigen::Vector3f draw_plane_normal_ = Eigen::Vector3f::UnitZ();
    float draw_plane_d_ = 0.0f;
};

}  // namespace sbox::editor

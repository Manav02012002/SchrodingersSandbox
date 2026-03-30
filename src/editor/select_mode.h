#pragma once

#include "editor/editor_mode.h"

namespace sbox::ui {
struct ContextMenuState;
}

namespace sbox::editor {

class SelectMode : public EditorMode {
public:
    void set_context_menu_state(sbox::ui::ContextMenuState* context_menu_state) { context_menu_state_ = context_menu_state; }

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

    const char* name() const override { return "Select"; }

private:
    bool dragging_ = false;
    bool is_drag_move_ = false;
    Eigen::Vector3f drag_start_ = Eigen::Vector3f::Zero();
    Ray drag_start_ray_{Eigen::Vector3f::Zero(), Eigen::Vector3f::UnitZ()};
    std::vector<int> drag_atom_indices_;
    std::vector<Eigen::Vector3d> drag_original_positions_;
    Eigen::Vector3f drag_plane_normal_ = Eigen::Vector3f::UnitZ();
    float drag_plane_d_ = 0.0f;
    sbox::ui::ContextMenuState* context_menu_state_ = nullptr;
};

}  // namespace sbox::editor

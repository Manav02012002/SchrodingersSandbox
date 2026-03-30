#pragma once

#include "editor/editor_mode.h"

namespace sbox::editor {

class EraseMode : public EditorMode {
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

    const char* name() const override { return "Erase"; }
    const char* cursor() const override { return "eraser"; }

private:
    int hover_atom_ = -1;
    int hover_bond_ = -1;
};

}  // namespace sbox::editor

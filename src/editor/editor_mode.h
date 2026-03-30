#pragma once

#include "editor/command.h"
#include "editor/picking.h"

#include <Eigen/Core>

#include <imgui.h>

namespace sbox::editor {

class EditorMode {
public:
    virtual ~EditorMode() = default;

    virtual void on_mouse_down(const Ray& ray, int button, bool shift,
                               sbox::chem::MolecularSystem& mol,
                               Selection& selection,
                               CommandStack& commands) = 0;

    virtual void on_mouse_up(const Ray& ray, int button,
                             sbox::chem::MolecularSystem& mol,
                             Selection& selection,
                             CommandStack& commands) = 0;

    virtual void on_mouse_move(const Ray& ray, float dx, float dy, bool dragging,
                               sbox::chem::MolecularSystem& mol,
                               Selection& selection,
                               CommandStack& commands) = 0;

    virtual void on_key(int key, bool ctrl, bool shift,
                        sbox::chem::MolecularSystem& mol,
                        Selection& selection,
                        CommandStack& commands) = 0;

    virtual void draw_overlay(ImDrawList* draw_list,
                              const sbox::chem::MolecularSystem& mol,
                              const Selection& selection,
                              const Eigen::Matrix4f& vp_matrix,
                              const ImVec2& viewport_pos,
                              const ImVec2& viewport_size) = 0;

    virtual const char* name() const = 0;
    virtual const char* cursor() const { return "arrow"; }
};

}  // namespace sbox::editor

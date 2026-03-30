#pragma once

#include "editor/editor_mode.h"
#include "editor/fragment_library.h"

namespace sbox::editor {

class FragmentMode : public EditorMode {
public:
    explicit FragmentMode(const FragmentLibrary* library);

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

    const char* name() const override { return "Fragment"; }

    void set_fragment(const Fragment* frag);
    const Fragment* selected_fragment() const;

private:
    const FragmentLibrary* library_ = nullptr;
    const Fragment* selected_fragment_ = nullptr;
    Eigen::Vector3f preview_position_ = Eigen::Vector3f::Zero();
    bool show_preview_ = false;
    int hover_atom_ = -1;
};

}  // namespace sbox::editor

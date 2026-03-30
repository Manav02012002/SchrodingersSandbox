#pragma once

#include "editor/command.h"
#include "editor/draw_mode.h"
#include "editor/erase_mode.h"
#include "editor/fragment_library.h"
#include "editor/fragment_mode.h"
#include "editor/measure_mode.h"
#include "editor/select_mode.h"
#include "ui/context_menu.h"

#include <memory>

namespace sbox::ui {

struct EditorState {
    enum class Mode { Select, Draw, Erase, Measure, Fragment };

    EditorState();

    Mode current_mode = Mode::Select;

    std::unique_ptr<sbox::editor::SelectMode> select_mode;
    std::unique_ptr<sbox::editor::DrawMode> draw_mode;
    std::unique_ptr<sbox::editor::EraseMode> erase_mode;
    std::unique_ptr<sbox::editor::MeasureMode> measure_mode;
    std::unique_ptr<sbox::editor::FragmentMode> fragment_mode;

    sbox::editor::Selection selection;
    sbox::editor::CommandStack commands;
    sbox::editor::FragmentLibrary fragment_library;
    ContextMenuState context_menu;

    sbox::editor::EditorMode* active_mode();
    const sbox::editor::EditorMode* active_mode() const;
};

void draw_editor_toolbar(EditorState& editor, sbox::chem::MolecularSystem& mol);

}  // namespace sbox::ui

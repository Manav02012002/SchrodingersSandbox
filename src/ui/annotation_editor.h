#pragma once

#include "core/molecular_system.h"
#include "editor/picking.h"
#include "ui/annotations.h"
#include "ui/app_state.h"

namespace sbox::ui {

void draw_annotation_editor(AnnotationManager& manager,
                            AppState& state,
                            const sbox::chem::MolecularSystem& mol,
                            const sbox::editor::Selection& selection);

}  // namespace sbox::ui

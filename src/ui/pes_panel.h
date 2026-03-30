#pragma once

#include "backend/backend_manager.h"
#include "core/molecular_system.h"
#include "editor/picking.h"
#include "ui/app_state.h"

namespace sbox::ui {

void draw_pes_panel(AppState& state,
                    sbox::backend::BackendManager& backend,
                    const sbox::chem::MolecularSystem& mol,
                    const sbox::editor::Selection& selection);

}  // namespace sbox::ui

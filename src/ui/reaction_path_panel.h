#pragma once

#include "backend/backend_manager.h"
#include "core/molecular_system.h"
#include "renderer/mol_renderer.h"
#include "ui/app_state.h"

namespace sbox::ui {

void draw_reaction_path_panel(AppState& state,
                              sbox::backend::BackendManager& backend,
                              sbox::chem::MolecularSystem& mol,
                              sbox::render::MolRenderer& renderer);

}  // namespace sbox::ui

#pragma once

#include "backend/job_types.h"
#include "core/molecular_system.h"
#include "renderer/mol_renderer.h"
#include "ui/app_state.h"

namespace sbox::ui {

void draw_optimization_panel(AppState& state,
                             const sbox::backend::JobResult& result,
                             sbox::chem::MolecularSystem& mol,
                             sbox::render::MolRenderer& renderer);

}  // namespace sbox::ui

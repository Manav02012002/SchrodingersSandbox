#pragma once

#include "backend/job_types.h"
#include "core/molecular_system.h"
#include "ui/app_state.h"

namespace sbox::ui {

void draw_bond_order_panel(AppState& state,
                           const sbox::backend::JobResult& result,
                           const sbox::chem::MolecularSystem& mol);

}  // namespace sbox::ui

#pragma once

#include "analysis/crystal_field.h"
#include "backend/job_types.h"
#include "core/molecular_system.h"
#include "ui/app_state.h"

namespace sbox::ui {

void draw_d_orbital_viewer(AppState& state,
                           const sbox::backend::JobResult& result,
                           const sbox::chem::MolecularSystem& mol,
                           const sbox::analysis::DOrbitalEnergies& d_orbs);

}  // namespace sbox::ui

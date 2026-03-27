#pragma once

#include "core/molecular_system.h"
#include "ui/app_state.h"

namespace sbox::ui {

void draw_molecule_info(const AppState& state, const sbox::chem::MolecularSystem& mol);

}  // namespace sbox::ui

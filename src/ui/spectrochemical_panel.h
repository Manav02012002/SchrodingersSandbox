#pragma once

#include "backend/backend_manager.h"
#include "chem/ligand_library.h"
#include "ui/app_state.h"

namespace sbox::ui {

void draw_spectrochemical_panel(AppState& state,
                                sbox::backend::BackendManager& backend,
                                const sbox::chem::LigandLibrary& ligand_library);

}  // namespace sbox::ui

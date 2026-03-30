#pragma once

#include "chem/ligand_library.h"
#include "editor/command.h"
#include "ui/app_state.h"

namespace sbox::ui {

void draw_complex_builder(AppState& state,
                          sbox::chem::MolecularSystem& mol,
                          sbox::editor::CommandStack& commands,
                          const sbox::chem::LigandLibrary& ligand_library);

}  // namespace sbox::ui

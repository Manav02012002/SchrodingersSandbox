#pragma once

#include "core/molecular_system.h"

#include <string>

namespace sbox::chem {

enum class PointGroup {
    C1, Ci, Cs,
    C2, C3, C4, C5, C6,
    C2v, C3v, C4v, C5v, C6v,
    C2h, C3h, C4h, C5h, C6h,
    D2, D3, D4, D5, D6,
    D2h, D3h, D4h, D5h, D6h,
    D2d, D3d, D4d, D5d, D6d,
    S4, S6, S8,
    T, Td, Th,
    O, Oh,
    Ih,
    Cinfv, Dinfh,
    Unknown
};

std::string point_group_name(PointGroup pg);
PointGroup detect_point_group(const MolecularSystem& mol, double tolerance = 0.1);

}  // namespace sbox::chem

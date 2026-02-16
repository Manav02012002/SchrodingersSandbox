#pragma once

#include "core/slater.h"

#include <array>
#include <string>

namespace sbox::elements {

struct Element {
    int Z;
    const char* symbol;
    const char* name;
    sbox::slater::ElectronConfig config;
    const char* category;
    double atomic_mass;
    double covalent_radius;
    double vdw_radius;
};

struct PTPosition {
    int Z;
    int row;
    int col;
};

extern const std::array<Element, 118> ALL_ELEMENTS;
extern const std::array<PTPosition, 118> PT_LAYOUT;

const Element& get_element(int Z);
const Element& get_element(const std::string& symbol);

}  // namespace sbox::elements

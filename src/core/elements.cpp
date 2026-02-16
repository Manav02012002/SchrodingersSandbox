#include "core/elements.h"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace sbox::elements {

namespace {

using Orbital = std::tuple<int, int, int>;

sbox::slater::ElectronConfig build_aufbau_config(int z) {
    static const std::array<Orbital, 19> order = {
        Orbital{1, 0, 2},  Orbital{2, 0, 2},  Orbital{2, 1, 6},  Orbital{3, 0, 2},  Orbital{3, 1, 6},
        Orbital{4, 0, 2},  Orbital{3, 2, 10}, Orbital{4, 1, 6},  Orbital{5, 0, 2},  Orbital{4, 2, 10},
        Orbital{5, 1, 6},  Orbital{6, 0, 2},  Orbital{4, 3, 14}, Orbital{5, 2, 10}, Orbital{6, 1, 6},
        Orbital{7, 0, 2},  Orbital{5, 3, 14}, Orbital{6, 2, 10}, Orbital{7, 1, 6},
    };

    sbox::slater::ElectronConfig config;
    int remaining = z;
    for (const auto& [n, l, cap] : order) {
        if (remaining <= 0) {
            break;
        }
        const int electrons = std::min(remaining, cap);
        if (electrons > 0) {
            config.push_back({n, l, electrons});
            remaining -= electrons;
        }
    }
    return config;
}

void set_subshell(sbox::slater::ElectronConfig& config, int n, int l, int electrons) {
    auto it = std::find_if(config.begin(), config.end(),
                           [n, l](const sbox::slater::SubshellConfig& s) {
                               return s.n == n && s.l == l;
                           });
    if (it != config.end()) {
        if (electrons <= 0) {
            config.erase(it);
        } else {
            it->electrons = electrons;
        }
        return;
    }

    if (electrons > 0) {
        config.push_back({n, l, electrons});
    }
}

void apply_exceptions(int z, sbox::slater::ElectronConfig& config) {
    switch (z) {
        case 24:
            set_subshell(config, 4, 0, 1);
            set_subshell(config, 3, 2, 5);
            break;
        case 29:
            set_subshell(config, 4, 0, 1);
            set_subshell(config, 3, 2, 10);
            break;
        case 41:
            set_subshell(config, 5, 0, 1);
            set_subshell(config, 4, 2, 4);
            break;
        case 42:
            set_subshell(config, 5, 0, 1);
            set_subshell(config, 4, 2, 5);
            break;
        case 44:
            set_subshell(config, 5, 0, 1);
            set_subshell(config, 4, 2, 7);
            break;
        case 45:
            set_subshell(config, 5, 0, 1);
            set_subshell(config, 4, 2, 8);
            break;
        case 46:
            set_subshell(config, 5, 0, 0);
            set_subshell(config, 4, 2, 10);
            break;
        case 47:
            set_subshell(config, 5, 0, 1);
            set_subshell(config, 4, 2, 10);
            break;
        case 57:
            set_subshell(config, 4, 3, 0);
            set_subshell(config, 5, 2, 1);
            break;
        case 58:
            set_subshell(config, 4, 3, 1);
            set_subshell(config, 5, 2, 1);
            break;
        case 64:
            set_subshell(config, 4, 3, 7);
            set_subshell(config, 5, 2, 1);
            break;
        case 78:
            set_subshell(config, 6, 0, 1);
            set_subshell(config, 5, 2, 9);
            break;
        case 79:
            set_subshell(config, 6, 0, 1);
            set_subshell(config, 5, 2, 10);
            break;
        case 89:
            set_subshell(config, 5, 3, 0);
            set_subshell(config, 6, 2, 1);
            break;
        case 90:
            set_subshell(config, 5, 3, 0);
            set_subshell(config, 6, 2, 2);
            break;
        case 91:
            set_subshell(config, 5, 3, 2);
            set_subshell(config, 6, 2, 1);
            break;
        case 92:
            set_subshell(config, 5, 3, 3);
            set_subshell(config, 6, 2, 1);
            break;
        case 93:
            set_subshell(config, 5, 3, 4);
            set_subshell(config, 6, 2, 1);
            break;
        case 96:
            set_subshell(config, 5, 3, 7);
            set_subshell(config, 6, 2, 1);
            break;
        case 103:
            set_subshell(config, 6, 2, 0);
            set_subshell(config, 7, 1, 1);
            break;
        default:
            break;
    }
}

const char* symbol_for_z(int z) {
    switch (z) {
        case 1: return "H"; case 2: return "He"; case 3: return "Li"; case 4: return "Be"; case 5: return "B";
        case 6: return "C"; case 7: return "N"; case 8: return "O"; case 9: return "F"; case 10: return "Ne";
        case 11: return "Na"; case 12: return "Mg"; case 13: return "Al"; case 14: return "Si"; case 15: return "P";
        case 16: return "S"; case 17: return "Cl"; case 18: return "Ar"; case 19: return "K"; case 20: return "Ca";
        case 21: return "Sc"; case 22: return "Ti"; case 23: return "V"; case 24: return "Cr"; case 25: return "Mn";
        case 26: return "Fe"; case 27: return "Co"; case 28: return "Ni"; case 29: return "Cu"; case 30: return "Zn";
        case 31: return "Ga"; case 32: return "Ge"; case 33: return "As"; case 34: return "Se"; case 35: return "Br";
        case 36: return "Kr"; case 37: return "Rb"; case 38: return "Sr"; case 39: return "Y"; case 40: return "Zr";
        case 41: return "Nb"; case 42: return "Mo"; case 43: return "Tc"; case 44: return "Ru"; case 45: return "Rh";
        case 46: return "Pd"; case 47: return "Ag"; case 48: return "Cd"; case 49: return "In"; case 50: return "Sn";
        case 51: return "Sb"; case 52: return "Te"; case 53: return "I"; case 54: return "Xe"; case 55: return "Cs";
        case 56: return "Ba"; case 57: return "La"; case 58: return "Ce"; case 59: return "Pr"; case 60: return "Nd";
        case 61: return "Pm"; case 62: return "Sm"; case 63: return "Eu"; case 64: return "Gd"; case 65: return "Tb";
        case 66: return "Dy"; case 67: return "Ho"; case 68: return "Er"; case 69: return "Tm"; case 70: return "Yb";
        case 71: return "Lu"; case 72: return "Hf"; case 73: return "Ta"; case 74: return "W"; case 75: return "Re";
        case 76: return "Os"; case 77: return "Ir"; case 78: return "Pt"; case 79: return "Au"; case 80: return "Hg";
        case 81: return "Tl"; case 82: return "Pb"; case 83: return "Bi"; case 84: return "Po"; case 85: return "At";
        case 86: return "Rn"; case 87: return "Fr"; case 88: return "Ra"; case 89: return "Ac"; case 90: return "Th";
        case 91: return "Pa"; case 92: return "U"; case 93: return "Np"; case 94: return "Pu"; case 95: return "Am";
        case 96: return "Cm"; case 97: return "Bk"; case 98: return "Cf"; case 99: return "Es"; case 100: return "Fm";
        case 101: return "Md"; case 102: return "No"; case 103: return "Lr"; case 104: return "Rf"; case 105: return "Db";
        case 106: return "Sg"; case 107: return "Bh"; case 108: return "Hs"; case 109: return "Mt"; case 110: return "Ds";
        case 111: return "Rg"; case 112: return "Cn"; case 113: return "Nh"; case 114: return "Fl"; case 115: return "Mc";
        case 116: return "Lv"; case 117: return "Ts"; case 118: return "Og";
        default: return "";
    }
}

const char* name_for_z(int z) {
    switch (z) {
        case 1: return "Hydrogen"; case 2: return "Helium"; case 3: return "Lithium"; case 4: return "Beryllium"; case 5: return "Boron";
        case 6: return "Carbon"; case 7: return "Nitrogen"; case 8: return "Oxygen"; case 9: return "Fluorine"; case 10: return "Neon";
        case 11: return "Sodium"; case 12: return "Magnesium"; case 13: return "Aluminium"; case 14: return "Silicon"; case 15: return "Phosphorus";
        case 16: return "Sulfur"; case 17: return "Chlorine"; case 18: return "Argon"; case 19: return "Potassium"; case 20: return "Calcium";
        case 21: return "Scandium"; case 22: return "Titanium"; case 23: return "Vanadium"; case 24: return "Chromium"; case 25: return "Manganese";
        case 26: return "Iron"; case 27: return "Cobalt"; case 28: return "Nickel"; case 29: return "Copper"; case 30: return "Zinc";
        case 31: return "Gallium"; case 32: return "Germanium"; case 33: return "Arsenic"; case 34: return "Selenium"; case 35: return "Bromine";
        case 36: return "Krypton"; case 37: return "Rubidium"; case 38: return "Strontium"; case 39: return "Yttrium"; case 40: return "Zirconium";
        case 41: return "Niobium"; case 42: return "Molybdenum"; case 43: return "Technetium"; case 44: return "Ruthenium"; case 45: return "Rhodium";
        case 46: return "Palladium"; case 47: return "Silver"; case 48: return "Cadmium"; case 49: return "Indium"; case 50: return "Tin";
        case 51: return "Antimony"; case 52: return "Tellurium"; case 53: return "Iodine"; case 54: return "Xenon"; case 55: return "Cesium";
        case 56: return "Barium"; case 57: return "Lanthanum"; case 58: return "Cerium"; case 59: return "Praseodymium"; case 60: return "Neodymium";
        case 61: return "Promethium"; case 62: return "Samarium"; case 63: return "Europium"; case 64: return "Gadolinium"; case 65: return "Terbium";
        case 66: return "Dysprosium"; case 67: return "Holmium"; case 68: return "Erbium"; case 69: return "Thulium"; case 70: return "Ytterbium";
        case 71: return "Lutetium"; case 72: return "Hafnium"; case 73: return "Tantalum"; case 74: return "Tungsten"; case 75: return "Rhenium";
        case 76: return "Osmium"; case 77: return "Iridium"; case 78: return "Platinum"; case 79: return "Gold"; case 80: return "Mercury";
        case 81: return "Thallium"; case 82: return "Lead"; case 83: return "Bismuth"; case 84: return "Polonium"; case 85: return "Astatine";
        case 86: return "Radon"; case 87: return "Francium"; case 88: return "Radium"; case 89: return "Actinium"; case 90: return "Thorium";
        case 91: return "Protactinium"; case 92: return "Uranium"; case 93: return "Neptunium"; case 94: return "Plutonium"; case 95: return "Americium";
        case 96: return "Curium"; case 97: return "Berkelium"; case 98: return "Californium"; case 99: return "Einsteinium"; case 100: return "Fermium";
        case 101: return "Mendelevium"; case 102: return "Nobelium"; case 103: return "Lawrencium"; case 104: return "Rutherfordium"; case 105: return "Dubnium";
        case 106: return "Seaborgium"; case 107: return "Bohrium"; case 108: return "Hassium"; case 109: return "Meitnerium"; case 110: return "Darmstadtium";
        case 111: return "Roentgenium"; case 112: return "Copernicium"; case 113: return "Nihonium"; case 114: return "Flerovium"; case 115: return "Moscovium";
        case 116: return "Livermorium"; case 117: return "Tennessine"; case 118: return "Oganesson";
        default: return "";
    }
}

const char* category_for_z(int z) {
    switch (z) {
        case 1:
        case 6:
        case 7:
        case 8:
        case 15:
        case 16:
        case 34:
            return "nonmetal";
        case 3:
        case 11:
        case 19:
        case 37:
        case 55:
        case 87:
            return "alkali";
        case 4:
        case 12:
        case 20:
        case 38:
        case 56:
        case 88:
            return "alkaline";
        case 5:
        case 14:
        case 32:
        case 33:
        case 51:
        case 52:
            return "metalloid";
        case 9:
        case 17:
        case 35:
        case 53:
        case 85:
        case 117:
            return "halogen";
        case 2:
        case 10:
        case 18:
        case 36:
        case 54:
        case 86:
        case 118:
            return "noble";
        case 13:
        case 31:
        case 49:
        case 50:
        case 81:
        case 82:
        case 83:
        case 84:
        case 113:
        case 114:
        case 115:
        case 116:
            return "post-transition";
        default:
            if (z >= 57 && z <= 71) {
                return "lanthanide";
            }
            if (z >= 89 && z <= 103) {
                return "actinide";
            }
            return "transition";
    }
}

double mass_for_z(int z) {
    switch (z) {
        case 1: return 1.008; case 2: return 4.002602; case 3: return 6.94; case 4: return 9.0121831; case 5: return 10.81;
        case 6: return 12.011; case 7: return 14.007; case 8: return 15.999; case 9: return 18.998403163; case 10: return 20.1797;
        case 11: return 22.98976928; case 12: return 24.305; case 13: return 26.9815385; case 14: return 28.085; case 15: return 30.973761998;
        case 16: return 32.06; case 17: return 35.45; case 18: return 39.948; case 19: return 39.0983; case 20: return 40.078;
        case 21: return 44.955908; case 22: return 47.867; case 23: return 50.9415; case 24: return 51.9961; case 25: return 54.938044;
        case 26: return 55.845; case 27: return 58.933194; case 28: return 58.6934; case 29: return 63.546; case 30: return 65.38;
        case 31: return 69.723; case 32: return 72.630; case 33: return 74.921595; case 34: return 78.971; case 35: return 79.904;
        case 36: return 83.798; case 37: return 85.4678; case 38: return 87.62; case 39: return 88.90584; case 40: return 91.224;
        case 41: return 92.90637; case 42: return 95.95; case 43: return 98.0; case 44: return 101.07; case 45: return 102.90550;
        case 46: return 106.42; case 47: return 107.8682; case 48: return 112.414; case 49: return 114.818; case 50: return 118.710;
        case 51: return 121.760; case 52: return 127.60; case 53: return 126.90447; case 54: return 131.293; case 55: return 132.90545196;
        case 56: return 137.327; case 57: return 138.90547; case 58: return 140.116; case 59: return 140.90766; case 60: return 144.242;
        case 61: return 145.0; case 62: return 150.36; case 63: return 151.964; case 64: return 157.25; case 65: return 158.92535;
        case 66: return 162.500; case 67: return 164.93033; case 68: return 167.259; case 69: return 168.93422; case 70: return 173.045;
        case 71: return 174.9668; case 72: return 178.49; case 73: return 180.94788; case 74: return 183.84; case 75: return 186.207;
        case 76: return 190.23; case 77: return 192.217; case 78: return 195.084; case 79: return 196.966569; case 80: return 200.592;
        case 81: return 204.38; case 82: return 207.2; case 83: return 208.98040; case 84: return 209.0; case 85: return 210.0;
        case 86: return 222.0; case 87: return 223.0; case 88: return 226.0; case 89: return 227.0; case 90: return 232.0377;
        case 91: return 231.03588; case 92: return 238.02891; case 93: return 237.0; case 94: return 244.0; case 95: return 243.0;
        case 96: return 247.0; case 97: return 247.0; case 98: return 251.0; case 99: return 252.0; case 100: return 257.0;
        case 101: return 258.0; case 102: return 259.0; case 103: return 266.0; case 104: return 267.0; case 105: return 268.0;
        case 106: return 269.0; case 107: return 270.0; case 108: return 269.0; case 109: return 278.0; case 110: return 281.0;
        case 111: return 282.0; case 112: return 285.0; case 113: return 286.0; case 114: return 289.0; case 115: return 290.0;
        case 116: return 293.0; case 117: return 294.0; case 118: return 294.0;
        default: return 0.0;
    }
}

double covalent_radius_for_z(int z) {
    switch (z) {
        case 1: return 0.31; case 2: return 0.28; case 3: return 1.28; case 4: return 0.96; case 5: return 0.84;
        case 6: return 0.76; case 7: return 0.71; case 8: return 0.66; case 9: return 0.57; case 10: return 0.58;
        case 11: return 1.66; case 12: return 1.41; case 13: return 1.21; case 14: return 1.11; case 15: return 1.07;
        case 16: return 1.05; case 17: return 1.02; case 18: return 1.06; case 19: return 2.03; case 20: return 1.76;
        case 21: return 1.70; case 22: return 1.60; case 23: return 1.53; case 24: return 1.39; case 25: return 1.39;
        case 26: return 1.32; case 27: return 1.26; case 28: return 1.24; case 29: return 1.32; case 30: return 1.22;
        case 31: return 1.22; case 32: return 1.20; case 33: return 1.19; case 34: return 1.20; case 35: return 1.20;
        case 36: return 1.16; case 37: return 2.20; case 38: return 1.95; case 39: return 1.90; case 40: return 1.75;
        case 41: return 1.64; case 42: return 1.54; case 43: return 1.47; case 44: return 1.46; case 45: return 1.42;
        case 46: return 1.39; case 47: return 1.45; case 48: return 1.44; case 49: return 1.42; case 50: return 1.39;
        case 51: return 1.39; case 52: return 1.38; case 53: return 1.39; case 54: return 1.40; case 55: return 2.44;
        case 56: return 2.15; case 57: return 2.07; case 58: return 2.04; case 59: return 2.03; case 60: return 2.01;
        case 61: return 1.99; case 62: return 1.98; case 63: return 1.98; case 64: return 1.96; case 65: return 1.94;
        case 66: return 1.92; case 67: return 1.92; case 68: return 1.89; case 69: return 1.90; case 70: return 1.87;
        case 71: return 1.87; case 72: return 1.75; case 73: return 1.70; case 74: return 1.62; case 75: return 1.51;
        case 76: return 1.44; case 77: return 1.41; case 78: return 1.36; case 79: return 1.36; case 80: return 1.32;
        case 81: return 1.45; case 82: return 1.46; case 83: return 1.48; case 84: return 1.40; case 85: return 1.50;
        case 86: return 1.50; case 87: return 2.60; case 88: return 2.21; case 89: return 2.15; case 90: return 2.06;
        case 91: return 2.00; case 92: return 1.96; case 93: return 1.90; case 94: return 1.87; case 95: return 1.80;
        case 96: return 1.69; case 97: return 1.68; case 98: return 1.68; case 99: return 1.65; case 100: return 1.67;
        case 101: return 1.73; case 102: return 1.76; case 103: return 1.61; case 104: return 1.57; case 105: return 1.49;
        case 106: return 1.43; case 107: return 1.41; case 108: return 1.34; case 109: return 1.29; case 110: return 1.28;
        case 111: return 1.21; case 112: return 1.22; case 113: return 1.36; case 114: return 1.43; case 115: return 1.62;
        case 116: return 1.75; case 117: return 1.65; case 118: return 1.57;
        default: return 1.50;
    }
}

double vdw_radius_for_z(int z) {
    switch (z) {
        case 1: return 1.20; case 2: return 1.40; case 3: return 1.82; case 4: return 1.53; case 5: return 1.92;
        case 6: return 1.70; case 7: return 1.55; case 8: return 1.52; case 9: return 1.47; case 10: return 1.54;
        case 11: return 2.27; case 12: return 1.73; case 13: return 1.84; case 14: return 2.10; case 15: return 1.80;
        case 16: return 1.80; case 17: return 1.75; case 18: return 1.88; case 19: return 2.75; case 20: return 2.31;
        case 21: return 2.11; case 22: return 2.00; case 23: return 2.00; case 24: return 2.00; case 25: return 2.00;
        case 26: return 2.00; case 27: return 2.00; case 28: return 1.63; case 29: return 1.40; case 30: return 1.39;
        case 31: return 1.87; case 32: return 2.11; case 33: return 1.85; case 34: return 1.90; case 35: return 1.85;
        case 36: return 2.02; case 37: return 3.03; case 38: return 2.49; case 39: return 2.00; case 40: return 2.00;
        case 41: return 2.00; case 42: return 2.00; case 43: return 2.00; case 44: return 2.00; case 45: return 2.00;
        case 46: return 1.63; case 47: return 1.72; case 48: return 1.58; case 49: return 1.93; case 50: return 2.17;
        case 51: return 2.06; case 52: return 2.06; case 53: return 1.98; case 54: return 2.16;
        default: return 2.00;
    }
}

std::array<Element, 118> build_elements() {
    std::array<Element, 118> elements{};
    for (int z = 1; z <= 118; ++z) {
        auto config = build_aufbau_config(z);
        apply_exceptions(z, config);
        elements[static_cast<std::size_t>(z - 1)] = {
            z,
            symbol_for_z(z),
            name_for_z(z),
            std::move(config),
            category_for_z(z),
            mass_for_z(z),
            covalent_radius_for_z(z),
            vdw_radius_for_z(z),
        };
    }
    return elements;
}

std::array<PTPosition, 118> build_layout() {
    std::array<PTPosition, 118> layout{};

    auto set = [&layout](int z, int row, int col) {
        layout[static_cast<std::size_t>(z - 1)] = {z, row, col};
    };

    set(1, 0, 0); set(2, 0, 17);

    set(3, 1, 0); set(4, 1, 1); set(5, 1, 12); set(6, 1, 13); set(7, 1, 14); set(8, 1, 15); set(9, 1, 16); set(10, 1, 17);
    set(11, 2, 0); set(12, 2, 1); set(13, 2, 12); set(14, 2, 13); set(15, 2, 14); set(16, 2, 15); set(17, 2, 16); set(18, 2, 17);

    set(19, 3, 0); set(20, 3, 1); set(21, 3, 2); set(22, 3, 3); set(23, 3, 4); set(24, 3, 5); set(25, 3, 6); set(26, 3, 7); set(27, 3, 8);
    set(28, 3, 9); set(29, 3, 10); set(30, 3, 11); set(31, 3, 12); set(32, 3, 13); set(33, 3, 14); set(34, 3, 15); set(35, 3, 16); set(36, 3, 17);

    set(37, 4, 0); set(38, 4, 1); set(39, 4, 2); set(40, 4, 3); set(41, 4, 4); set(42, 4, 5); set(43, 4, 6); set(44, 4, 7); set(45, 4, 8);
    set(46, 4, 9); set(47, 4, 10); set(48, 4, 11); set(49, 4, 12); set(50, 4, 13); set(51, 4, 14); set(52, 4, 15); set(53, 4, 16); set(54, 4, 17);

    set(55, 5, 0); set(56, 5, 1); set(72, 5, 3); set(73, 5, 4); set(74, 5, 5); set(75, 5, 6); set(76, 5, 7); set(77, 5, 8); set(78, 5, 9);
    set(79, 5, 10); set(80, 5, 11); set(81, 5, 12); set(82, 5, 13); set(83, 5, 14); set(84, 5, 15); set(85, 5, 16); set(86, 5, 17);

    set(87, 6, 0); set(88, 6, 1); set(104, 6, 3); set(105, 6, 4); set(106, 6, 5); set(107, 6, 6); set(108, 6, 7); set(109, 6, 8); set(110, 6, 9);
    set(111, 6, 10); set(112, 6, 11); set(113, 6, 12); set(114, 6, 13); set(115, 6, 14); set(116, 6, 15); set(117, 6, 16); set(118, 6, 17);

    set(57, 8, 2); set(58, 8, 3); set(59, 8, 4); set(60, 8, 5); set(61, 8, 6); set(62, 8, 7); set(63, 8, 8); set(64, 8, 9); set(65, 8, 10);
    set(66, 8, 11); set(67, 8, 12); set(68, 8, 13); set(69, 8, 14); set(70, 8, 15); set(71, 8, 16);

    set(89, 9, 2); set(90, 9, 3); set(91, 9, 4); set(92, 9, 5); set(93, 9, 6); set(94, 9, 7); set(95, 9, 8); set(96, 9, 9); set(97, 9, 10);
    set(98, 9, 11); set(99, 9, 12); set(100, 9, 13); set(101, 9, 14); set(102, 9, 15); set(103, 9, 16);

    return layout;
}

}  // namespace

const std::array<Element, 118> ALL_ELEMENTS = build_elements();
const std::array<PTPosition, 118> PT_LAYOUT = build_layout();

const Element& get_element(int z) {
    if (z < 1 || z > 118) {
        throw std::out_of_range("Atomic number out of range");
    }
    return ALL_ELEMENTS[static_cast<std::size_t>(z - 1)];
}

const Element& get_element(const std::string& symbol) {
    const auto it = std::find_if(ALL_ELEMENTS.begin(), ALL_ELEMENTS.end(),
                                 [&symbol](const Element& e) {
                                     return symbol == e.symbol;
                                 });
    if (it == ALL_ELEMENTS.end()) {
        throw std::out_of_range("Element symbol not found");
    }
    return *it;
}

}  // namespace sbox::elements

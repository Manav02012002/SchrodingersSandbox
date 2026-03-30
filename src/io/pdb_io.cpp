#include "io/pdb_io.h"

#include "core/elements.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>

namespace sbox::io {
namespace {

constexpr double kAngstromToBohr = 1.8897259886;

const std::set<std::string> kStandardResidues = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "ASX","GLX","SEC","PYL","ACE","NME","HOH","WAT"
};

std::string trim(const std::string& value) {
    const std::size_t begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) {
        return "";
    }
    const std::size_t end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
}

std::string upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::toupper(ch));
    });
    return value;
}

std::string canonical_symbol(std::string value) {
    value = trim(value);
    if (value.empty()) {
        return value;
    }
    value = upper(value);
    if (value.size() > 2) {
        value.resize(2);
    }
    if (value.size() == 2) {
        value[1] = static_cast<char>(std::tolower(static_cast<unsigned char>(value[1])));
    }
    return value;
}

int symbol_to_z(const std::string& symbol) {
    if (trim(symbol).empty()) {
        return 0;
    }
    try {
        return sbox::elements::get_element(canonical_symbol(symbol)).Z;
    } catch (const std::exception&) {
        return 0;
    }
}

std::string infer_element_from_atom_name(const std::string& atom_name, const std::string& residue_name) {
    std::string name = trim(atom_name);
    if (name.empty()) {
        return "";
    }

    if (!name.empty() && std::isdigit(static_cast<unsigned char>(name.front()))) {
        name.erase(name.begin());
    }
    if (name.empty()) {
        return "";
    }

    const std::string residue = upper(trim(residue_name));
    const bool standard_residue = kStandardResidues.count(residue) > 0;

    if (standard_residue && upper(name) == "CA") {
        return "C";
    }

    if (name.size() >= 2) {
        const std::string first_two = canonical_symbol(name.substr(0, 2));
        if (symbol_to_z(first_two) != 0 && !standard_residue) {
            return first_two;
        }
        const std::string one = canonical_symbol(name.substr(0, 1));
        if (symbol_to_z(one) != 0) {
            return one;
        }
        if (symbol_to_z(first_two) != 0) {
            return first_two;
        }
    }

    return canonical_symbol(name.substr(0, 1));
}

int parse_int_field(const std::string& line, std::size_t start, std::size_t length) {
    if (start >= line.size()) {
        return 0;
    }
    return std::stoi(trim(line.substr(start, std::min(length, line.size() - start))));
}

double parse_double_field(const std::string& line, std::size_t start, std::size_t length) {
    if (start >= line.size()) {
        return 0.0;
    }
    const std::string token = trim(line.substr(start, std::min(length, line.size() - start)));
    if (token.empty()) {
        return 0.0;
    }
    return std::stod(token);
}

std::string field(const std::string& line, std::size_t start, std::size_t length) {
    if (start >= line.size()) {
        return "";
    }
    return trim(line.substr(start, std::min(length, line.size() - start)));
}

}  // namespace

sbox::chem::MolecularSystem PDBData::to_molecular_system() const {
    sbox::chem::MolecularSystem mol;
    mol.set_name(title);

    std::unordered_map<int, int> serial_to_index;
    serial_to_index.reserve(atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const PDBAtom& atom = atoms[i];
        mol.add_atom({atom.Z,
                      atom.position * kAngstromToBohr,
                      atom.name,
                      0});
        serial_to_index[atom.serial] = static_cast<int>(i);
    }

    bool has_explicit_bonds = false;
    for (const auto& [a_serial, b_serial] : conect_bonds) {
        const auto a_it = serial_to_index.find(a_serial);
        const auto b_it = serial_to_index.find(b_serial);
        if (a_it == serial_to_index.end() || b_it == serial_to_index.end()) {
            continue;
        }
        if (!mol.has_bond(a_it->second, b_it->second)) {
            mol.add_bond(a_it->second, b_it->second, sbox::chem::BondOrder::Single);
            has_explicit_bonds = true;
        }
    }

    if (!has_explicit_bonds) {
        mol.perceive_bonds();
    }
    return mol;
}

PDBData read_pdb(const std::string& filepath) {
    std::ifstream input(filepath);
    if (!input) {
        throw std::runtime_error("Could not open PDB file: " + filepath);
    }

    PDBData data;
    std::string line;
    std::map<std::tuple<std::string, int, std::string>, int> residue_lookup;
    std::map<std::string, int> chain_lookup;

    while (std::getline(input, line)) {
        if (line.size() < 6) {
            continue;
        }
        const std::string record = trim(line.substr(0, 6));
        if (record == "TITLE") {
            const std::string piece = field(line, 10, line.size() > 10 ? line.size() - 10 : 0);
            if (!piece.empty()) {
                if (!data.title.empty()) {
                    data.title += " ";
                }
                data.title += piece;
            }
            continue;
        }

        if (record == "ATOM" || record == "HETATM") {
            PDBAtom atom;
            atom.serial = parse_int_field(line, 6, 5);
            atom.name = field(line, 12, 4);
            atom.alt_loc = field(line, 16, 1);
            atom.residue_name = field(line, 17, 3);
            atom.chain_id = field(line, 21, 1);
            atom.residue_seq = parse_int_field(line, 22, 4);
            atom.position.x() = parse_double_field(line, 30, 8);
            atom.position.y() = parse_double_field(line, 38, 8);
            atom.position.z() = parse_double_field(line, 46, 8);
            atom.occupancy = parse_double_field(line, 54, 6);
            atom.b_factor = parse_double_field(line, 60, 6);
            atom.element = canonical_symbol(field(line, 76, 2));
            if (atom.element.empty()) {
                atom.element = infer_element_from_atom_name(atom.name, atom.residue_name);
            }
            atom.Z = symbol_to_z(atom.element);
            if (atom.Z == 0) {
                throw std::runtime_error("Unknown element in PDB atom record: " + atom.name);
            }

            const int atom_index = static_cast<int>(data.atoms.size());
            data.atoms.push_back(atom);

            const auto residue_key = std::make_tuple(atom.chain_id, atom.residue_seq, atom.residue_name);
            int residue_index = -1;
            auto residue_it = residue_lookup.find(residue_key);
            if (residue_it == residue_lookup.end()) {
                residue_index = static_cast<int>(data.residues.size());
                data.residues.push_back({atom.residue_name, atom.residue_seq, atom.chain_id, {}});
                residue_lookup[residue_key] = residue_index;
            } else {
                residue_index = residue_it->second;
            }
            data.residues[static_cast<std::size_t>(residue_index)].atom_indices.push_back(atom_index);

            int chain_index = -1;
            auto chain_it = chain_lookup.find(atom.chain_id);
            if (chain_it == chain_lookup.end()) {
                chain_index = static_cast<int>(data.chains.size());
                data.chains.push_back({atom.chain_id, {}});
                chain_lookup[atom.chain_id] = chain_index;
            } else {
                chain_index = chain_it->second;
            }
            std::vector<int>& residue_indices = data.chains[static_cast<std::size_t>(chain_index)].residue_indices;
            if (residue_indices.empty() || residue_indices.back() != residue_index) {
                residue_indices.push_back(residue_index);
            }
            continue;
        }

        if (record == "CONECT") {
            const int serial = parse_int_field(line, 6, 5);
            const std::array<std::size_t, 4> starts = {11, 16, 21, 26};
            for (std::size_t start : starts) {
                const int bonded = parse_int_field(line, start, 5);
                if (bonded > 0 && serial > 0 && bonded != serial) {
                    data.conect_bonds.emplace_back(serial, bonded);
                }
            }
        }
    }

    if (data.title.empty()) {
        data.title = "PDB Structure";
    }
    return data;
}

void write_pdb(const std::string& filepath, const PDBData& data) {
    std::ofstream output(filepath);
    if (!output) {
        throw std::runtime_error("Could not open PDB file for writing: " + filepath);
    }

    if (!data.title.empty()) {
        output << "TITLE     " << data.title << '\n';
    }

    for (const PDBAtom& atom : data.atoms) {
        const std::string record = (upper(atom.residue_name) == "HOH" || upper(atom.residue_name) == "WAT") ? "HETATM" : "ATOM  ";
        output << std::left << std::setw(6) << record
               << std::right << std::setw(5) << atom.serial << ' '
               << std::setw(4) << atom.name
               << std::setw(1) << (atom.alt_loc.empty() ? " " : atom.alt_loc)
               << std::setw(3) << atom.residue_name << ' '
               << std::setw(1) << atom.chain_id
               << std::setw(4) << atom.residue_seq
               << "    "
               << std::fixed << std::setprecision(3)
               << std::setw(8) << atom.position.x()
               << std::setw(8) << atom.position.y()
               << std::setw(8) << atom.position.z()
               << std::setprecision(2)
               << std::setw(6) << atom.occupancy
               << std::setw(6) << atom.b_factor
               << "          "
               << std::setw(2) << atom.element
               << '\n';
    }

    std::set<std::pair<int, int>> unique_bonds;
    for (const auto& bond : data.conect_bonds) {
        const int a = std::min(bond.first, bond.second);
        const int b = std::max(bond.first, bond.second);
        unique_bonds.insert({a, b});
    }

    std::map<int, std::vector<int>> adjacency;
    for (const auto& [a, b] : unique_bonds) {
        adjacency[a].push_back(b);
        adjacency[b].push_back(a);
    }

    for (auto& [serial, neighbors] : adjacency) {
        std::sort(neighbors.begin(), neighbors.end());
        output << std::left << std::setw(6) << "CONECT" << std::right << std::setw(5) << serial;
        for (int bonded : neighbors) {
            output << std::setw(5) << bonded;
        }
        output << '\n';
    }

    output << "END\n";
}

}  // namespace sbox::io

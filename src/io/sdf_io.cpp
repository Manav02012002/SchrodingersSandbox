#include "io/sdf_io.h"

#include "core/elements.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace sbox::io {
namespace {

constexpr double kAngstromToBohr = 1.8897259886;
constexpr double kBohrToAngstrom = 1.0 / kAngstromToBohr;

std::string to_lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

int symbol_to_z(const std::string& symbol) {
    for (int z = 1; z <= 118; ++z) {
        const char* candidate = sbox::elements::get_element(z).symbol;
        if (symbol == candidate) {
            return z;
        }
    }

    const std::string lowered = to_lower(symbol);
    for (int z = 1; z <= 118; ++z) {
        const char* candidate = sbox::elements::get_element(z).symbol;
        if (lowered == to_lower(candidate)) {
            return z;
        }
    }

    throw std::runtime_error("Unknown element symbol in SDF: " + symbol);
}

sbox::chem::BondOrder bond_order_from_int(int bond_type) {
    switch (bond_type) {
        case 1: return sbox::chem::BondOrder::Single;
        case 2: return sbox::chem::BondOrder::Double;
        case 3: return sbox::chem::BondOrder::Triple;
        case 4: return sbox::chem::BondOrder::Aromatic;
        default:
            throw std::runtime_error("Unsupported bond type in SDF: " + std::to_string(bond_type));
    }
}

int bond_order_to_int(sbox::chem::BondOrder order) {
    switch (order) {
        case sbox::chem::BondOrder::Single: return 1;
        case sbox::chem::BondOrder::Double: return 2;
        case sbox::chem::BondOrder::Triple: return 3;
        case sbox::chem::BondOrder::Aromatic: return 4;
        case sbox::chem::BondOrder::Unknown: return 1;
    }
    return 1;
}

sbox::chem::MolecularSystem parse_sdf_record(const std::vector<std::string>& lines) {
    if (lines.size() < 4) {
        throw std::runtime_error("Malformed SDF: missing header or counts line");
    }

    sbox::chem::MolecularSystem mol;
    mol.set_name(lines[0]);

    int num_atoms = 0;
    int num_bonds = 0;
    {
        std::istringstream counts(lines[3]);
        std::string version;
        if (!(counts >> num_atoms >> num_bonds)) {
            throw std::runtime_error("Malformed SDF: invalid counts line");
        }

        while (counts >> version) {
        }

        if (version != "V2000") {
            throw std::runtime_error("Malformed SDF: only V2000 is supported");
        }
    }

    const std::size_t required_lines =
        4U + static_cast<std::size_t>(num_atoms) + static_cast<std::size_t>(num_bonds);
    if (lines.size() < required_lines) {
        throw std::runtime_error("Malformed SDF: atom/bond block shorter than counts line declares");
    }

    std::size_t line_index = 4;
    for (int i = 0; i < num_atoms; ++i, ++line_index) {
        std::istringstream atom_stream(lines[line_index]);
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        std::string symbol;
        if (!(atom_stream >> x >> y >> z >> symbol)) {
            throw std::runtime_error("Malformed SDF: invalid atom line");
        }
        mol.add_atom({
            symbol_to_z(symbol),
            Eigen::Vector3d(x * kAngstromToBohr, y * kAngstromToBohr, z * kAngstromToBohr),
            "",
            0,
        });
    }

    for (int i = 0; i < num_bonds; ++i, ++line_index) {
        std::istringstream bond_stream(lines[line_index]);
        int atom1 = 0;
        int atom2 = 0;
        int bond_type = 0;
        if (!(bond_stream >> atom1 >> atom2 >> bond_type)) {
            throw std::runtime_error("Malformed SDF: invalid bond line");
        }
        mol.add_bond(atom1 - 1, atom2 - 1, bond_order_from_int(bond_type));
    }

    bool saw_end = false;
    for (; line_index < lines.size(); ++line_index) {
        if (lines[line_index] == "M  END") {
            saw_end = true;
            break;
        }
    }

    if (!saw_end) {
        throw std::runtime_error("Malformed SDF: missing M  END");
    }

    return mol;
}

std::vector<sbox::chem::MolecularSystem> parse_sdf_records(std::istream& input) {
    std::vector<sbox::chem::MolecularSystem> molecules;
    std::vector<std::string> current_record;
    std::string line;

    while (std::getline(input, line)) {
        if (line == "$$$$") {
            if (!current_record.empty()) {
                molecules.push_back(parse_sdf_record(current_record));
                current_record.clear();
            }
            continue;
        }
        current_record.push_back(line);
    }

    if (!current_record.empty()) {
        molecules.push_back(parse_sdf_record(current_record));
    }

    return molecules;
}

}  // namespace

sbox::chem::MolecularSystem read_sdf(const std::string& filepath) {
    std::ifstream input(filepath);
    if (!input) {
        throw std::runtime_error("Could not open SDF file: " + filepath);
    }

    const std::vector<sbox::chem::MolecularSystem> molecules = parse_sdf_records(input);
    if (molecules.empty()) {
        throw std::runtime_error("Malformed SDF: no molecules found in " + filepath);
    }
    return molecules.front();
}

sbox::chem::MolecularSystem read_sdf_string(const std::string& content) {
    std::istringstream input(content);
    const std::vector<sbox::chem::MolecularSystem> molecules = parse_sdf_records(input);
    if (molecules.empty()) {
        throw std::runtime_error("Malformed SDF: no molecules found in input string");
    }
    return molecules.front();
}

std::vector<sbox::chem::MolecularSystem> read_sdf_multi(const std::string& filepath) {
    std::ifstream input(filepath);
    if (!input) {
        throw std::runtime_error("Could not open SDF file: " + filepath);
    }
    return parse_sdf_records(input);
}

void write_sdf(const std::string& filepath, const sbox::chem::MolecularSystem& mol) {
    std::ofstream output(filepath);
    if (!output) {
        throw std::runtime_error("Could not open SDF file for writing: " + filepath);
    }

    output << mol.name() << '\n';
    output << '\n';
    output << '\n';
    output << std::setw(3) << mol.num_atoms()
           << std::setw(3) << mol.num_bonds()
           << "  0  0  0  0  0  0  0999 V2000\n";

    output << std::fixed << std::setprecision(4);
    for (const sbox::chem::Atom& atom : mol.atoms()) {
        const char* symbol = sbox::elements::get_element(atom.Z).symbol;
        output << std::setw(10) << atom.position.x() * kBohrToAngstrom
               << std::setw(10) << atom.position.y() * kBohrToAngstrom
               << std::setw(10) << atom.position.z() * kBohrToAngstrom
               << ' ' << std::left << std::setw(3) << symbol << std::right
               << " 0  0  0  0  0  0  0  0  0  0  0  0\n";
    }

    for (const sbox::chem::Bond& bond : mol.bonds()) {
        output << std::setw(3) << (bond.atom_i + 1)
               << std::setw(3) << (bond.atom_j + 1)
               << std::setw(3) << bond_order_to_int(bond.order)
               << "  0  0  0  0\n";
    }

    output << "M  END\n";
}

}  // namespace sbox::io

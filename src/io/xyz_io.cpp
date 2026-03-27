#include "io/xyz_io.h"

#include "core/elements.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

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

    throw std::runtime_error("Unknown element symbol in XYZ: " + symbol);
}

sbox::chem::MolecularSystem read_xyz_stream(std::istream& input) {
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("Malformed XYZ: missing atom count line");
    }

    std::istringstream count_stream(line);
    int atom_count = 0;
    if (!(count_stream >> atom_count) || atom_count < 0) {
        throw std::runtime_error("Malformed XYZ: invalid atom count");
    }

    std::string comment;
    if (!std::getline(input, comment)) {
        throw std::runtime_error("Malformed XYZ: missing comment line");
    }

    sbox::chem::MolecularSystem mol;
    mol.set_name(comment);

    for (int i = 0; i < atom_count; ++i) {
        if (!std::getline(input, line)) {
            throw std::runtime_error("Malformed XYZ: atom count does not match number of coordinate lines");
        }

        std::istringstream atom_stream(line);
        std::string symbol;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        if (!(atom_stream >> symbol >> x >> y >> z)) {
            throw std::runtime_error("Malformed XYZ: invalid atom line");
        }

        mol.add_atom({
            symbol_to_z(symbol),
            Eigen::Vector3d(x * kAngstromToBohr, y * kAngstromToBohr, z * kAngstromToBohr),
            "",
            0,
        });
    }

    mol.perceive_bonds();
    return mol;
}

}  // namespace

sbox::chem::MolecularSystem read_xyz(const std::string& filepath) {
    std::ifstream input(filepath);
    if (!input) {
        throw std::runtime_error("Could not open XYZ file: " + filepath);
    }
    return read_xyz_stream(input);
}

sbox::chem::MolecularSystem read_xyz_string(const std::string& content) {
    std::istringstream input(content);
    return read_xyz_stream(input);
}

void write_xyz(const std::string& filepath, const sbox::chem::MolecularSystem& mol) {
    std::ofstream output(filepath);
    if (!output) {
        throw std::runtime_error("Could not open XYZ file for writing: " + filepath);
    }

    output << mol.num_atoms() << '\n';
    output << mol.name() << '\n';
    output << std::fixed << std::setprecision(10);

    for (const sbox::chem::Atom& atom : mol.atoms()) {
        const char* symbol = sbox::elements::get_element(atom.Z).symbol;
        output << symbol << ' '
               << atom.position.x() * kBohrToAngstrom << ' '
               << atom.position.y() * kBohrToAngstrom << ' '
               << atom.position.z() * kBohrToAngstrom << '\n';
    }
}

}  // namespace sbox::io

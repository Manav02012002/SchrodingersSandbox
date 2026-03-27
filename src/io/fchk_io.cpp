#include "io/fchk_io.h"

#include <Eigen/Core>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace sbox::io {
namespace {

std::string trim(const std::string& input) {
    const std::size_t start = input.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
        return "";
    }
    const std::size_t end = input.find_last_not_of(" \t\r\n");
    return input.substr(start, end - start + 1);
}

struct EntryHeader {
    std::string label;
    char type = '\0';
    bool is_array = false;
    int count = 0;
    std::string scalar_token;
};

EntryHeader parse_header_line(const std::string& line) {
    if (line.size() < 44) {
        throw std::runtime_error("Malformed FCHK entry header: " + line);
    }

    EntryHeader header;
    header.label = trim(line.substr(0, 43));

    std::istringstream rest(line.substr(43));
    std::string type_token;
    if (!(rest >> type_token) || type_token.size() != 1) {
        throw std::runtime_error("Malformed FCHK entry type for label: " + header.label);
    }
    header.type = type_token[0];

    std::string next_token;
    if (!(rest >> next_token)) {
        throw std::runtime_error("Malformed FCHK entry payload for label: " + header.label);
    }

    if (next_token.rfind("N=", 0) == 0) {
        header.is_array = true;
        header.count = std::stoi(next_token.substr(2));
    } else {
        header.scalar_token = next_token;
    }

    return header;
}

std::vector<int> read_int_array(std::istream& input, int count, const std::string& label) {
    std::vector<int> values;
    values.reserve(static_cast<std::size_t>(count));
    int value = 0;
    for (int i = 0; i < count; ++i) {
        if (!(input >> value)) {
            throw std::runtime_error("Malformed FCHK integer array for label: " + label);
        }
        values.push_back(value);
    }
    return values;
}

std::vector<double> read_real_array(std::istream& input, int count, const std::string& label) {
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(count));
    double value = 0.0;
    for (int i = 0; i < count; ++i) {
        if (!(input >> value)) {
            throw std::runtime_error("Malformed FCHK real array for label: " + label);
        }
        values.push_back(value);
    }
    return values;
}

void skip_array(std::istream& input, char type, int count, const std::string& label) {
    if (type == 'I') {
        (void)read_int_array(input, count, label);
        return;
    }
    if (type == 'R') {
        (void)read_real_array(input, count, label);
        return;
    }
    if (type == 'C') {
        std::string token;
        for (int i = 0; i < count; ++i) {
            if (!(input >> token)) {
                throw std::runtime_error("Malformed FCHK character array for label: " + label);
            }
        }
        return;
    }
    throw std::runtime_error("Unsupported FCHK entry type for label: " + label);
}

}  // namespace

FchkData read_fchk(const std::string& filepath) {
    std::ifstream input(filepath);
    if (!input) {
        throw std::runtime_error("Could not open FCHK file: " + filepath);
    }

    FchkData data;
    if (!std::getline(input, data.title)) {
        throw std::runtime_error("Malformed FCHK: missing title line");
    }

    std::string method_line;
    if (!std::getline(input, method_line)) {
        throw std::runtime_error("Malformed FCHK: missing method line");
    }
    {
        std::istringstream iss(method_line);
        iss >> data.method >> data.basis_name;
    }

    std::string line;
    while (std::getline(input, line)) {
        if (trim(line).empty()) {
            continue;
        }

        const EntryHeader header = parse_header_line(line);

        if (header.is_array) {
            if (header.label == "Atomic numbers") {
                data.atomic_numbers = read_int_array(input, header.count, header.label);
            } else if (header.label == "Shell types") {
                data.shell_types = read_int_array(input, header.count, header.label);
            } else if (header.label == "Shell to atom map") {
                data.shell_to_atom_map = read_int_array(input, header.count, header.label);
            } else if (header.label == "Number of primitives per shell") {
                data.primitives_per_shell = read_int_array(input, header.count, header.label);
            } else if (header.label == "Current cartesian coordinates") {
                data.coordinates = read_real_array(input, header.count, header.label);
            } else if (header.label == "Primitive exponents") {
                data.primitive_exponents = read_real_array(input, header.count, header.label);
            } else if (header.label == "Contraction coefficients") {
                data.contraction_coefficients = read_real_array(input, header.count, header.label);
            } else if (header.label == "P(S=P) Contraction coefficients") {
                data.sp_contraction_coefficients = read_real_array(input, header.count, header.label);
            } else if (header.label == "Alpha Orbital Energies") {
                const std::vector<double> values = read_real_array(input, header.count, header.label);
                data.num_mo = static_cast<int>(values.size());
                data.mo_energies = Eigen::Map<const Eigen::VectorXd>(values.data(), data.num_mo);
            } else if (header.label == "Alpha MO coefficients") {
                const std::vector<double> values = read_real_array(input, header.count, header.label);
                if (data.num_basis > 0 && data.num_mo > 0 &&
                    static_cast<int>(values.size()) == data.num_basis * data.num_mo) {
                    data.mo_coefficients =
                        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(
                            values.data(), data.num_basis, data.num_mo);
                }
            } else if (header.label == "Dipole Moment") {
                const std::vector<double> values = read_real_array(input, header.count, header.label);
                if (values.size() >= 3U) {
                    data.dipole_moment = Eigen::Vector3d(values[0], values[1], values[2]);
                }
            } else if (header.label == "Mulliken Charges") {
                data.mulliken_charges = read_real_array(input, header.count, header.label);
            } else {
                skip_array(input, header.type, header.count, header.label);
            }
            continue;
        }

        if (header.label == "Number of atoms" && header.type == 'I') {
            data.num_atoms = std::stoi(header.scalar_token);
        } else if (header.label == "Charge" && header.type == 'I') {
            data.charge = std::stoi(header.scalar_token);
        } else if (header.label == "Multiplicity" && header.type == 'I') {
            data.multiplicity = std::stoi(header.scalar_token);
        } else if (header.label == "Number of basis functions" && header.type == 'I') {
            data.num_basis = std::stoi(header.scalar_token);
        } else if (header.label == "Total Energy" && header.type == 'R') {
            data.total_energy = std::stod(header.scalar_token);
        }
    }

    if (data.num_atoms <= 0) {
        throw std::runtime_error("Malformed FCHK: no atoms found");
    }

    if (data.num_mo == 0 && data.mo_energies.size() > 0) {
        data.num_mo = static_cast<int>(data.mo_energies.size());
    }

    if (data.occupations.size() == 0 && data.num_mo > 0) {
        const int num_electrons =
            std::accumulate(data.atomic_numbers.begin(), data.atomic_numbers.end(), 0) - data.charge;
        const int num_doubly_occupied = std::max(0, num_electrons / 2);
        data.occupations = Eigen::VectorXd::Zero(data.num_mo);
        for (int i = 0; i < std::min(num_doubly_occupied, data.num_mo); ++i) {
            data.occupations(i) = 2.0;
        }
    }

    return data;
}

}  // namespace sbox::io

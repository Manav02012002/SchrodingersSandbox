#include "core/molden_parser.h"

#include <Eigen/Core>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace sbox::molden {
namespace {

constexpr double kBohrPerAngstrom = 1.8897261254578281;

std::string trim(const std::string& input) {
    const std::size_t start = input.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
        return "";
    }
    const std::size_t end = input.find_last_not_of(" \t\r\n");
    return input.substr(start, end - start + 1);
}

std::string to_upper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char ch) {
        return static_cast<char>(std::toupper(ch));
    });
    return s;
}

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return s;
}

std::vector<std::string> split_ws(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream iss(line);
    std::string token;
    while (iss >> token) {
        out.push_back(token);
    }
    return out;
}

double parse_double_token(std::string token) {
    for (char& ch : token) {
        if (ch == 'D' || ch == 'd') {
            ch = 'E';
        }
    }
    return std::stod(token);
}

int parse_int_token(const std::string& token) {
    std::size_t pos = 0;
    const int value = std::stoi(token, &pos);
    if (pos != token.size()) {
        throw std::runtime_error("Invalid integer token: " + token);
    }
    return value;
}

int shell_l_from_label(const std::string& label) {
    const std::string lower = to_lower(label);
    if (lower == "s") return 0;
    if (lower == "p") return 1;
    if (lower == "d") return 2;
    if (lower == "f") return 3;
    return -1;
}

void renormalize_shell_contractions(sbox::basis::MOData& data) {
    for (sbox::basis::BasisShell& shell : data.basis.shells) {
        double sum_sq = 0.0;
        for (const sbox::basis::GaussianPrimitive& prim : shell.primitives) {
            sum_sq += prim.coefficient * prim.coefficient;
        }
        if (sum_sq <= 0.0) {
            continue;
        }
        const double inv_norm = 1.0 / std::sqrt(sum_sq);
        for (sbox::basis::GaussianPrimitive& prim : shell.primitives) {
            prim.coefficient *= inv_norm;
        }
    }
}

void finalize_mo_block(std::vector<double>& energies,
                       std::vector<double>& occupations,
                       std::vector<std::vector<double>>& coeff_columns,
                       bool& has_active_orbital,
                       bool& has_content,
                       double current_energy,
                       double current_occup,
                       std::vector<double>& current_coeffs) {
    if (!has_active_orbital || !has_content) {
        return;
    }

    energies.push_back(current_energy);
    occupations.push_back(current_occup);
    coeff_columns.push_back(current_coeffs);

    has_active_orbital = false;
    has_content = false;
    current_energy = 0.0;
    current_occup = 0.0;
}

}  // namespace

sbox::basis::MOData parse_molden_file(const std::string& filepath) {
    return parse_molden_file(filepath, ParseOptions{});
}

sbox::basis::MOData parse_molden_file(const std::string& filepath, const ParseOptions& options) {
    std::ifstream in(filepath);
    if (!in) {
        throw std::runtime_error("Could not open molden file: " + filepath);
    }

    std::vector<std::string> lines;
    std::string raw;
    while (std::getline(in, raw)) {
        lines.push_back(raw);
    }

    sbox::basis::MOData result;
    result.total_energy = 0.0;

    enum class Section { None, Atoms, GTO, MO };
    Section section = Section::None;

    bool atoms_in_au = false;
    bool atoms_seen = false;
    bool gto_seen = false;
    bool mo_seen = false;

    std::vector<double> mo_energies;
    std::vector<double> mo_occupations;
    std::vector<std::vector<double>> mo_coeff_columns;
    bool mo_has_active_orbital = false;
    bool mo_has_content = false;
    double mo_current_energy = 0.0;
    double mo_current_occup = 0.0;
    std::vector<double> mo_current_coeffs;

    int current_gto_atom = -1;

    for (std::size_t i = 0; i < lines.size(); ++i) {
        const std::string line = trim(lines[i]);
        if (line.empty()) {
            continue;
        }

        if (line.front() == '[') {
            if (section == Section::MO) {
                finalize_mo_block(mo_energies,
                                  mo_occupations,
                                  mo_coeff_columns,
                                  mo_has_active_orbital,
                                  mo_has_content,
                                  mo_current_energy,
                                  mo_current_occup,
                                  mo_current_coeffs);
            }

            const std::size_t close = line.find(']');
            if (close == std::string::npos) {
                continue;
            }

            const std::string sec_name = to_upper(trim(line.substr(1, close - 1)));
            const std::string sec_suffix = to_upper(trim(line.substr(close + 1)));

            if (sec_name == "ATOMS") {
                section = Section::Atoms;
                atoms_seen = true;
                atoms_in_au = sec_suffix.find("AU") != std::string::npos;
            } else if (sec_name == "GTO") {
                section = Section::GTO;
                gto_seen = true;
                current_gto_atom = -1;
            } else if (sec_name == "MO") {
                section = Section::MO;
                mo_seen = true;
                mo_has_active_orbital = false;
                mo_has_content = false;
                mo_current_energy = 0.0;
                mo_current_occup = 0.0;
                mo_current_coeffs.assign(static_cast<std::size_t>(result.basis.num_basis_functions()), 0.0);
            } else {
                section = Section::None;
            }

            if (sec_name == "5D" || sec_name == "5D7F" || sec_name == "7F" || sec_name == "9G") {
                result.basis.spherical = true;
            }
            if (sec_name == "6D" || sec_name == "6D10F" || sec_name == "10F" || sec_name == "15G") {
                result.basis.spherical = false;
            }

            continue;
        }

        if (section == Section::Atoms) {
            const std::vector<std::string> tokens = split_ws(line);
            if (tokens.size() < 6) {
                continue;
            }
            const double x = parse_double_token(tokens[tokens.size() - 3]);
            const double y = parse_double_token(tokens[tokens.size() - 2]);
            const double z = parse_double_token(tokens[tokens.size() - 1]);
            const double scale = atoms_in_au ? 1.0 : kBohrPerAngstrom;
            result.atom_positions.emplace_back(x * scale, y * scale, z * scale);
            continue;
        }

        if (section == Section::GTO) {
            if (line == "****") {
                current_gto_atom = -1;
                continue;
            }

            const std::vector<std::string> tokens = split_ws(line);
            if (tokens.empty()) {
                continue;
            }

            if (current_gto_atom < 0) {
                try {
                    current_gto_atom = parse_int_token(tokens[0]) - 1;
                } catch (const std::exception&) {
                    current_gto_atom = -1;
                }
                continue;
            }

            const std::string shell_label = to_lower(tokens[0]);
            if (tokens.size() < 2) {
                continue;
            }
            const int nprim = parse_int_token(tokens[1]);
            const double shell_scale = tokens.size() >= 3 ? parse_double_token(tokens[2]) : 1.0;

            if (shell_label == "sp") {
                sbox::basis::BasisShell shell_s;
                shell_s.atom_index = current_gto_atom;
                shell_s.angular_momentum = 0;
                sbox::basis::BasisShell shell_p;
                shell_p.atom_index = current_gto_atom;
                shell_p.angular_momentum = 1;

                for (int p = 0; p < nprim && i + 1 < lines.size(); ++p) {
                    const std::vector<std::string> prim_toks = split_ws(trim(lines[++i]));
                    if (prim_toks.size() < 3) {
                        continue;
                    }
                    const double exponent = parse_double_token(prim_toks[0]);
                    shell_s.primitives.push_back({exponent, parse_double_token(prim_toks[1]) * shell_scale});
                    shell_p.primitives.push_back({exponent, parse_double_token(prim_toks[2]) * shell_scale});
                }

                if (!shell_s.primitives.empty()) {
                    result.basis.shells.push_back(shell_s);
                }
                if (!shell_p.primitives.empty()) {
                    result.basis.shells.push_back(shell_p);
                }
                continue;
            }

            const int l = shell_l_from_label(shell_label);
            if (l < 0) {
                for (int p = 0; p < nprim && i + 1 < lines.size(); ++p) {
                    ++i;
                }
                continue;
            }

            sbox::basis::BasisShell shell;
            shell.atom_index = current_gto_atom;
            shell.angular_momentum = l;
            shell.primitives.reserve(static_cast<std::size_t>(nprim));

            for (int p = 0; p < nprim && i + 1 < lines.size(); ++p) {
                const std::vector<std::string> prim_toks = split_ws(trim(lines[++i]));
                if (prim_toks.size() < 2) {
                    continue;
                }
                shell.primitives.push_back({parse_double_token(prim_toks[0]),
                                            parse_double_token(prim_toks[1]) * shell_scale});
            }

            if (!shell.primitives.empty()) {
                result.basis.shells.push_back(shell);
            }

            continue;
        }

        if (section == Section::MO) {
            if (mo_current_coeffs.empty()) {
                mo_current_coeffs.assign(static_cast<std::size_t>(result.basis.num_basis_functions()), 0.0);
            }

            const std::size_t eq_pos = line.find('=');
            if (eq_pos != std::string::npos) {
                const std::string key = to_upper(trim(line.substr(0, eq_pos)));
                const std::string value = trim(line.substr(eq_pos + 1));

                if (key == "ENE") {
                    finalize_mo_block(mo_energies,
                                      mo_occupations,
                                      mo_coeff_columns,
                                      mo_has_active_orbital,
                                      mo_has_content,
                                      mo_current_energy,
                                      mo_current_occup,
                                      mo_current_coeffs);
                    mo_has_active_orbital = true;
                    mo_has_content = true;
                    mo_current_energy = parse_double_token(value);
                    mo_current_occup = 0.0;
                    mo_current_coeffs.assign(static_cast<std::size_t>(result.basis.num_basis_functions()), 0.0);
                } else if (key == "OCCUP") {
                    if (!mo_has_active_orbital) {
                        mo_has_active_orbital = true;
                        mo_current_coeffs.assign(static_cast<std::size_t>(result.basis.num_basis_functions()), 0.0);
                    }
                    mo_has_content = true;
                    mo_current_occup = parse_double_token(value);
                } else if (key == "SPIN" || key == "SYM") {
                    if (!mo_has_active_orbital) {
                        mo_has_active_orbital = true;
                        mo_current_coeffs.assign(static_cast<std::size_t>(result.basis.num_basis_functions()), 0.0);
                    }
                } else if (key == "ENERGY") {
                    // Some producers store total SCF energy here; keep as optional metadata.
                    result.total_energy = parse_double_token(value);
                }
                continue;
            }

            const std::vector<std::string> coeff_toks = split_ws(line);
            if (coeff_toks.size() >= 2) {
                if (!mo_has_active_orbital) {
                    mo_has_active_orbital = true;
                    mo_current_coeffs.assign(static_cast<std::size_t>(result.basis.num_basis_functions()), 0.0);
                }
                mo_has_content = true;

                const int basis_index = parse_int_token(coeff_toks[0]) - 1;
                if (basis_index >= 0 && basis_index < static_cast<int>(mo_current_coeffs.size())) {
                    mo_current_coeffs[static_cast<std::size_t>(basis_index)] = parse_double_token(coeff_toks[1]);
                }
            }
        }
    }

    if (section == Section::MO) {
        finalize_mo_block(mo_energies,
                          mo_occupations,
                          mo_coeff_columns,
                          mo_has_active_orbital,
                          mo_has_content,
                          mo_current_energy,
                          mo_current_occup,
                          mo_current_coeffs);
    }

    if (!atoms_seen || !gto_seen || !mo_seen) {
        throw std::runtime_error("Molden file missing required sections [Atoms], [GTO], or [MO]: " + filepath);
    }

    if (!options.contraction_coefficients_include_shell_normalization) {
        renormalize_shell_contractions(result);
    }

    const int n_basis = result.basis.num_basis_functions();
    const int n_mo = static_cast<int>(mo_coeff_columns.size());
    result.coefficients = Eigen::MatrixXd::Zero(n_basis, n_mo);
    result.energies = Eigen::VectorXd::Zero(n_mo);
    result.occupations = Eigen::VectorXd::Zero(n_mo);

    for (int col = 0; col < n_mo; ++col) {
        result.energies(col) = mo_energies[static_cast<std::size_t>(col)];
        result.occupations(col) = mo_occupations[static_cast<std::size_t>(col)];

        const std::vector<double>& coeffs = mo_coeff_columns[static_cast<std::size_t>(col)];
        const int limit = std::min<int>(n_basis, static_cast<int>(coeffs.size()));
        for (int row = 0; row < limit; ++row) {
            result.coefficients(row, col) = coeffs[static_cast<std::size_t>(row)];
        }
    }

    return result;
}

}  // namespace sbox::molden

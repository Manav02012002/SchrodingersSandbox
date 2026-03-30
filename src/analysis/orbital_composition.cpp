#include "analysis/orbital_composition.h"

#include "core/elements.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace sbox::analysis {

namespace {

constexpr double kHartreeToEv = 27.2114;

int shell_basis_count(const sbox::basis::BasisSet& basis, int l) {
    if (basis.spherical) {
        return 2 * l + 1;
    }
    return (l + 1) * (l + 2) / 2;
}

std::string l_to_letter(int l) {
    switch (l) {
    case 0: return "s";
    case 1: return "p";
    case 2: return "d";
    case 3: return "f";
    default: return "?";
    }
}

int approximate_principal_n(int Z, int l, int shell_count_for_l) {
    if (Z <= 1 && l == 0) {
        return 1;
    }
    if (Z <= 2) {
        return l + 1;
    }
    if (Z <= 10) {
        if (l == 0) {
            return shell_count_for_l == 0 ? 1 : 2;
        }
        if (l == 1) {
            return 2;
        }
        return l + 1;
    }
    if (Z <= 18) {
        return std::max(l + 1, shell_count_for_l + 1);
    }
    return std::max(l + 1, shell_count_for_l + 1);
}

}  // namespace

OrbitalComposition analyze_orbital_composition(const sbox::basis::MOData& mo_data,
                                               const sbox::chem::MolecularSystem& mol,
                                               int mo_index) {
    OrbitalComposition composition;
    composition.mo_index = mo_index;

    if (mo_index < 0 || mo_index >= mo_data.coefficients.cols()) {
        return composition;
    }

    if (mo_index < mo_data.energies.size()) {
        composition.energy_eV = mo_data.energies(mo_index) * kHartreeToEv;
    }
    if (mo_index < mo_data.occupations.size()) {
        composition.occupation = mo_data.occupations(mo_index);
    }

    const Eigen::VectorXd coeffs = mo_data.coefficients.col(mo_index);
    const std::vector<int>& atomic_numbers = !mo_data.atomic_numbers.empty() ? mo_data.atomic_numbers : std::vector<int>{};

    struct AtomAccum {
        double total = 0.0;
        std::map<std::string, double> ao;
    };

    std::map<int, AtomAccum> atom_map;
    std::map<int, std::map<int, int>> atom_l_counts;

    int basis_offset = 0;
    double total_weight = 0.0;
    for (const sbox::basis::BasisShell& shell : mo_data.basis.shells) {
        const int atom_index = shell.atom_index;
        if (atom_index < 0 || atom_index >= mol.num_atoms()) {
            basis_offset += shell_basis_count(mo_data.basis, shell.angular_momentum);
            continue;
        }

        const int count = shell_basis_count(mo_data.basis, shell.angular_momentum);
        double shell_weight = 0.0;
        for (int i = 0; i < count && basis_offset + i < coeffs.size(); ++i) {
            const double c = coeffs(basis_offset + i);
            shell_weight += c * c;
        }

        const int shell_count_for_l = atom_l_counts[atom_index][shell.angular_momentum]++;
        const int Z = atom_index < static_cast<int>(atomic_numbers.size()) ? atomic_numbers[static_cast<std::size_t>(atom_index)] : mol.atom(atom_index).Z;
        const int n = approximate_principal_n(Z, shell.angular_momentum, shell_count_for_l);
        const std::string ao_label = std::to_string(n) + l_to_letter(shell.angular_momentum);

        atom_map[atom_index].total += shell_weight;
        atom_map[atom_index].ao[ao_label] += shell_weight;
        total_weight += shell_weight;
        basis_offset += count;
    }

    if (total_weight <= 1.0e-12) {
        return composition;
    }

    for (auto& [atom_index, accum] : atom_map) {
        AtomContribution contribution;
        contribution.atom_index = atom_index;
        contribution.Z = mol.atom(atom_index).Z;
        contribution.element = sbox::elements::get_element(contribution.Z).symbol;
        contribution.total_weight = accum.total / total_weight;

        std::vector<std::pair<std::string, double>> ao_parts;
        ao_parts.reserve(accum.ao.size());
        for (const auto& [label, weight] : accum.ao) {
            ao_parts.emplace_back(label, weight / total_weight);
        }
        std::sort(ao_parts.begin(), ao_parts.end(), [](const auto& a, const auto& b) {
            return a.second > b.second;
        });
        contribution.ao_contributions = std::move(ao_parts);
        composition.atom_contributions.push_back(std::move(contribution));
    }

    std::sort(composition.atom_contributions.begin(), composition.atom_contributions.end(), [](const AtomContribution& a, const AtomContribution& b) {
        return a.total_weight > b.total_weight;
    });

    struct SummaryPart {
        std::string text;
        double weight = 0.0;
    };
    std::vector<SummaryPart> parts;
    for (const AtomContribution& atom : composition.atom_contributions) {
        for (const auto& [label, weight] : atom.ao_contributions) {
            if (weight > 0.05) {
                std::ostringstream oss;
                oss << static_cast<int>(std::round(weight * 100.0)) << "% " << atom.element << " " << label;
                parts.push_back({oss.str(), weight});
            }
        }
    }
    std::sort(parts.begin(), parts.end(), [](const SummaryPart& a, const SummaryPart& b) {
        return a.weight > b.weight;
    });

    std::ostringstream summary;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (i > 0) {
            summary << " + ";
        }
        summary << parts[i].text;
    }
    if (parts.empty()) {
        summary << "Mixed composition";
    }
    composition.summary = summary.str();

    return composition;
}

}  // namespace sbox::analysis

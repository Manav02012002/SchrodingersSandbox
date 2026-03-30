#include "analysis/crystal_field.h"

#include "core/elements.h"

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
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

struct OrbitalSample {
    int mo_index = -1;
    double energy_ev = 0.0;
    double total_d_weight = 0.0;
    std::map<std::string, double> component_weights;
};

std::vector<std::string> spherical_d_labels() {
    return {"dz2", "dxz", "dyz", "dx2y2", "dxy"};
}

std::vector<std::string> cartesian_d_labels() {
    return {"dx2y2", "dx2y2", "dz2", "dxy", "dxz", "dyz"};
}

double get_orbital_energy(const DOrbitalEnergies& d_orbs, const std::string& label) {
    if (label == "dxy") return d_orbs.dxy;
    if (label == "dxz") return d_orbs.dxz;
    if (label == "dyz") return d_orbs.dyz;
    if (label == "dz2") return d_orbs.dz2;
    return d_orbs.dx2y2;
}

int* get_orbital_mo_index(DOrbitalEnergies& d_orbs, const std::string& label) {
    if (label == "dxy") return &d_orbs.mo_dxy;
    if (label == "dxz") return &d_orbs.mo_dxz;
    if (label == "dyz") return &d_orbs.mo_dyz;
    if (label == "dz2") return &d_orbs.mo_dz2;
    return &d_orbs.mo_dx2y2;
}

double average_group(const DOrbitalEnergies& d_orbs, const std::vector<std::string>& labels) {
    double sum = 0.0;
    for (const std::string& label : labels) {
        sum += get_orbital_energy(d_orbs, label);
    }
    return labels.empty() ? 0.0 : sum / static_cast<double>(labels.size());
}

bool is_transition_metal_Z(int Z) {
    return (Z >= 21 && Z <= 30) || (Z >= 39 && Z <= 48) || (Z >= 57 && Z <= 80) || (Z >= 89 && Z <= 112);
}

}  // namespace

double DOrbitalEnergies::delta_oct() const {
    return average_group(*this, {"dz2", "dx2y2"}) - average_group(*this, {"dxy", "dxz", "dyz"});
}

double DOrbitalEnergies::delta_tet() const {
    return average_group(*this, {"dxy", "dxz", "dyz"}) - average_group(*this, {"dz2", "dx2y2"});
}

double DOrbitalEnergies::mean_energy() const {
    return (dxy + dxz + dyz + dz2 + dx2y2) / 5.0;
}

DOrbitalEnergies extract_d_orbitals(const sbox::basis::MOData& mo_data,
                                    const sbox::chem::MolecularSystem& mol,
                                    int metal_atom_index) {
    DOrbitalEnergies out;
    if (metal_atom_index < 0 || metal_atom_index >= mol.num_atoms() || !is_transition_metal_Z(mol.atom(metal_atom_index).Z) ||
        mo_data.coefficients.cols() == 0) {
        return out;
    }

    const std::vector<std::string> labels = mo_data.basis.spherical ? spherical_d_labels() : cartesian_d_labels();
    std::vector<OrbitalSample> samples;
    samples.reserve(static_cast<std::size_t>(mo_data.coefficients.cols()));

    for (int mo = 0; mo < mo_data.coefficients.cols(); ++mo) {
        OrbitalSample sample;
        sample.mo_index = mo;
        sample.energy_ev = mo < mo_data.energies.size() ? mo_data.energies(mo) * kHartreeToEv : 0.0;

        int basis_offset = 0;
        for (const auto& shell : mo_data.basis.shells) {
            const int count = shell_basis_count(mo_data.basis, shell.angular_momentum);
            if (shell.atom_index == metal_atom_index && shell.angular_momentum == 2) {
                for (int i = 0; i < count && i < static_cast<int>(labels.size()) && basis_offset + i < mo_data.coefficients.rows(); ++i) {
                    const double coeff = mo_data.coefficients(basis_offset + i, mo);
                    const double w = coeff * coeff;
                    sample.total_d_weight += w;
                    sample.component_weights[labels[static_cast<std::size_t>(i)]] += w;
                }
            }
            basis_offset += count;
        }

        if (sample.total_d_weight > 0.0) {
            samples.push_back(std::move(sample));
        }
    }

    std::sort(samples.begin(), samples.end(), [](const OrbitalSample& a, const OrbitalSample& b) {
        return a.total_d_weight > b.total_d_weight;
    });
    if (samples.size() > 5) {
        samples.resize(5);
    }

    std::map<std::string, double*> target_map = {
        {"dxy", &out.dxy},
        {"dxz", &out.dxz},
        {"dyz", &out.dyz},
        {"dz2", &out.dz2},
        {"dx2y2", &out.dx2y2},
    };

    std::map<std::string, bool> assigned;
    for (const auto& [label, ptr] : target_map) {
        (void)ptr;
        assigned[label] = false;
    }

    for (const OrbitalSample& sample : samples) {
        std::vector<std::pair<std::string, double>> comps(sample.component_weights.begin(), sample.component_weights.end());
        std::sort(comps.begin(), comps.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
        for (const auto& [label, weight] : comps) {
            (void)weight;
            auto it = target_map.find(label);
            if (it != target_map.end() && !assigned[label]) {
                *it->second = sample.energy_ev;
                *get_orbital_mo_index(out, label) = sample.mo_index;
                assigned[label] = true;
                break;
            }
        }
    }

    std::vector<double> fallback_energies;
    for (const OrbitalSample& sample : samples) {
        fallback_energies.push_back(sample.energy_ev);
    }
    std::sort(fallback_energies.begin(), fallback_energies.end());
    std::size_t fallback_idx = 0;
    for (auto& [label, ptr] : target_map) {
        if (!assigned[label] && fallback_idx < fallback_energies.size()) {
            *ptr = fallback_energies[fallback_idx++];
        }
    }

    return out;
}

void identify_splitting(DOrbitalEnergies& d_orbs, sbox::chem::CoordinationGeometry geom) {
    d_orbs.groups.clear();
    d_orbs.ordering_warning = false;

    auto add_group = [&](const std::string& label, const std::vector<std::string>& orbitals) {
        DOrbitalEnergies::OrbitalGroup group;
        group.label = label;
        group.orbitals = orbitals;
        group.average_energy = average_group(d_orbs, orbitals);
        group.degeneracy = static_cast<int>(orbitals.size());
        d_orbs.groups.push_back(group);
    };

    switch (geom) {
    case sbox::chem::CoordinationGeometry::Octahedral:
        add_group("t2g", {"dxy", "dxz", "dyz"});
        add_group("e_g", {"dz2", "dx2y2"});
        if (d_orbs.groups[0].average_energy > d_orbs.groups[1].average_energy) {
            d_orbs.ordering_warning = true;
        }
        break;
    case sbox::chem::CoordinationGeometry::Tetrahedral:
        add_group("e", {"dz2", "dx2y2"});
        add_group("t2", {"dxy", "dxz", "dyz"});
        if (d_orbs.groups[0].average_energy > d_orbs.groups[1].average_energy) {
            d_orbs.ordering_warning = true;
        }
        break;
    case sbox::chem::CoordinationGeometry::SquarePlanar:
        add_group("dxz/dyz", {"dxz", "dyz"});
        add_group("dz2", {"dz2"});
        add_group("dxy", {"dxy"});
        add_group("dx2y2", {"dx2y2"});
        for (std::size_t i = 1; i < d_orbs.groups.size(); ++i) {
            if (d_orbs.groups[i - 1].average_energy > d_orbs.groups[i].average_energy) {
                d_orbs.ordering_warning = true;
                break;
            }
        }
        break;
    default:
        add_group("dxy", {"dxy"});
        add_group("dxz", {"dxz"});
        add_group("dyz", {"dyz"});
        add_group("dz2", {"dz2"});
        add_group("dx2y2", {"dx2y2"});
        break;
    }
}

int d_electron_count(int metal_Z, int oxidation_state) {
    int core = 18;
    if (metal_Z >= 39 && metal_Z <= 48) {
        core = 36;
    } else if (metal_Z >= 57 && metal_Z <= 80) {
        core = 54;
    } else if (metal_Z >= 89) {
        core = 86;
    }
    return std::max(0, metal_Z - core - oxidation_state);
}

double octahedral_cfse_dq(int d_electrons, bool high_spin) {
    d_electrons = std::clamp(d_electrons, 0, 10);
    int t2g = 0;
    int eg = 0;
    if (high_spin) {
        if (d_electrons <= 3) {
            t2g = d_electrons;
        } else if (d_electrons <= 5) {
            t2g = 3;
            eg = d_electrons - 3;
        } else if (d_electrons <= 8) {
            t2g = 3 + (d_electrons - 5);
            eg = 2;
        } else {
            t2g = 6;
            eg = d_electrons - 6;
        }
    } else {
        if (d_electrons <= 6) {
            t2g = d_electrons;
        } else {
            t2g = 6;
            eg = d_electrons - 6;
        }
    }
    return -4.0 * static_cast<double>(t2g) + 6.0 * static_cast<double>(eg);
}

}  // namespace sbox::analysis

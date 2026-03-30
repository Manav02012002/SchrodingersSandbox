#include "ui/bond_order_panel.h"

#include "core/elements.h"

#include <imgui.h>
#include <implot.h>

#include <algorithm>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

namespace sbox::ui {

namespace {

struct BondEntry {
    int atom_i = -1;
    int atom_j = -1;
    double order = 0.0;
};

const char* classify_bond_order(double order) {
    if (order < 0.5) {
        return "weak/nonbonding";
    }
    if (order < 1.3) {
        return "single";
    }
    if (order < 1.8) {
        return "aromatic / partial double";
    }
    if (order < 2.3) {
        return "double";
    }
    if (order < 2.8) {
        return "partial triple";
    }
    return "triple";
}

ImVec4 classification_color(double order) {
    if (order < 0.5) {
        return ImVec4(0.60f, 0.60f, 0.62f, 1.0f);
    }
    if (order < 1.3) {
        return ImVec4(0.80f, 0.80f, 0.84f, 1.0f);
    }
    if (order < 1.8) {
        return ImVec4(0.95f, 0.66f, 0.24f, 1.0f);
    }
    if (order < 2.3) {
        return ImVec4(0.24f, 0.78f, 0.36f, 1.0f);
    }
    if (order < 2.8) {
        return ImVec4(0.42f, 0.64f, 0.95f, 1.0f);
    }
    return ImVec4(0.55f, 0.42f, 0.95f, 1.0f);
}

double topological_numeric_order(sbox::chem::BondOrder order) {
    switch (order) {
    case sbox::chem::BondOrder::Single: return 1.0;
    case sbox::chem::BondOrder::Double: return 2.0;
    case sbox::chem::BondOrder::Triple: return 3.0;
    case sbox::chem::BondOrder::Aromatic: return 1.5;
    case sbox::chem::BondOrder::Unknown: return 0.0;
    }
    return 0.0;
}

const char* topological_name(sbox::chem::BondOrder order) {
    switch (order) {
    case sbox::chem::BondOrder::Single: return "Single";
    case sbox::chem::BondOrder::Double: return "Double";
    case sbox::chem::BondOrder::Triple: return "Triple";
    case sbox::chem::BondOrder::Aromatic: return "Aromatic";
    case sbox::chem::BondOrder::Unknown: return "Unknown";
    }
    return "Unknown";
}

std::string atom_ref(const sbox::chem::MolecularSystem& mol, int atom_index) {
    return std::string(sbox::elements::get_element(mol.atom(atom_index).Z).symbol) + std::to_string(atom_index + 1);
}

std::vector<BondEntry> collect_bond_orders(const sbox::backend::JobResult& result, const sbox::chem::MolecularSystem& mol) {
    std::vector<BondEntry> bonds;
    if (result.mayer_bond_orders.rows() == 0 || result.mayer_bond_orders.cols() == 0) {
        return bonds;
    }
    for (int i = 0; i < mol.num_atoms(); ++i) {
        for (int j = i + 1; j < mol.num_atoms(); ++j) {
            if (i < result.mayer_bond_orders.rows() && j < result.mayer_bond_orders.cols()) {
                const double order = result.mayer_bond_orders(i, j);
                if (order > 0.3) {
                    bonds.push_back({i, j, order});
                }
            }
        }
    }
    std::sort(bonds.begin(), bonds.end(), [](const BondEntry& a, const BondEntry& b) {
        return a.order > b.order;
    });
    return bonds;
}

}  // namespace

void draw_bond_order_panel(AppState& state,
                           const sbox::backend::JobResult& result,
                           const sbox::chem::MolecularSystem& mol) {
    if (!result.converged() || result.mayer_bond_orders.rows() == 0 || result.mayer_bond_orders.cols() == 0) {
        return;
    }

    if (!ImGui::Begin("Bond Orders")) {
        ImGui::End();
        return;
    }

    if (ImGui::Button("<- Back to Dashboard")) {
        state.property_view = PropertyView::Dashboard;
    }
    ImGui::Separator();

    const std::vector<BondEntry> bonds = collect_bond_orders(result, mol);

    if (ImGui::BeginTable("BondOrderTable", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("Bond", ImGuiTableColumnFlags_WidthFixed, 48.0f);
        ImGui::TableSetupColumn("Atoms", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Mayer Order", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Type", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();

        for (std::size_t idx = 0; idx < bonds.size(); ++idx) {
            const BondEntry& bond = bonds[idx];
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("B%zu", idx + 1);
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%s-%s", atom_ref(mol, bond.atom_i).c_str(), atom_ref(mol, bond.atom_j).c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%.2f", bond.order);
            ImGui::TableSetColumnIndex(3);
            ImGui::TextColored(classification_color(bond.order), "%s", classify_bond_order(bond.order));
        }

        ImGui::EndTable();
    }

    ImGui::Separator();
    const int n = mol.num_atoms();
    if (n <= 30) {
        std::vector<float> data(static_cast<std::size_t>(n * n), 0.0f);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                data[static_cast<std::size_t>(i * n + j)] = static_cast<float>(result.mayer_bond_orders(i, j));
            }
        }

        if (ImPlot::BeginPlot("Bond Order Matrix", ImVec2(-1, 300))) {
            ImPlot::SetupAxes("Atom i", "Atom j", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::PlotHeatmap("##bo_heatmap",
                                data.data(),
                                n,
                                n,
                                0.0,
                                3.0,
                                nullptr,
                                ImPlotPoint(0, 0),
                                ImPlotPoint(n, n));
            ImPlot::EndPlot();
        }
    } else {
        ImGui::TextUnformatted("Heatmap not shown for molecules with > 30 atoms.");
    }

    ImGui::Separator();
    int total_bonds = 0;
    double sum_orders = 0.0;
    BondEntry strongest{};
    BondEntry weakest{};
    bool have_strong = false;
    bool have_weak = false;
    for (const BondEntry& bond : bonds) {
        if (bond.order > 0.5) {
            ++total_bonds;
            sum_orders += bond.order;
            if (!have_strong || bond.order > strongest.order) {
                strongest = bond;
                have_strong = true;
            }
            if (!have_weak || bond.order < weakest.order) {
                weakest = bond;
                have_weak = true;
            }
        }
    }
    ImGui::Text("Total bonds (order > 0.5): %d", total_bonds);
    ImGui::Text("Average bond order: %.2f", total_bonds > 0 ? sum_orders / static_cast<double>(total_bonds) : 0.0);
    if (have_strong) {
        ImGui::Text("Strongest bond: %s-%s = %.2f",
                    atom_ref(mol, strongest.atom_i).c_str(),
                    atom_ref(mol, strongest.atom_j).c_str(),
                    strongest.order);
    }
    if (have_weak) {
        ImGui::Text("Weakest bond (>0.5): %s-%s = %.2f",
                    atom_ref(mol, weakest.atom_i).c_str(),
                    atom_ref(mol, weakest.atom_j).c_str(),
                    weakest.order);
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Topological Comparison");
    bool any_warning = false;
    for (const sbox::chem::Bond& bond : mol.bonds()) {
        if (bond.atom_i >= result.mayer_bond_orders.rows() || bond.atom_j >= result.mayer_bond_orders.cols()) {
            continue;
        }
        const double mayer = result.mayer_bond_orders(bond.atom_i, bond.atom_j);
        const double topo = topological_numeric_order(bond.order);
        if ((bond.order == sbox::chem::BondOrder::Single && mayer > 1.5) ||
            (bond.order == sbox::chem::BondOrder::Double && mayer < 1.4) ||
            (bond.order == sbox::chem::BondOrder::Triple && mayer < 2.4) ||
            (bond.order == sbox::chem::BondOrder::Aromatic && (mayer < 1.1 || mayer > 1.9)) ||
            (bond.order == sbox::chem::BondOrder::Unknown && mayer > 0.8) ||
            (topo > 0.0 && std::abs(mayer - topo) > 0.6)) {
            any_warning = true;
            ImGui::TextColored(ImVec4(0.98f, 0.72f, 0.24f, 1.0f),
                               "Warning Bond %s-%s: topological=%s, Mayer=%.2f",
                               atom_ref(mol, bond.atom_i).c_str(),
                               atom_ref(mol, bond.atom_j).c_str(),
                               topological_name(bond.order),
                               mayer);
        }
    }
    if (!any_warning) {
        ImGui::TextDisabled("No significant discrepancies between topological and Mayer bond orders.");
    }

    ImGui::End();
}

}  // namespace sbox::ui

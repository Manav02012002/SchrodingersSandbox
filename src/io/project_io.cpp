#include "io/project_io.h"

#include <fstream>
#include <stdexcept>
#include <string>

namespace sbox::io {

void save_project(const std::string& filepath,
                  const sbox::chem::MolecularSystem& mol,
                  const nlohmann::json& extra_state) {
    nlohmann::json root;
    root["version"] = 1;
    root["name"] = mol.name();
    root["charge"] = mol.charge();
    root["multiplicity"] = mol.multiplicity();

    root["atoms"] = nlohmann::json::array();
    for (const sbox::chem::Atom& atom : mol.atoms()) {
        root["atoms"].push_back({
            {"Z", atom.Z},
            {"position", {atom.position.x(), atom.position.y(), atom.position.z()}},
            {"label", atom.label},
            {"formal_charge", atom.formal_charge},
        });
    }

    root["bonds"] = nlohmann::json::array();
    for (const sbox::chem::Bond& bond : mol.bonds()) {
        root["bonds"].push_back({
            {"i", bond.atom_i},
            {"j", bond.atom_j},
            {"order", static_cast<int>(bond.order)},
        });
    }

    root["state"] = extra_state.is_object() ? extra_state : nlohmann::json::object();

    std::ofstream output(filepath);
    if (!output) {
        throw std::runtime_error("Could not open project file for writing: " + filepath);
    }
    output << root.dump(2);
}

sbox::chem::MolecularSystem load_project(const std::string& filepath, nlohmann::json* out_state) {
    std::ifstream input(filepath);
    if (!input) {
        throw std::runtime_error("Could not open project file: " + filepath);
    }

    nlohmann::json root;
    try {
        input >> root;
    } catch (const std::exception& ex) {
        throw std::runtime_error(std::string("Invalid project JSON: ") + ex.what());
    }

    if (!root.contains("version")) {
        throw std::runtime_error("Invalid project JSON: missing version");
    }
    if (root.at("version").get<int>() != 1) {
        throw std::runtime_error("Unsupported project version");
    }

    sbox::chem::MolecularSystem mol;
    mol.set_name(root.value("name", ""));
    mol.set_charge(root.value("charge", 0));
    mol.set_multiplicity(root.value("multiplicity", 1));

    if (root.contains("atoms")) {
        for (const nlohmann::json& atom_json : root.at("atoms")) {
            const auto& position = atom_json.at("position");
            mol.add_atom({
                atom_json.at("Z").get<int>(),
                Eigen::Vector3d(position.at(0).get<double>(),
                                position.at(1).get<double>(),
                                position.at(2).get<double>()),
                atom_json.value("label", std::string()),
                atom_json.value("formal_charge", 0),
            });
        }
    }

    if (root.contains("bonds")) {
        for (const nlohmann::json& bond_json : root.at("bonds")) {
            mol.add_bond(bond_json.at("i").get<int>(),
                         bond_json.at("j").get<int>(),
                         static_cast<sbox::chem::BondOrder>(bond_json.at("order").get<int>()));
        }
    }

    if (out_state != nullptr) {
        if (root.contains("state") && root.at("state").is_object()) {
            *out_state = root.at("state");
        } else {
            *out_state = nlohmann::json::object();
        }
    }

    return mol;
}

}  // namespace sbox::io

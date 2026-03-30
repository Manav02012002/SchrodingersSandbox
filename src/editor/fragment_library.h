#pragma once

#include "core/molecular_system.h"

#include <Eigen/Core>

#include <string>
#include <vector>

namespace sbox::editor {

struct Fragment {
    std::string name;
    std::string category;
    sbox::chem::MolecularSystem molecule;
    int attachment_atom = -1;
};

class FragmentLibrary {
public:
    FragmentLibrary();

    const std::vector<Fragment>& all() const;
    std::vector<const Fragment*> by_category(const std::string& category) const;
    const Fragment* find(const std::string& name) const;
    std::vector<std::string> categories() const;

    sbox::chem::MolecularSystem place(const Fragment& fragment,
                                      const Eigen::Vector3d& position) const;

private:
    std::vector<Fragment> fragments_;
    void build_library();
};

}  // namespace sbox::editor

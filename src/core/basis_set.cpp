#include "core/basis_set.h"

namespace sbox::basis {

int BasisSet::num_basis_functions() const {
    int total = 0;

    for (const BasisShell& shell : shells) {
        const int l = shell.angular_momentum;
        if (l < 0) {
            continue;
        }

        if (spherical) {
            total += 2 * l + 1;
        } else {
            total += (l + 1) * (l + 2) / 2;
        }
    }

    return total;
}

}  // namespace sbox::basis

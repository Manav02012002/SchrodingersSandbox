#pragma once

#include "core/basis_set.h"

#include <Eigen/Core>

#include <vector>

namespace sbox::basis {

double evaluate_primitive(double alpha,
                          double coeff,
                          int lx,
                          int ly,
                          int lz,
                          double dx,
                          double dy,
                          double dz);

double evaluate_shell(const BasisShell& shell,
                      const Eigen::Vector3d& center,
                      const Eigen::Vector3d& point,
                      int basis_offset,
                      std::vector<double>& out_values);

void evaluate_basis_at_point(const MOData& mo_data,
                             const Eigen::Vector3d& point,
                             Eigen::VectorXd& basis_values);

double evaluate_mo_at_point(const MOData& mo_data, int mo_index, const Eigen::Vector3d& point);

double evaluate_mo_density_at_point(const MOData& mo_data, int mo_index, const Eigen::Vector3d& point);

Eigen::VectorXd evaluate_mo_on_grid(const MOData& mo_data,
                                    int mo_index,
                                    const Eigen::Vector3d& origin,
                                    const Eigen::Vector3d& step,
                                    int nx,
                                    int ny,
                                    int nz);

}  // namespace sbox::basis

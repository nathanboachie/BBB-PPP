/*  __       __          __
 * |__)||\/||__) /\ |\/||__)/  \|\/|
 * |__)||  ||__)/--\|  ||__)\__/|  |
 *
 * This file is part of BIMBAMBUM.
 *
 * -----------------------------------------------------------------------------
 * Copyright (C) 2023 Armand Sieber
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 * -----------------------------------------------------------------------------
 *
 */

/*! \file BoundaryData.hpp
    \brief Header for BoundaryData.cpp
*/

/*! \class BoundaryData
    \brief Base class containing function declaration, properties and quantities associated with the fluid-fluid
    interface.
*/

#ifndef BOUNDARYDATA_HPP
#define BOUNDARYDATA_HPP

#include <vector>
#include <string>

class BoundaryData {
public:
    template<typename Input>
    BoundaryData(Input &data) {

        Ns = data.Ns;
        delta_phi = data.delta_phi;

        gamma = data.gamma;
        alpha = data.alpha;

        r_nodes.resize(Ns + 1);
        z_nodes.resize(Ns + 1);
        ur_nodes.resize(Ns + 1);
        uz_nodes.resize(Ns + 1);
        u_nodes1.resize(Ns + 1);
        u_nodes2.resize(Ns + 1);
        u_vec_prod.resize(Ns + 1);
        F_nodes.resize(Ns + 1);
        curv_nodes.resize(Ns + 1);

        r_nodes1.resize(Ns + 1);
        z_nodes1.resize(Ns + 1);
        F_nodes1.resize(Ns + 1);
        dr1.resize(Ns + 1);
        dz1.resize(Ns + 1);
        dF1.resize(Ns + 1);
        dr2.resize(Ns + 1);
        dz2.resize(Ns + 1);
        dF2.resize(Ns + 1);


    }

    virtual ~BoundaryData() = default;


    //Input data;
    int Ns;
    double delta_phi;

    // Bubble and interface properties
    double R0;
    double gamma;
    double alpha;

    // Solution vectors
    std::vector<double> r_nodes;
    std::vector<double> z_nodes;
    std::vector<double> ur_nodes;
    std::vector<double> uz_nodes;
    std::vector<double> u_nodes1;
    std::vector<double> u_nodes2;
    std::vector<double> u_vec_prod;
    std::vector<double> F_nodes;
    std::vector<double> curv_nodes;

    // Endpoints derivatives (for cubic spline interpolation)
    double drds1;
    double drds2;
    double dzds1;
    double dzds2;

    double dphi1ds1;
    double dphi1ds2;
    double dphi2ds1;
    double dphi2ds2;

    // Time integration supplement storage
    std::vector<double> r_nodes1;
    std::vector<double> z_nodes1;
    std::vector<double> F_nodes1;
    std::vector<double> dr1;
    std::vector<double> dz1;
    std::vector<double> dF1;
    std::vector<double> dr2;
    std::vector<double> dz2;
    std::vector<double> dF2;

    virtual void initialize() = 0; // initialize fluid-fluid interface node points

    virtual void remesh_boundary() = 0; // re-grid bubble fluid-fluid interface

    virtual void filter_boundary() = 0; // filters fluid-fluid interface

    virtual void boundary_endpoints_derivatives() = 0; // estimates the derivatives at the endpoints on the fluid-fluid interface

    virtual void boundary_curvature() = 0; // computes the curvature on the fluid-fluid interface

    double time_step_boundary(); // compute time step based on fluid-fluid interface dynamics

};


#endif //BOUNDARYDATA_H

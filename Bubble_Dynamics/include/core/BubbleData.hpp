/*  __       __          __
 * |__)||\/||__) /\ |\/||__)/  \|\/|
 * |__)||  ||__)/--\|  ||__)\__/|  |
 *
 * This file is part of BIMBAMBUM.
 *
 * BubbleData.hpp
 *
 * Description:
 * Base class containing function declaration, properties and quantities
 * associated with the bubble dynamics.
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

/*! \file BubbleData.hpp
    \brief Header for BubbleData.cpp
*/

/*! \class BubbleData
    \brief Base class containing function declaration, properties and quantities associated with the bubble dynamics.
*/

#ifndef BUBBLEDATA_HPP
#define BUBBLEDATA_HPP

#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "init_R0.hpp"

using namespace std;


class BubbleData {
public:
    template<typename Input>
    BubbleData(Input &data) {

        Nb = data.Nb;
        r_nodes.resize(Nb + 1);
        z_nodes.resize(Nb + 1);
        ur_nodes.resize(Nb + 1);
        uz_nodes.resize(Nb + 1);
        u_nodes.resize(Nb + 1);
        phi_nodes.resize(Nb + 1);

        r_nodes1.resize(Nb + 1);
        z_nodes1.resize(Nb + 1);
        phi_nodes1.resize(Nb + 1);
        dr1.resize(Nb + 1);
        dz1.resize(Nb + 1);
        dphi1.resize(Nb + 1);
        dr2.resize(Nb + 1);
        dz2.resize(Nb + 1);
        dphi2.resize(Nb + 1);

        zeta = data.zeta;
        gamma = data.gamma;
        epsilon = data.epsilon;
        k = data.k;

        bubble_dynamics = data.bubble_dynamics;
        delta_phi = data.delta_phi;

    }

    virtual ~BubbleData() = default;


    //Input data;
    int Nb;
    std::vector<double> r_nodes;
    std::vector<double> z_nodes;
    std::vector<double> ur_nodes;
    std::vector<double> uz_nodes;
    std::vector<double> u_nodes;
    std::vector<double> phi_nodes;

    // Time integration supplement storage
    std::vector<double> r_nodes1;
    std::vector<double> z_nodes1;
    std::vector<double> phi_nodes1;
    std::vector<double> dr1;
    std::vector<double> dz1;
    std::vector<double> dphi1;
    std::vector<double> dr2;
    std::vector<double> dz2;
    std::vector<double> dphi2;


    // Bubble properties
    double R0;
    double V0;
    double V;
    double zeta;
    double gamma;
    double epsilon;
    double k;

    // Simulation case properties
    std::string bubble_dynamics;
    double delta_phi;
    double t0;

    virtual void initialize() = 0;

    void filter_bubble(); // filters bubble surface

    void remesh_bubble(); // re-grid bubble surface

    double compute_volume(); // compute bubble volume

    int intersect(); //check for bubble intersection

    double time_step_bubble(double epsilon, double k); // compute time step based on bubble's dynamics

};

#endif // BUBBLEDATA_HPP

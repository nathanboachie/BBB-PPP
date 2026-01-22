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

/*! \file Case_RayleighPlessetBubble.cpp
    \brief Derived class: handles the dynamics of a bubble that has been initiated with the Rayleigh-Plesset model.

    This derived class inherits members from BubbleData.cpp and BubbleData.cpp
 */

#include "Case_RayleighPlessetBubble.hpp"

#include <armadillo>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <numeric>

#include "init_R0.hpp"

using namespace std;
using namespace arma;

/**
 * nk2b 22/01/26
 * @brief Used to read initial arrays of r_nodes and z_nodes
 * @param vec Vector to be entered and changed (pass by ref)
 * @param filename full path to text file containing initital position
 **/
void readInitialArray(vector<double>& vec, string filename, const double length) {
	ifstream file(filename);
	if(!file) {
		throw runtime_error("Can't open file");
	}
	for(auto& x: vec) {
		file>>x;
	}
	if(vec.size() !=length+1) {
		throw runtime_error("Read in vector is not equal to Nb+1, Check sizes of entred arrays!");
	}
}

/**
 * nk2b 22/01/26
 * @brief Calculate an effective radius of a deformed bubble
 * @param xvec Vector containing x positions
 * @param yvec Vector containing y positions
 * @return effective radius
 **/
double effectiveRadius(vector<double>& xvec, vector<double>& yvec) {
	if(xvec.size() !=yvec.size()) {
		throw runtime_error("Vectors aren't of equal sizes");
	}
	double ycentroid=accumulate(yvec.begin(),yvec.end(),0.0)/yvec.size();
	vector<double> radiusVec(xvec.size());
	for(int i=0; i<xvec.size();++i) {
		radiusVec[i]=sqrt(xvec[i]*xvec[i]+(yvec[i]-ycentroid)*(yvec[i]-ycentroid));
	}
	double effRad=accumulate(radiusVec.begin(),radiusVec.end(),0.0)/radiusVec.size();
	return effRad;
}


/*! Initialize the position of the bubble nodes and the value of the potentials at those nodes based on the Rayleigh-Plesset model.*/
void Case_RayleighPlessetBubble::initialize() {

    const double pi = 3.14159265358979323846264338328;
		// nk2b 22/01/26
		// nk2b HARDCODE
		readInitialArray(r_nodes,"/home/exy214/Documents/cavitation/code/bbb-dir/BIMBAMBUM/SIMS/PerturbRayleigh1P0/rInitial.txt",Nb);
		readInitialArray(z_nodes,"/home/exy214/Documents/cavitation/code/bbb-dir/BIMBAMBUM/SIMS/PerturbRayleigh1P0/zInitial.txt",Nb);
    for (int i = 0; i < Nb + 1; ++i) {
        ur_nodes[i] = 0.0;
        uz_nodes[i] = 0.0;
        phi_nodes[i] = 0.0;
    }
    // impose first and last node to be exactly on the axis of symmetry
    r_nodes.front() = 0.0;
    r_nodes.back() = 0.0;
		t0 = 0.0;
    R0 = effectiveRadius(r_nodes,z_nodes); // compute initial radius based on the values of epsilon and k
    std::cout << "Bubble initial radius R0: " << "\t" << R0 << std::endl;
    V0 = 4.0 / 3.0 * pi * pow(R0, 3.0);
    V = V0;
}

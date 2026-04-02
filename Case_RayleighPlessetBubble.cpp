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

#include <armadillo>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <numeric>

using namespace std;
using namespace arma;

// Test Class
// Purpose is to initialise a bubble based on a group of initial parameters
// AKA a deformed bubble initial state

class Case_RayleighPlessetBubble {
	private:
		double R0;
		const int Nb;
		const double gamma;
		string rPositions;
		string zPositions;
		vector<double> r_nodes;
		vector<double> z_nodes;
		vector<double> ur_nodes;
		vector<double> uz_nodes;
		vector<double> phi_nodes;

		void readInitialArray(vector<double>& vec,string filename) {
        ifstream file(filename);
				if(!file) {
					throw runtime_error("Can't open file");
				}
        for(auto& x: vec) {
                file>>x;
        }
		}

		double effectiveRadius(vector<double>& xvec,vector<double>& yvec) {
			if(xvec.size() != yvec.size()) {
				throw runtime_error("Vectors aren't equal size");
			}
			double xcentroid=0;
			double ycentroid=std::accumulate(yvec.begin(),yvec.end(),0.0)/yvec.size();

			vector<double> radiusVec(xvec.size());
			for(int i=0; i<xvec.size();++i) {
				radiusVec[i]=sqrt(xvec[i]*xvec[i]+(yvec[i]-ycentroid)*(yvec[i]-ycentroid));
			}
		return accumulate(radiusVec.begin(),radiusVec.end(),0.0)/radiusVec.size();
		}	
				


			
	public:
		Case_RayleighPlessetBubble(double r0, int nb, double g,string r,string z)
			: R0(r0)
			, Nb(nb)
			, gamma(g)
			, rPositions(r)
			, zPositions(z)
			, r_nodes(nb+1)
			, z_nodes(nb+1)
			, ur_nodes(nb+1)
			, uz_nodes(nb+1)
			,	phi_nodes(nb+1)
			{}
		
		void initializeCircle() {

				const double pi = 3.14159265358979323846264338328;
				double t0 = 0.0;
				std::cout << "Bubble initial radius R0: " << "\t" << R0 << std::endl;
				double V0 = 4.0 / 3.0 * pi * pow(R0, 3.0);
				double V = V0;

				for (int i = 0; i < Nb + 1; ++i) {
						double alpha_i = i * pi / Nb;
						r_nodes[i] = R0 * sin(alpha_i);
						z_nodes[i] = R0 * cos(alpha_i) - gamma;
						ur_nodes[i] = 0.0;
						uz_nodes[i] = 0.0;
						phi_nodes[i] = 0.0;
				}

				// impose first and last node to be exactly on the axis of symmetry
				r_nodes.front() = 0.0;
				r_nodes.back() = 0.0;

		}

		void initializeArbitrary() {
			readInitialArray(r_nodes,rPositions);
			readInitialArray(z_nodes,zPositions);	
			const double pi = 3.14159265358979323846264338328;
			double t0 = 0.0;
			for (int i = 0; i < Nb + 1; ++i) {
						double alpha_i = i * pi / Nb;
						ur_nodes[i] = 0.0;
						uz_nodes[i] = 0.0;
						phi_nodes[i] = 0.0;
				}
				// impose first and last node to be exactly on the axis of symmetry
			r_nodes.front() = 0.0;
			r_nodes.back() = 0.0;
			R0=effectiveRadius(r_nodes,z_nodes);
			std::cout << "Bubble initial radius R0: " << "\t" << R0 << std::endl;
			double V0 = 4.0 / 3.0 * pi * pow(R0, 3.0);
			double V = V0;

		}


		void writeSolution(const std::vector<double> &vec, std::string filename) {
        std::ofstream file(filename);
        for(auto &x: vec) {
                file<<x<<'\n';
        }
		}
		void Run() {
			initializeArbitrary();
			writeSolution(r_nodes,"r.txt");
			writeSolution(z_nodes,"z.txt");
		}
};


int main(){
	Case_RayleighPlessetBubble Bubble(4e-5,100,1.5,"/home/exy214/Documents/cavitation/code/bbb-dir/BIMBAMBUM/SIMS/rInitial.txt",
																								 "/home/exy214/Documents/cavitation/code/bbb-dir/BIMBAMBUM/SIMS/zInitial.txt");
	Bubble.Run();	
	return 0;
}

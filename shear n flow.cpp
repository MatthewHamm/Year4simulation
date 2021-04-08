// shear n flow.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include <omp.h>
#include <iostream>
#include <boost/math/special_functions/bessel.hpp> 
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;
struct Postion {
	vector<double> x;
	vector<double> y;
};
Postion r;
string filename;

const int nw = 13;

Postion channel_new;

Postion vave;
double vysum[nw];
double vxsum[nw];


double* p;
const int no = 30;
double w = (nw + 1) * sqrt(3) / 2;
double a = 1;
double b = a * w;
double dt = 0.01;


double F_current[2 * nw * no];
double F_tick[2 * nw * no];
double d;
int Count;
double T_tot = 2000;
double F_tot;

const double pi = 3.14159265358979323846;
double F_x;
double F_y;
double q;
double x;
double y;
double v_shear = 100000000000000;
Postion F_all;
int s = static_cast <int>(T_tot / dt) - 1;

vector <double > pot_table;

double step = 0.0001;
double ostep = 1.0 / step;
double olambda = 1.0;
double rcut = 9 * a;
double Fmax = 0;
double ftab(double r);
double LinearInterp(double a, double b, double f);


// calculates force at r
// using a linear interpolation of
// a the tabulated potential
double ftab(double r)
{
	double f = 0.0;
	if (r < step) f = pot_table[0]; // special case if r < h
	else if (r > rcut) f = Fmax; // special case if r > rcut
	else
	{
		int pot_lindex = static_cast <int>(floor(r / step)) - 1;
		f = LinearInterp(pot_table[pot_lindex],
			pot_table[pot_lindex + 1],
			r * ostep - (pot_lindex + 1));

	}
	return f;

}

// performs linear interpolation

double LinearInterp(double flow, double fhigh, double d)
{
	return flow + d * (fhigh - flow);
}

void poscsv(Postion r, int no, double T_tot, double dt) {

	ofstream myfile;
	int nosteps = static_cast <int>(T_tot / dt);
	myfile.open("pos.csv");
	myfile.clear();
	for (int i = 0; i <= 2 * no * nosteps - 1; i++) {

		if (i % (2 * no) == 2 * no - 1) {
			myfile << r.x[i] << ',';
			myfile << r.y[i] << endl;
		}
		else {
			myfile << r.x[i] << ',';
			myfile << r.y[i] << ',';
		}





	}
	myfile.close();
}
double* interactions(Postion channel_pos, double v_shear, double t) {
	static double F_current[2 * nw * no];
#pragma omp parallel for num_threads(4) private(F_x,F_y,q,d,x,y,F_tot)
	for (int i = 0; i <= nw * no - 1; i++) {
		//loops through all the vortices

		F_x = 0.0;
		F_y = 0.0;


		for (int j = 0; j <= nw * no - 1; j++) {
			//loops trhough all the vortices
			if (i != j) {
				//stops vortices interacting with themselves
				y = channel_pos.y[i] - channel_pos.y[j];
				x = channel_pos.x[i] - (channel_pos.x[j]);
				if (x >= no * a / 2) {
					//imposes periodic boundary conditions 
					x = x - no * a;
				}
				else if (x <= -no * a / 2) {
					x = no * a + x;
				}


				d = pow(pow(x, 2) + pow(y, 2), 0.5);
				if (d <= 9 * a) {
					//sets a cutoff
					F_tot = ftab(d);

					F_x += (F_tot * (x / d));

					F_y += (F_tot * (y / d));
				}

			}

		}

		//force for boundary
		F_y += (-pi / a) * (-exp(-channel_pos.y[i]) + exp(-(b - channel_pos.y[i]))) / ((1 - exp(-sqrt(3) / 2)));
		for (int n = 1; n <= 2; n++) {

			q = sqrt(1 + pow(2 * pi * n / a, 2));

			F_x += (2.0 * n * pi / a) * (2 * pi / a) * (sin(2.0 * n * pi * (channel_pos.x[i] + v_shear * t) / a) * exp(-q * channel_pos.y[i]) + sin(2.0 * n * pi * channel_pos.x[i] / a) * exp(-q * (b - channel_pos.y[i]))) / (q * (1 - pow(-1, n) * exp(-q * sqrt(3) / 2)));

			F_y += (-2 * pi / a) * (cos(2.0 * n * pi * (channel_pos.x[i] + v_shear * t) / a) * -q * exp(-q * channel_pos.y[i]) + cos(2.0 * n * pi * channel_pos.x[i] / a) * q * exp(-q * (b - channel_pos.y[i]))) / (q * (1 - pow(-1, n) * exp(-q * sqrt(3) / 2)));


		}

		F_current[2 * i] = F_x;

		F_current[2 * i + 1] = F_y;

	}
	return F_current;
}



int main()
{

	for (double v = 0.1; v <= 1; v += 0.1) {
		T_tot = 200;

		v_shear = v;
		//gets number of time steps
		s = static_cast <int>(T_tot / dt) - 1;

		//clears postions of previous values
		channel_new.x.clear();
		channel_new.y.clear();

		vave.x.clear();
		vave.y.clear();
		//sets inital postion of vortices
		for (int n = 1; n <= nw; n++) {



			for (int k = 0; k <= no - 1; k++) {
				channel_new.x.push_back(a / 4 - pow(-1, n) * a / 4 + a * k);
				channel_new.y.push_back(n * b / (nw + 1));

			}

		}
		cout << channel_new.x[0] << ',' << channel_new.x[no] << endl;
		std::cout << s << endl;

		for (double t = 0.00; t <= T_tot; t = t + dt) {
			//loop through time steps


			F_all.x.clear();
			F_all.y.clear();
			//gets interaction forces
			p = interactions(channel_new, 1, cos(v_shear * t));
			for (int i = 0; i < nw * no * 2; i++) {
				F_tick[i] = *(p + i);

			}
			//set sum used for average to 0

			for (int i = 0; i <= nw - 1; i++) {
				vxsum[i] = 0;
				vysum[i] = 0;
			}

			for (int o = 0; o <= nw * no - 1; o++) {
				//loops through votrices
				F_all.x.push_back(F_tick[2 * o]);
				F_all.y.push_back(F_tick[2 * o + 1]);
				//updates postions in channel
				channel_new.x[o] = channel_new.x[o] + F_all.x[o] * dt;

				channel_new.y[o] = channel_new.y[o] + F_all.y[o] * dt;
				//boundary condtions
				if (channel_new.x[o] > a* no) {
					channel_new.x[o] = channel_new.x[o] - a * no;
				}
				if (channel_new.x[o] < 0) {
					channel_new.x[o] = channel_new.x[o] + a * no;
				}
				if (channel_new.y[o] > b) {
					channel_new.y[o] = channel_new.y[o] - b;
				}
				if (channel_new.y[o] < 0) {
					channel_new.y[o] = channel_new.y[o] + b;
				}
				//sums up the vx and vy values for vortices 
				for (int i = 0; i <= nw - 1; i++) {
					if ((i * no <= o) && (o < (i + 1) * no)) {

						vxsum[i] += abs(F_tick[2 * o]);
						vysum[i] += abs(F_tick[2 * o + 1]);
					}

				}



			}
			//stores averages
			for (int i = 0; i <= nw - 1; i++) {
				vave.x.push_back(vxsum[i] / (no));
				vave.y.push_back(vysum[i] / (no));
			}




		}



		//saves value to csv
		ofstream myfile;



		filename = "vosc" + to_string(nw) + " func dt 0.01 i t 200flow " + to_string(v) + ".csv";
		for (int i = 0; i <= nw - 1; i++) {
			vxsum[i] = 0;
			vysum[i] = 0;
		}
		Count = 0;
		std::cout << filename << endl;
		myfile.open(filename);
		myfile.clear();


		for (int n = 0; n <= s; n++) {
			for (int i = 0; i <= nw - 1; i++) {
				vxsum[i] += vave.x[nw * n + i];
				vysum[i] += vave.y[nw * n + i];

			}

		}
		for (int i = 1; i <= nw; i++) {
			myfile << i * sqrt(3) / 2 << ',' << vxsum[i] / (s * nw) << ',' << vysum[i] / (s * nw) << endl;

		}








		myfile.close();
	}

	std::cout << s << endl;
	for (long i = 1; i <= static_cast <int>(rcut / step) + 1; ++i)
	{
		double r = step * i;
		double f = boost::math::cyl_bessel_k(1, r * olambda);
		pot_table.push_back(f);
	}



	for (double v = 0.1; v <= 1; v += 0.1) {
		T_tot = 200;

		v_shear = v;

		s = static_cast <int>(T_tot / dt) - 1;


		channel_new.x.clear();
		channel_new.y.clear();

		vave.x.clear();
		vave.y.clear();
		for (int n = 1; n <= nw; n++) {



			for (int k = 0; k <= no - 1; k++) {
				channel_new.x.push_back(a / 4 - pow(-1, n) * a / 4 + a * k);
				channel_new.y.push_back(n * b / (nw + 1));

			}

		}
		cout << channel_new.x[0] << ',' << channel_new.x[no] << endl;
		std::cout << s << endl;

		for (double t = 0.00; t <= T_tot; t = t + dt) {



			F_all.x.clear();
			F_all.y.clear();
			p = interactions(channel_new, v_shear, t);
			for (int i = 0; i < nw * no * 2; i++) {
				F_tick[i] = *(p + i);

			}

			for (int i = 0; i <= nw - 1; i++) {
				vxsum[i] = 0;
				vysum[i] = 0;
			}

			for (int o = 0; o <= nw * no - 1; o++) {
				F_all.x.push_back(F_tick[2 * o]);
				F_all.y.push_back(F_tick[2 * o + 1]);

				channel_new.x[o] = channel_new.x[o] + F_all.x[o] * dt;

				channel_new.y[o] = channel_new.y[o] + F_all.y[o] * dt;

				if (channel_new.x[o] > a* no) {
					channel_new.x[o] = channel_new.x[o] - a * no;
				}
				if (channel_new.x[o] < 0) {
					channel_new.x[o] = channel_new.x[o] + a * no;
				}
				if (channel_new.y[o] > b) {
					channel_new.y[o] = channel_new.y[o] - b;
				}
				if (channel_new.y[o] < 0) {
					channel_new.y[o] = channel_new.y[o] + b;
				}
				for (int i = 0; i <= nw - 1; i++) {
					if ((i * no <= o)  && (o< (i + 1) * no)) {

						vxsum[i] += abs(F_tick[2 * o]);
						vysum[i] += abs(F_tick[2 * o + 1]);
					}

				}



			}
			for (int i = 0; i <= nw - 1; i++) {
				vave.x.push_back(vxsum[i] / ( no));
				vave.y.push_back(vysum[i] / (no));
			}




		}




		ofstream myfile;



		filename = "vlin"+to_string(nw)+" func dt 0.01 i t 200flow " + to_string(v) + ".csv";
		for (int i = 0; i <= nw - 1; i++) {
			vxsum[i] = 0;
			vysum[i] = 0;
		}
		Count = 0;
		std::cout << filename << endl;
		myfile.open(filename);
		myfile.clear();


		for (int n = 0; n <= s; n++) {
			for (int i = 0; i <= nw - 1; i++) {
				vxsum[i] += vave.x[nw * n + i];
				vysum[i] += vave.y[nw * n + i];

			}
			
		}
		for (int i = 1; i <= nw ; i++) {
			myfile <<i*sqrt(3)/2 <<','<<vxsum[i] / (s * nw) << ',' << vysum[i] / (s * nw) << endl;

		}
	







		myfile.close();
	}

}

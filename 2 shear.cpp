// 2 shear.cpp : This file contains the 'main' function. Program execution begins and ends there.
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

string filename;

const int nw = 2;

Postion channel_new;


double* p;
const int no = 30;
double w = (nw + 1) * sqrt(3) / 2;
double a = 1;
double w = a * w;
double dt = 0.01;
double xsum[nw];
Postion rave;
double ysum[nw];

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
double LinearInterp(double a, double w, double f);


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


double* interactions(Postion channel_pos, double v_shear, double t) {

	static double F_current[2 * nw * no];
#pragma omp parallel for num_threads(4) private(F_x,F_y,q,d,x,y,F_tot)
	for (int i = 0; i <= nw * no - 1; i++) {
		//loops through vortices

		F_x = 0.0;
		F_y = 0.0;


		for (int j = 0; j <= nw * no - 1; j++) {
			//loops through vortices
			if (i != j) {
				//ensure vortices don't interact with self
				y = channel_pos.y[i] - channel_pos.y[j];
				x = channel_pos.x[i] - (channel_pos.x[j]);
				if (x >= no * a / 2) {
					//imposes periodic boundary conditions
					x = x - no * a;
				}
				else if (x <= -no * a / 2) {
					x = no * a + x;
				}

				//finds distance bewteen vortices
				d = pow(pow(x, 2) + pow(y, 2), 0.5);
				if (d <= 9 * a) {
					F_tot = ftab(d);

					F_x += (F_tot * (x / d));

					F_y += (F_tot * (y / d));
				}

			}

		}

		//force for boundary
		F_y += (-pi / a) * (-exp(-channel_pos.y[i]) + exp(-(w - channel_pos.y[i]))) / ((1 - exp(-sqrt(3) / 2)));
		for (int n = 1; n <= 2; n++) {




			q = sqrt(1 + pow(2 * pi * n / a, 2));

			F_x += (2.0 * n * pi / a) * (2 * pi / a) * (sin(2.0 * n * pi * (channel_pos.x[i] + v_shear * t) / a) * exp(-q * channel_pos.y[i]) + pow(-1, n)* sin(2.0 * n * pi * channel_pos.x[i] / a) * exp(-q * (w - channel_pos.y[i]))) / (q * (1 - pow(-1, n) * exp(-q * sqrt(3) / 2)));




			F_y += (-2 * pi / a) * (cos(2.0 * n * pi * (channel_pos.x[i] + v_shear * t) / a) * -q * exp(-q * channel_pos.y[i]) + pow(-1, n)*cos(2.0 * n * pi * channel_pos.x[i] / a) * q * exp(-q * (w - channel_pos.y[i]))) / (q * (1 - pow(-1, n) * exp(-q * sqrt(3) / 2)));






		}







		F_current[2 * i] = F_x;

		F_current[2 * i + 1] = F_y;



	}
	return F_current;
}



int main()
{



	//creates table of bessel function values
	for (long i = 1; i <= static_cast <int>(rcut / step) + 1; ++i)
	{
		double r = step * i;
		double f = boost::math::cyl_bessel_k(1, r * olambda);
		pot_table.push_back(f);
	}



	for (double v = 1; v <= 1; v += 0.0001) {
		//selects T_tot such that an integer numbe rof lattice sapcing is moved
		T_tot = 100;
		for (int i = 10; i >= 2; --i) {
			if (i / v >= 100) {
				T_tot = i / v;

			}
		}
		v_shear = v;
		//gets the number of lattice spacings
		s = static_cast <int>(T_tot / dt) - 1;

		//Clears postion struct of previous values
		channel_new.x.clear();
		channel_new.y.clear();
		rave.x.clear();
		rave.y.clear();

		//places vortices in intial postions
		for (int n = 1; n <= nw; n++) {



			for (int k = 0; k <= no - 1; k++) {
				channel_new.x.push_back(a / 4 - pow(-1, n) * a / 4 + a * k);
				channel_new.y.push_back(n * w / (nw + 1));

			}

		}

		std::cout << s << endl;

		for (double t = 0.00; t <= T_tot; t = t + dt) {
			//loops through time steps


			F_all.x.clear();
			F_all.y.clear();
			p = interactions(channel_new, v_shear, t);
			for (int i = 0; i < nw * no * 2; i++) {
				//gets value form interaction func
				F_tick[i] = *(p + i);
				
			}
			//set sum used for average to 0
			xsum[0] = 0;
			ysum[0] = 0;
			xsum[1] = 0;
			ysum[1] = 0;

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
				if (channel_new.y[o] > w) {
					channel_new.y[o] = channel_new.y[o] - w;
				}
				if (channel_new.y[o] < 0) {
					channel_new.y[o] = channel_new.y[o] + w;
				}
				//sums up the x and y values for vortices 
				if (o>=no) {

					xsum[1] += channel_new.x[o] - (a / 4 -pow(-1, (o / no) + 1) * a / 4 + a * (o % no));
					ysum[1] += channel_new.y[o] - (o / no + 1) * w / (nw + 1);

				}
				else {
			

					xsum[0] += channel_new.x[o] - (a / 4 - pow(-1, (o / no) + 1) * a / 4 + a * (o % no));
					ysum[0] += channel_new.y[o] - (o / no + 1) * w / (nw + 1);
				}


			}
			//stores the averages
			rave.x.push_back(xsum[0] / (no));
			rave.y.push_back(ysum[0] / (no));
			rave.x.push_back(xsum[1] / (no));
			rave.y.push_back(ysum[1] / (no));

		}



		//saves values to csv
		ofstream myfile;



		filename = "vlin new func dt 0.01 i t 200pos " + to_string(v) + ".csv";
		xsum[0] = 0;
		ysum[0] = 0;
		xsum[1] = 0;
		ysum[1] = 0;
		Count = 0;
		std::cout << filename << endl;
		myfile.open(filename);
		myfile.clear();


		for (int n = 0; n <= s; n++) {
			if (0.5 < rave.x[2 * n + 1]) {
				rave.x[2 * n + 1] = rave.x[2 * n + 1] - 1;
			}
			xsum[0] += rave.x[2*n];
			ysum[0] += rave.y[2*n];
			xsum[1] += rave.x[2 * n+1];
			ysum[1] += rave.y[2 * n+1];
			Count += 1;
			myfile << rave.x[2*n] << ',' << rave.y[2*n] <<','<< rave.x[2 * n+1] << ',' << rave.y[2 * n+1]<<endl;
		}
		myfile << '0' << ',' << '0' << ',' << '0' << ',' << '0' << ',' << '0' << ',' << '0' << ',' << '0' << ',' << '0' << endl;

		myfile << xsum[0] / Count << ',' << ysum[0] / Count <<','<< xsum[1] / Count << ',' << ysum[1] / Count << endl;








		myfile.close();
	}
	//gets equilibrum postions
	for (double v = 0.5; v <= 0.5; v += 0.000001) {
		T_tot = 200;

		v_shear = v;

		s = static_cast <int>(T_tot / dt) - 1;


		channel_new.x.clear();
		channel_new.y.clear();
		rave.x.clear();
		rave.y.clear();

		for (int n = 1; n <= nw; n++) {



			for (int k = 0; k <= no - 1; k++) {
				channel_new.x.push_back(a / 4 - pow(-1, n) * a / 4 + a * k);
				channel_new.y.push_back(n * w / (nw + 1));

			}

		}
		std::cout << channel_new.x[0] << ',' << channel_new.x[no] << endl;
		std::cout << s << endl;

		for (double t = 0.00; t <= T_tot; t = t + dt) {



			F_all.x.clear();
			F_all.y.clear();
			p = interactions(channel_new, v_shear, 1);
			for (int i = 0; i < nw * no * 2; i++) {
				F_tick[i] = *(p + i);

			}
			xsum[0] = 0;
			ysum[0] = 0;
			xsum[1] = 0;
			ysum[1] = 0;

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
				if (channel_new.y[o] > w) {
					channel_new.y[o] = channel_new.y[o] - w;
				}
				if (channel_new.y[o] < 0) {
					channel_new.y[o] = channel_new.y[o] + w;
				}
				if (o >= no) {

					xsum[1] += channel_new.x[o] - (a / 4 - pow(-1, (o / no) + 1) * a / 4 + a * (o % no));
					ysum[1] += channel_new.y[o] - (o / no + 1) * w / (nw + 1);

				}
				else {


					xsum[0] += channel_new.x[o] - (a / 4 - pow(-1, (o / no) + 1) * a / 4 + a * (o % no));
					ysum[0] += channel_new.y[o] - (o / no + 1) * w / (nw + 1);
				}


			}

			rave.x.push_back(xsum[0] / (no));
			rave.y.push_back(ysum[0] / (no));
			rave.x.push_back(xsum[1] / (no));
			rave.y.push_back(ysum[1] / (no));


		}




		ofstream myfile;



		filename = "equ " + to_string(v) + ".csv";
		xsum[0] = 0;
		ysum[0] = 0;
		xsum[1] = 0;
		ysum[1] = 0;
		Count = 0;
		std::cout << filename << endl;
		myfile.open(filename);
		myfile.clear();


		for (int n = 0; n <= s; n++) {
			if (0.5 < rave.x[2 * n + 1]) {
				rave.x[2 * n + 1] = rave.x[2 * n + 1] - 1;
			}
			xsum[0] += rave.x[2 * n];
			ysum[0] += rave.y[2 * n];
			xsum[1] += rave.x[2 * n + 1];
			ysum[1] += rave.y[2 * n + 1];
			Count += 1;
			myfile << rave.x[2 * n] << ',' << rave.y[2 * n] << ',' << rave.x[2 * n + 1] << ',' << rave.y[2 * n + 1] << endl;
		}








		myfile.close();
	}


	//The equilibrum position is calculated for the time avergaed wall
	channel_new.x.clear();
	channel_new.y.clear();
	rave.x.clear();
	rave.y.clear();


	for (int n = 1; n <= nw; n++) {



		for (int k = 0; k <= no - 1; k++) {
			channel_new.x.push_back(a / 4 - pow(-1, n) * a / 4 + a * k);
			channel_new.y.push_back(n * w / (nw + 1));

		}

	}
	std::cout << channel_new.x[0] << ',' << channel_new.x[no] << endl;
	std::cout << s << endl;

	std::cout << "parallel" << endl;
	for (double t = 0.00; t <= T_tot; t = t + dt) {



		F_all.x.clear();
		F_all.y.clear();
#pragma omp parallel for num_threads(4) private(F_x,F_y,q,d,x,y,F_tot)
		for (int i = 0; i <= nw * no - 1; i++) {

			F_x = 0.0;
			F_y = 0.0;


			for (int j = 0; j <= nw * no - 1; j++) {

				if (i != j) {
					y = channel_new.y[i] - channel_new.y[j];
					x = channel_new.x[i] - (channel_new.x[j]);
					if (x >= no * a / 2) {
						x = x - no * a;
					}
					else if (x <= -no * a / 2) {
						x = no * a + x;
					}


					d = pow(pow(x, 2) + pow(y, 2), 0.5);
					if (d <= 9 * a) {
						F_tot = ftab(d);

						F_x += (F_tot * (x / d));

						F_y += (F_tot * (y / d));
					}

				}

			}

			//force for boundary
			F_y += (-pi / a) * (-exp(-channel_new.y[i]) + exp(-(w - channel_new.y[i]))) / ((1 - exp(-sqrt(3) / 2)));
			for (int n = 1; n <= 1; n++) {




				q = sqrt(1 + pow(2 * pi * n / a, 2));

				F_x += (2.0 * n * pi / a) * (2 * pi / a) * pow(-1, n) * (sin(2.0 * n * pi * channel_new.x[i] / a) * exp(-q * (w - channel_new.y[i]))) / (q * (1 - pow(-1, n) * exp(-q * sqrt(3) / 2)));




				F_y += (-2 * pi / a) * pow(-1, n) * (cos(2.0 * n * pi * channel_new.x[i] / a) * q * exp(-q * (w - channel_new.y[i]))) / (q * (1 - pow(-1, n) * exp(-q * sqrt(3) / 2)));






			}







			F_current[2 * i] = F_x;

			F_current[2 * i + 1] = F_y;

			F_x = 0.0;
			F_y = 0.0;

		}
		xsum[0] = 0;
		ysum[0] = 0;
		xsum[1] = 0;
		ysum[1] = 0;

		for (int o = 0; o <= nw * no - 1; o++) {
			F_all.x.push_back(F_current[2 * o]);
			F_all.y.push_back(F_current[2 * o + 1]);

			channel_new.x[o] = channel_new.x[o] + F_all.x[o] * dt;

			channel_new.y[o] = channel_new.y[o] + F_all.y[o] * dt;

			if (channel_new.x[o] > a* no) {
				channel_new.x[o] = channel_new.x[o] - a * no;
			}
			if (channel_new.x[o] < 0) {
				channel_new.x[o] = channel_new.x[o] + a * no;
			}
			if (channel_new.y[o] > w) {
				channel_new.y[o] = channel_new.y[o] - w;
			}
			if (channel_new.y[o] < 0) {
				channel_new.y[o] = channel_new.y[o] + w;
			}
			if (o >= no) {

				xsum[1] += channel_new.x[o] - (a / 4 - pow(-1, (o / no) + 1) * a / 4 + a * (o % no));
				ysum[1] += channel_new.y[o] - (o / no + 1) * w / (nw + 1);

			}
			else {


				xsum[0] += channel_new.x[o] - (a / 4 - pow(-1, (o / no) + 1) * a / 4 + a * (o % no));
				ysum[0] += channel_new.y[o] - (o / no + 1) * w / (nw + 1);
			}


		}

		rave.x.push_back(xsum[0] / (no));
		rave.y.push_back(ysum[0] / (no));
		rave.x.push_back(xsum[1] / (no));
		rave.y.push_back(ysum[1] / (no));



	}




	ofstream myfile;



	filename = "time average .csv";
	xsum[0] = 0;
	ysum[0] = 0;
	xsum[1] = 0;
	ysum[1] = 0;
	Count = 0;
	std::cout << filename << endl;
	myfile.open(filename);
	myfile.clear();


	for (int n = 0; n <= s; n++) {
		if (0.5 < rave.x[2 * n + 1]) {
			rave.x[2 * n + 1] = rave.x[2 * n + 1] - 1;
		}
		xsum[0] += rave.x[2 * n];
		ysum[0] += rave.y[2 * n];
		xsum[1] += rave.x[2 * n + 1];
		ysum[1] += rave.y[2 * n + 1];
		Count += 1;
		myfile << rave.x[2 * n] << ',' << rave.y[2 * n] << ',' << rave.x[2 * n + 1] << ',' << rave.y[2 * n + 1] << endl;
	}








	myfile.close();
}


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

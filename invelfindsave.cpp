#include <iostream>
#include <cmath>
#include<fstream>

using namespace std;

int gravity(double q[], double t, double par[], double deriv[]){
	double gm = par[0];

	double r = sqrt(q[0]*q[0] + q[1] *q[1]);

	deriv[0] = q[2];
	deriv[1] = q[3];

	
	deriv[2] = -gm * q[0] /r/r/r;
	deriv[3] = -gm * q[1] /r/r/r;

	return 0;
}

int rk4( int n, double q[], double t, double tau, 
int(*derivative)(double q[], double t, double par[], double deriv[]),
double par[] ){

	int c1, c2, c3, c4;
	double  *k1, *k2, *k3, *k4;
	double *qtmp;
	double t2;
	k1 = new double[n];
	k2 = new double[n];
	k3 = new double[n];
	k4 = new double[n];
	qtmp = new double[n];

	c1 = derivative( q, t, par, k1);
	for(int i = 0; i<n; i++){
		qtmp[i] = q[i] + k1[i]*tau/2.0;
	}
	t2 = t+tau/2.0;

	c2 = derivative(qtmp, t2, par, k2);
	for(int i = 0; i<n; i++){
		qtmp[i] = q[i] + k2[i]*tau/2.0;
	}

	c3 = derivative(qtmp, t2, par, k3);
	for(int i = 0; i<n; i++){
		qtmp[i] = q[i] + k3[i]*tau;
	}
	t2=t + tau;
	c4 = derivative(qtmp, t2, par, k4);
	

	for(int i = 0; i<n; i++){
		q[i] += (k1[i]/6. +k2[i]/3. +k3[i]/3. + k4[i]/6.)*tau;
	}

	delete[] qtmp;
	delete[] k4;
	delete[] k3;
	delete[] k2;
	delete[] k1;

	return c1 + c2 + c3 + c4;
}


int main(int argc, char *argv[]){

	double theta;
	double x0,y0,vx0,vy0, v;
	double t, tau, T, target_T;
	double A, P, a, e;
	double *par;
	double PI = acos(-1);
	par = new double[1];
	par[0] = 4*PI*PI;
	theta = 0;

	double period = 0.0; 
	bool period_found = false; 
	T = 1000;

	double *q, *q0;
	q = new double[4];
	q0 = new double[4];

	if(argc < 6){
		cout << "Usage: " << argv[0] << " A P e target_T tau" << endl;
		cout << "A and P are the Aphilion and Perihilion of the planet in 10^6 Km, e is the eccentricity of the orbit, and target_T is the period" << endl;
		
		return -1;
	}

	A = atof(argv[1]);
	P = atof(argv[2]);
	e = atof(argv[3]);
	target_T = atof(argv[4]);
	tau = atof(argv[5]);
	
	
	A = A * pow(10,6) /(1.496 * pow(10,8));
	P = P * pow(10,6) /(1.496 * pow(10,8));
	a = 0.5*(A+P);
	cout << A << endl;
	cout << P << endl;
	q0[0] = A;
	q0[1] = 0.;
	v = sqrt(par[0]*(1./a)*((1.+e)/(1.-e)));
	q0[2] = v * sin(theta);
	q0[3] = v * cos(theta);
	
	T = 200; //none of the periods are supposed to be greater than 200 years

	for(int i =0; i<4; i++){

		q[i] = q0[i];
	}

t = 0.0; 
bool period_found1 = false, period_found2 = false;
double int1 = 0.0, int2 = 0.0;
double last_y = q[1]; 
// bool yc = true;
double total_period;
bool dope = true;
double delta_theta = 0.001;
ofstream fp;
fp.open("invelfind.dat");
// while(dope && theta < 2*PI){
// 	t=0;
// 	int1 = int2 = 0;
// 	q[0] = A;
// 	q[1] = 0.;
// 	v = sqrt(par[0]*(1./a)*((1.+e)/(1.-e)));
// 	q[2] = v * sin(theta);
// 	q[3] = v * cos(theta);
// 	// cout << q[2] << "\t" <<q[3]<< endl;
// 	// cout << theta << endl;

while (t < T && (!period_found1 || !period_found2)) {
    rk4(4, q, t, tau, gravity, par);
    t += tau;

    if ((period_found1 && !period_found2) && (last_y > 0 && q[1] <= 0)) {
            int2 = t ;
            period_found2 = true;
    }

    if (!period_found1 && (last_y > 0 && q[1] <= 0)) {
        int1 = t;
        period_found1 = true;
    }

    last_y = q[1];

    double r = sqrt(q[0]*q[0] + q[1] *q[1]);
    // cout << t << endl;
     fp << t << "\t" << q[0] << "\t" << q[1] << "\t" << q[2] << "\t" << q[3] << endl;
}
total_period = int2 - int1;
cout << total_period<< endl;
// cout << "THIS IS THE DIFFERNECE:  " << abs(total_period-target_T) << endl;
// cout << "period: " << total_period << endl;


if(abs(total_period-target_T) < 0.1){
	cout << "THIS IS THE DIFFERNECE:  " << abs(total_period-target_T) << endl;
	cout << "period: " << total_period << endl;
	dope = false;
	q[2] = v * cos(theta);
	q[3] = v * sin(theta);
	cout << "This is the theta: " << theta << endl;
	cout << "these are the vx: "<< q[2] << "these are the vy: " << q[3] << endl;
}
else{
	theta+=delta_theta;
}	


// cout << int2 << endl;

// if (period_found1 && period_found2) {
//     total_period = int2 - int1;
//     cout << "The orbital period is: " << total_period << " years"<< endl;
// }
// else{
// 	cout <<" Either the orbit is non-stable or the total time allowed (T) did not permit the completion of 1 full period. Try increasing the value of T" <<endl;
// }

// if(abs(total_period-target_T) < 0.01){
// 	cout << "YOOOOOOO" << endl;
// 	dope = false;
// }
// else{
// 	theta+=delta_theta;
// }

// }

fp.close();
	     /* Now to plot */
  FILE *gnuplot = popen("gnuplot", "w");
  fprintf(gnuplot, "set out 'invelfind.ps'\n");
  fprintf(gnuplot, "set term post land \n");

  // fprintf(gnuplot, "set size square\n");
  // fprintf(gnuplot, "set xlabel 'x'\n");
  // fprintf(gnuplot, "set ylabel 'y'\n");


  // First plot
fprintf(gnuplot, "set xlabel 'X Position'\n");  // Label x-axis
fprintf(gnuplot, "set ylabel 'Y Position'\n");  // Label y-axis
fprintf(gnuplot, "set title 'Position'\n");  // Set title for first plot
fprintf(gnuplot, "plot 'invelfind.dat' u 2:3 w l title 'Position'\n");
fprintf(gnuplot, "plot 'invelfind.dat' u 1:2 w l title 'Position'\n");
fprintf(gnuplot, "plot 'invelfind.dat' u 1:3 w l title 'Position'\n");

  fflush(gnuplot);// Just to ensure writing out 
  pclose(gnuplot);
  // // The next command converts the gnuplot output to PDF
  system((char *)"ps2pdf invelfind.ps");

  delete[] q;
  delete[] q0;
  delete[] par;


	return 0;
}
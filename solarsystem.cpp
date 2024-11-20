#include <iostream>
#include <cmath>
#include<fstream>

using namespace std;


const int planet = 5; //number of planets
// double q[planet*4]; // array for planet position and velocity
// double mass[planet]; // masses for all the planets

int print(double print){

	cout << print << endl;

	return 0;
}

int gravity(int pn,double q[], double t, double par[], double deriv[]){
	double gm = par[0];
	double ax = 0;
	double ay = 0;


// this for loop adds to all of the forces from each planet 
// we can do this because that would be the sum of forces

for (int j = 1; j < planet; ++j){
	ax = 0;
	ay = 0;
	for (int i = 0; i < planet; i++) {
			if (i == j) continue; 
      double x = (q[0 + 4 * i] - q[0 + 4*j]); 
      double y = (q[1 + 4 * i] - q[1 + 4*j]);
      // double x = - q[0 + 4*j]; 
      // double y = - q[1 + 4*j];
     	
     	// double x = -q[0+4*j];
     	// print(x);
     	// print(q[4]);
     	// double y = -q[1 + 4*j];
			// print(y);
     	// print(q[5]);
      double r = sqrt(x*x + y*y);
      
      double f = gm * par[i+1] / (r*r*r);
      // print(par[i]);
      // double f = gm / (r*r*r);
      // cout << par[i] << endl;
      // if(r==0)
      // {
      // 	// print(4);
      // 	print(r);
    	// 	print(t);
  		// }
  		// print(x);
  		// print(y);
      // print(f);
      // print(r);
  		// cout << ax << endl;

      ax += f * x;
      ay += f * y;
      // print (ax);
      // print (ay);

      // cout << ay << endl;
  	

  }
  	deriv[j*4 + 0] = q[j*4 + 2];
    deriv[j*4 + 1] = q[j*4 + 3];
    deriv[j*4 + 2] = ax;
    deriv[j*4 + 3] = ay;
}
		

	return 0;
}

int rk4( int pn, int n, double q[], double t, double tau, 
int(*derivative)(int pn, double q[], double t, double par[], double deriv[]),
double par[] ){

	int c1, c2, c3, c4;
	// double  *k1, *k2, *k3, *k4;
	// double *qtmp;
	double t2;
	double k1[planet*4];
	double k2[planet*4];
	double k3[planet*4];
	double k4[planet*4];
	// k2 = new double[planet][4];
	// k3 = new double[planet][4];
	// k4 = new double[planet][4];
	double qtmp[planet*4];

	c1 = derivative( pn, q, t, par, k1);
	for(int i = 4; i<(planet*4); i++){
		qtmp[i] = q[i] + k1[i]*tau/2.0;
	}
	// print(n);
	t2 = t+tau/2.0;
	

	c2 = derivative(pn, qtmp, t2, par, k2);
	for(int i = 4; i<planet*4; i++){
		qtmp[i] = q[i] + k2[i]*tau/2.0;
	}

	c3 = derivative(pn, qtmp, t2, par, k3);
	for(int i = 4; i<planet*4; i++){
		qtmp[i] = q[i] + k3[i]*tau;
	}
	t2=t + tau;
	c4 = derivative(pn, qtmp, t2, par, k4);
	
	// for (int i = 0; i < planet *4; ++i)
	// 	{
	// 		print(k4[i]);

	// 	}	

	for(int i = 4; i<planet*4; i++){
		q[i] += (k1[i]/6. +k2[i]/3. +k3[i]/3. + k4[i]/6.)*tau;
	}

	// delete[] qtmp;
	// delete[] k4;
	// delete[] k3;
	// delete[] k2;
	// delete[] k1;

	return c1 + c2 + c3 + c4;
}

int main(int argc, char *argv[]){ ////////////////////////////////////

    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <T_value>" << endl;
        return 1;  
    }

 
   

//initailizations //initailizations //initailizations //initailizations 
	double q[planet*4];
	double x0,y0,vx0,vy0;
	double t, tau, T;
	double ke, pe, te, lz;
	double save;
	t = 0;
	double *par;
	double PI = acos(-1);
	par = new double[1 + planet];
	par[0] = 4*PI*PI;
	tau = 0.001;
 	T = std::stod(argv[1]);

	ifstream mass_dat("mass.dat"); //masses
	ifstream ssinitial("ssinitial.dat"); //solar system parameters

	
	
	for(int i = 1; i < (1 + planet);i++){
		mass_dat >> par[i]; //mass.dat has masses in terms of earth masses
		par[i] = par[i]/333030.; // converts them to sun masses
		// print(par[i]);
	}
	for (int i = 0; i < (planet+1); ++i)
	{
		print(par[i]);
	}

	for(int i = 0; i<planet; i++){
		for(int j =0; j < 4; j++){
			ssinitial >> save;
			q[j + 4 * i] = save;
			print(save);
		}
	}


	ofstream fp;
	fp.open("solarsystem.dat");
	fp << t << "\t";
	for(int i = 1; i<planet; i++){
			for(int j =0; j < 4; j++){
				fp << q[j + 4 * i]<< "\t";
			}
		}
	fp << endl;

	while(t < T){				//WHILE//WHILE//WHILE//WHILE//WHILE//WHILE//WHILE//
		fp << t <<"\t";
		
			rk4( 0, 4, q,  t, tau, gravity, par);
				t += tau;
		for (int i = 1; i < planet; i++){	
			for(int j =0; j < 4; j++){
				fp << q[j + 4 * i]<< "\t";
			}
		}
		fp << endl;
	}
	fp.close();
  delete[] par;

// FILE *gnuplot = popen("gnuplot", "w");
// fprintf(gnuplot, "set output 'solarsystem.ps'\n");
// fprintf(gnuplot, "set terminal postscript landscape enhanced color\n");
// fprintf(gnuplot, "plot 'solarsystem.dat' using 2:3 with lines, \n");
// fprintf(gnuplot, "plot 'solarsystem.dat' using 1:2 with lines, \n");
// fprintf(gnuplot, "plot 'solarsystem.dat' using 1:3 with lines, \n");
// fflush(gnuplot); 
// pclose(gnuplot);
// system("ps2pdf solarsystem.ps");

FILE *gnuplot = popen("gnuplot", "w");
fprintf(gnuplot, "set output 'solarsystem.ps'\n");
fprintf(gnuplot, "set terminal postscript landscape enhanced color\n");
fprintf(gnuplot, "set size square\n");
fprintf(gnuplot, "set xlabel 'x'\n");
fprintf(gnuplot, "set ylabel 'y'\n");
// Start the plot command
fprintf(gnuplot, "plot ");

for (int i = 0; i < planet-1; ++i) {    
    fprintf(gnuplot, "'solarsystem.dat' using %d:%d with lines title 'Planet %d'",2+i*4,3+i*4 , i);
    if (i < planet - 1) {
        fprintf(gnuplot, ", ");
    }
}

// Close the plot command
fprintf(gnuplot, "\n");

fflush(gnuplot);
pclose(gnuplot);
system("ps2pdf solarsystem.ps");


	return 0;
}
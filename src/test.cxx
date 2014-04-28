/*
 * test.cxx
 * This file is part of optim
 *
 * Copyright (C) 2012 - Christian Diener
 *
 * optim is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * optim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with optim. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>

#include "lsoda.h"
#include "system.h"
#include "optim.h"

#define DIM 8

using namespace std;

// Global variables needed for fit

//TODO: there must be a solver instance here
std::vector<double> fit_time;		// time points used for fitting (where t_0 = 0)
std::vector<double> init_vars;		// initial variable values

void fitness(double* x, double* f, double* g)
{
	(*f) = 0.0;
	for(unsigned int i=0; i<DIM; i++) (*f) += x[i]*sin( sqrt( fabs(x[i]) ) );

	(*g) = 0.0;
}



int main (int argc, char const* argv[])
{		
	double t[10];
	for(unsigned int i=0; i<10; i++) t[i] = i/0.4;

	auto start = chrono::high_resolution_clock::now();
	Lsoda solver(sys, jac, p, ystart, neq, np, 0);
	solver.timecourse_to_file(t, 10, "test.txt");
	auto end = chrono::high_resolution_clock::now();

	double* y = solver.get_y();
	cout<<"y("<<solver.get_t()<<") = "<<y[0]<<"\t"<<y[1]<<endl;
	cout<<"Needed "<<chrono::duration<double, milli>(end-start).count()<<" ms."<<endl;
	cout<<solver;

	
	int myid;
	double mt1, mt2;
	
	MPI_SRES* opt = new MPI_SRES(fitness, DIM, 0);
	
	opt->toFile = 0;
	opt->logFreq = 100;
	for(unsigned int i=0; i<DIM; i++) 
	{
		opt->lb[i] = -500.0;
		opt->ub[i] = 500.0;
	}
	opt->init();
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	
	mt1 = MPI_Wtime();	
	*opt+=1000;	

	mt2 = MPI_Wtime();
	
	if(opt->converged() && myid==0) std::cout<<"Converged!"<<std::endl;	
	
	delete opt;
	if(myid == 0) std::cout<<"Needed "<<mt2-mt1<<" s."<<std::endl;	

	return 0;
}

 

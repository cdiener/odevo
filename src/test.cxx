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
#include <vector>
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

int main (int argc, char* argv[])
{	
	MPI_Init(&argc, &argv); 
	int pID;
	double start, stop;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &pID);
		
	vector<double> t;
	for(unsigned int i=0; i<10; i++) t.push_back( i/0.4 );

	if(pID == 0)
	{
		start = MPI_Wtime();
		lsoda solver(sys, jac, p, ystart, neq, np, 0);
		solver.timecourse_to_file(t, "test.txt");
		stop = MPI_Wtime();

		double* y = solver.get_y();
		cout<<"y("<<solver.get_t()<<") = "<<y[0]<<"\t"<<y[1]<<endl;
		cout<<"Needed "<<(stop-start)*1.0e3<<" ms."<<endl;
		cout<<solver;
	}
	
	mpi_sres opt(fitness, DIM, 0);
	
	opt.toFile = 0;
	opt.logFreq = 100;
	
	for(unsigned int i=0; i<DIM; i++) 
	{
		opt.lb[i] = -500.0;
		opt.ub[i] = 500.0;
	}
	opt.init();
	
	start = MPI_Wtime();	
	opt+=1000;	

	stop = MPI_Wtime();
	
	if(opt.converged() && pID==0) std::cout<<"Converged!"<<std::endl;	
	
	if(pID == 0) std::cout<<"Needed "<<stop-start<<" s."<<std::endl;	

	return 0;
}

 

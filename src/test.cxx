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
#include <random>
#include <algorithm>

#include "lsoda.h"
#include "system.h"
#include "optim.h"

// Parameters for the fitting test
const double noise_sd = 0.1;
const double tmax = 10.0;

using namespace std;

// Global variables needed for fit

unsigned int nt = 16;				// number of time points 
vector<double> fit_time;			// time points used for fitting (where t_0 = 0)
//vector<double> init_vars;			// initial variable values
double *pars;						// the parameters
vector<int> fit_idx;				// indices of the variables to be fitted
vector<vector<double> > data_means;	// means of the measured data points
vector<vector<double> > data_vars;	// variances of the measured data points
double data_mean_var = 1.0;			// mean variance of the data, used if no variance available

// the ode solver
lsoda solver(sys, jac, pars, ystart, neq, np, use_jac);	

// the fitness function we use
void fitness(double* x, double* f, double* g)
{
	for(unsigned int i=0; i<np; i++) p[i] = x[i];
	vector<vector<double> > tc = solver.timecourse(fit_time);
	double var = 0.0;
	
	(*f) = 0.0;
	for(unsigned int i=0; i<fit_time.size(); i++)
		for(unsigned int j=0; j<fit_idx.size(); j++) 
		{
			var = data_vars[i][j];
			if( var < 1.0e-6 ) var = data_mean_var;
			(*f) += pow(tc[i][ fit_idx[j] ] - data_means[i][j], 2)/var;
		}

	(*g) = 0.0;
}

int main (int argc, char* argv[])
{	
	MPI_Init(&argc, &argv); 
	int pID;
	double start, stop;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &pID);
	
	double *t = new double[nt];
	double *x = new double[nt];
	double *y = new double[nt];
	
	// Generate some random data points in parallel
	if(pID == 0)
	{
		cout<<"Generating random data...";	
		start = MPI_Wtime();
		random_device rd;
		mt19937 gen(rd());
		std::normal_distribution<> rnorm(0, noise_sd);
		std::uniform_real_distribution<> runif(0, tmax);
		
		for(unsigned int i=0; i<nt; i++)
		{
			t[i] = runif(gen);
			x[i] = rnorm(gen);
		} 
		cout<<"Done."<<endl;
	}
	else if(pID == 1)
	{
		random_device rd;
		mt19937 gen(rd());
		std::normal_distribution<> rnorm(0, noise_sd);
		
		for(unsigned int i=0; i<nt; i++) y[i] = rnorm(gen);
		cout<<"Bla"<<endl;
	}
	MPI_Bcast(t, nt, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(x, nt, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(y, nt, MPI_DOUBLE, 1, MPI_COMM_WORLD);
	// ============================================================
	
	// Assemble data set to look like the one in the odevo programs
	if( pID==0 ) cout<<"Assembling data...";
	fit_time.insert(fit_time.begin(), t, t+nt);
	sort( fit_time.begin(), fit_time.end() );
	vector<vector<double> > tc = solver.timecourse(fit_time);
	
	for(unsigned int i=0; i<nt; i++)
	{
		data_means.push_back(vector<double>(2));
		data_means[i][0] = x[i] + tc[i][0];
		data_means[i][1] = y[i] + tc[i][1];
		cout<<data_means[i][0]<<" "<<data_means[i][1]<<" "<<tc[i][0]<<endl;;
	}
	
	fit_idx.emplace_back(0);
	fit_idx.emplace_back(1);
	data_vars = vector<vector<double> >(nt, vector<double>(2, 0.0) );
	
	delete[] t, x, y;
	if( pID==0 ) cout<<"Done."<<endl;
	// ============================================================
	
	// Initialize the fit
	mpi_sres opt(fitness, np, 0);
	
	opt.toFile = 0;
	opt.logFreq = 100;
	
	for(unsigned int i=0; i<np; i++) 
	{
		opt.lb[i] = 0.0;
		opt.ub[i] = 10.0;
	}
	opt.init();
		
	opt+=1000;	

	stop = MPI_Wtime();
	
	if(opt.converged() && pID==0) std::cout<<"Converged!"<<std::endl;	
	
	if(pID == 0) std::cout<<"Needed "<<stop-start<<" s."<<std::endl;	
	
	MPI_Finalize();
	
	return 0;
}

 

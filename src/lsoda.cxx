/*
 * lsoda.cxx
 * This file is part of odevo
 *
 * Copyright (C) 2014 - Christian Diener
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

#include <fstream>
#include "lsoda.h"

// The error codes for the ode solver
const char* solver_errors[] = {
	"excessive amount of work done (more than 500 steps)",
	"too much accuracy was requested for the precision",
	"illegal input was detected",
	"repeated error test failures on one attempted step",
	"repeated convergence test failures on one attempted step",
	"error weight became zero, did you set ATOL=0.0?",
	"the length of RWORK and/or IWORK was too small, you should not see that :("
};

// Initializes all the stuff for LSODA
// call everytime you change the system
void lsoda::init()
{
	itol = 1; //scalar Tolerances
	rtol = RTOL;
	atol = ATOL;
	itask = 1;
	istate = 1;
	iopt = 1;
	jt = (use_jac) ? 1 : 2;
	lrw = 22 + neq * std::max(16, neq + 9);
    liw = neq + 20;
	
	// Check if we initialized before and clean up 
	if(!init_once)
	{
		rwork = new double[lrw];
		iwork = new int[liw];
		y = new double[neq];
	}
	
	for(unsigned int i=4; i<=9; i++)
	{
        	rwork[i] = 0.0;
        	iwork[i] = 0;
    	}
	iwork[5] = 1e5; //number of max internal iterations
	iwork[6] = 1;
	t = 0.0;
	
	for(unsigned int i=0; i<neq; i++) y[i] = y0[i];
	init_once = 1;
}

// Minimal constructor, you need to set sys before calling the solver
lsoda::lsoda(int n, unsigned int use_jac)
{
	init_once = 0;
	neq = n;
	this->use_jac = use_jac;
	init(); 
}


lsoda::lsoda(void (*sys_f)(int *neq, double *t, double *y, double *ydot), \
		void (*sys_jac)(int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd), \
		double* pars, double* ystart, int neq, int np, unsigned int use_jac)
{
	init_once = 0;
	f = sys_f;
	jac = sys_jac;
	this->neq = neq;
	this->np = np;
	p = pars;
	y0 = ystart;
	this->use_jac = use_jac;
	init();
}

lsoda::~lsoda()
{
	// Check if we initialized before and clean up 
	if(init_once)
	{
		delete[] rwork;
		delete[] iwork;
		delete[] y;
	}
}


int lsoda::step(const double to)
{
	tout = to;	

	dlsoda_(f, &neq, y, &t, &tout, &itol, &rtol, &atol, &itask, &istate, &iopt, \
                rwork, &lrw, iwork, &liw, jac, &jt);

	if(istate < 1)
	{
		int state = -(istate+1); 
		std::cout<<"Solver error: "<<solver_errors[state]<<std::endl;
		return istate;
	} 
	t = tout;
	return istate;
}

std::vector<std::vector<double> > lsoda::timecourse(const std::vector<double> to)
{
	std::vector<std::vector<double> > out;
	
	init();
	for(unsigned int i=0; i<to.size(); i++)
	{
		if(to[i] <= t)
		{
			std::cout<<"Output times must be strictly monotonously increasing and larger than 0.0!\n"<<std::endl;
			return out;
		}
		
		step(to[i]);
		
		out.push_back( std::vector<double>(y, y+neq-1) );
	}
	return out;
}

int lsoda::timecourse_to_file(const std::vector<double> to, const char* out_file)
{
	std::ofstream ofs(out_file);
	
	ofs<<"Time";
	for(unsigned int i=1; i<=neq; i++) ofs<<"\t\ty"<<i<<" ";
	ofs<<std::endl;

	init();
	for(unsigned int i=0; i<to.size(); i++)
	{
		if(to[i] < t)
		{
			std::cout<<"Output times must be strictly monotonously increasing and larger than 0.0!\n"<<std::endl;
			return istate;
		}
		
		step(to[i]);
		ofs<<t;
		for(unsigned int j=1; j<=neq; j++) ofs<<"\t\t"<<y[j-1]<<" ";
		ofs<<std::endl;
	}
	
	ofs.close();

	return istate;
}

std::ostream& operator<<(std::ostream& out, lsoda& ls)
{
	out<<"lsoda solver diagnostics:"<<std::endl;
	out<<"> step size = "<<ls.rwork[10]<<" | "<<"t = "<<ls.rwork[12]<<" | "<<"f evals = "<<ls.iwork[11]<<std::endl;
	
	if(ls.iwork[19]==2) 
	{
		out<<"> using the bdf method (last switch at t="
		   <<ls.rwork[15]<<")."<<std::endl;
	}
	else 
	{
		out<<"> using the adams method (last switch t="<<ls.rwork[15]
		   <<")."<<std::endl;
	}

	return out;
}

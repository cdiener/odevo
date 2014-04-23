/*
 * lsoda.h
 * This file is part of optim
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

#ifndef __LSODA_H__
#define __LSODA_H__


#include <vector>
#include <iostream>

#define ATOL 1e-12
#define RTOL 1e-8

// External call masking the FORTRAN implementation
extern "C"
{
	extern void dlsoda_(void f(int *neq, double *t, double *y, double *ydot),
        int *neq, double *y, double *t, double *tout,
        int *itol, double *rtol, double *atol, int *itask, int *istate,
        int *iopt, double *rwork, int *lrw, int *iwork, int *liw,
        void jac(int *neq, double *t, double *y, int *ml,
        int *mu, double *pd, int *nrowpd), int *jt);
}

// The "model" (ode system) and optional jacobian
typedef void (*Model)(int*, double*, double*, double*);
typedef void (*Jacobian)(int*, double*, double*, int*, int*, double*, int*);

/***
 * The lsoda solver class describing an adaptive fast solver with automatic
 * stiff/non-stiff switching.
 */ 
class Lsoda
{
	private:
		// Solver properties
		double* y;
		int itol; //whether absolute errors are scalar (1) or per component (2)
		int itask;
		int istate;
		int iopt;
		int lrw;
		int liw;
		int jt;
		double* rwork;
		int* iwork;
		double t;
		double tout;	
		unsigned int init_once;			
		
	public:
		// Model properties
		Model f;
		Jacobian jac;
		double* p;
		double* y0;
		double rtol;
		double atol;
		int neq, np;
		unsigned int use_jac;

		// constructors
		lsoda(int n, unsigned int use_jac);
		
		lsoda(Model sys_f, Jacobian sys_jac, double* pars, double* ystart, int neq, \
			int np, unsigned int use_jac);
				
		~lsoda();
		
		// helper functions
		void init();
		double* get_y() { return y; }
		double get_t() { return t; }
		
		// solver functions
		int step(const double tout);
		std::vector<std::vector<double> > 
			timecourse(const double* to, const unsigned int n_t);
		int timecourse_to_file(const double* to, const unsigned int n_t, 
									const char* out_file);	

		// output functions
		friend std::ostream& operator<<(std::ostream& out, Lsoda& ls);	
};


#endif /* __LSODA_H__ */

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

// default error tolerances
const double ATOL  = 1e-12;
const double RTOL = 1e-8;

// External call masking the FORTRAN lsoda implementation
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
class lsoda
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
		// solver settings
		Model f;				// the ode system to use
		Jacobian jac;			// the (optional) jacobian of the system
		double* p;				// the parameters of the ode system
		double* y0;				// the initial values of the ode system
		double rtol;			// the relative error tolerance
		double atol;			// the absolute error tolerance
		int neq, np;			// number of equations and parameters
		unsigned int use_jac;	// whether to use a jacobian

		// constructors
		
		/***
		 * Basic empty constructor. You will still have to set most of the 
		 * solver settings by hand.
		 * 
		 * @param n The number of equations.
		 * @param use_jac Whether the jacobian is used or no.
		 */ 
		lsoda(int n, unsigned int use_jac);
		
		/***
		 * Complete constructor.
		 * 
		 * @param sys_f The ode system f(y,t) to be used
		 * @param sys_jac The jacobian df/dy to be used. Can be an empty function.
		 * @param pars The parameters of the system.
		 * @param ystart The initial y values.
		 * @param neq The number of equations in the system.
		 * @param np The number of parameters in the system.
		 * @param use_jac Whether the jacobain is used or approximated internally.
		 */ 
		lsoda(Model sys_f, Jacobian sys_jac, double* pars, double* ystart, 
				int neq, int np, unsigned int use_jac);
		
		/***
		 * The destructor.
		 */ 		
		~lsoda();
		
		// helper functions
		
		/***
		 * Initializes the solver. Must be called before stepping. Only the
		 * complete contructor calls this function automatically. The function 
		 * must be recalled when changing the system.
		 */
		void init();
		
		/***
		 * Gets the current y vakues
		 */ 
		double* get_y() { return y; }
		
		/***
		 * Gets the current time.
		 */ 
		double get_t() { return t; }
		
		// solver functions
		
		/***
		 * Steps a single time step.
		 * 
		 * @param tout The next time to step to. Must be larger than current t.
		 * @return The solver state: >0 = success, <0 = failure.
		 */ 
		int step(const double tout);
		
		/***
		 * Gets an entire timecourse.
		 * 
		 * @param to A vector of strictly increasing time points (t0 = 0.0).
		 * @return A matrix of outputs mit M[T][Y]
		 */ 
		std::vector<std::vector<double> > timecourse(const std::vector<double> to);
		
		/*** 
		 * Solve a timecourse and save it directly to a tab delimited file.
		 * 
		 * @param to A vector of strictly increasing time points (t0 = 0.0).
		 * @param out_file Name of the output text file.
		 * @return The solver state: >0 = success, <0 = failure.
		 */ 
		int timecourse_to_file(const std::vector<double> to, const char* out_file);	

		// output functions
		
		/***
		 * Gives some nice solver diagnostics.
		 */
		friend std::ostream& operator<<(std::ostream& out, lsoda& ls);	
};


#endif /* __LSODA_H__ */

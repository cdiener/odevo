/*
 * optim.h
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

#ifndef __OPTIM_H__
#define __OPTIM_H__

#include <fstream>
#include <iostream>
#include <vector>
#include <mpi.h>

#include "SRES/sharefunc.h"
#include "SRES/ESSRSort.h"
#include "SRES/ESES.h"

typedef void (*FitnessFunction)(double*, double*, double*);

class mpi_sres
{
	private:
		FitnessFunction f;

		ESParameter *param;
		ESPopulation *population;
		ESStatistics *stats;

		int n_const, n_p, gen_max;
		unsigned int seed;

		std::ofstream out_file;	

		double* strat_stat();	
		
	public:
		// settings for evolutionary algorithm		
		ESfcnTrsfm *trsfm;					// parameter transformation function
		int es, miu, lambda, retry;			// integer algorithm parameters
		double gamma, alpha, varphi, pf;	// double algorithm parameters
		int toFile, logFreq;				// whether logging to file and how often
		char* logFile;						// name of the log file
		double *ub, *lb;					// upper and lower parameter bounds
		
		// Parameters for convergence check
		double rtol, atol, stol;
		int gen_min;

		//Convergence tester
		int converged();
		
		
		//Stepping functions
		mpi_sres& operator++();
		friend mpi_sres& operator+=(mpi_sres& opt, int steps);

		//Initializer
		void init();

		//Nice output functions
		friend std::ostream& operator<<(std::ostream& out, mpi_sres& sres);

		//Con- and Destructors
		mpi_sres(FitnessFunction ff, int n_p, int init_now = 1, int n_const = 1, int gmax = 1e5);
		~mpi_sres();
};

#endif /* __OPTIM_H__ */


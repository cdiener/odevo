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
#include <mpi.h>

#include "SRES/sharefunc.h"
#include "SRES/ESSRSort.h"
#include "SRES/ESES.h"

//Maximal parameter columns in terminal output
#define MAXCOL 5
#define SEP "--------------------------------------------------------------------------------"

typedef void (*FitnessFunction)(double*, double*, double*);

class MPI_SRES
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
		//Parameters for evolutionary algorithm		
		ESfcnTrsfm *trsfm;
		int es, miu, lambda, retry;
		double gamma, alpha, varphi, pf;
		int toFile, logFreq;
		char* logFile;
		double *ub, *lb;
		
		// Parameters for convergence check
		double rtol, atol, stol;
		int gen_min;

		//Convergence tester
		int converged();
		
		
		//Stepping functions
		MPI_SRES& operator++();
		friend MPI_SRES& operator+=(MPI_SRES& opt, int steps);

		//Initializer
		void init();

		//Nice output functions
		friend std::ostream& operator<<(std::ostream& out, MPI_SRES& sres);

		//Con- and Destructors
		MPI_SRES(FitnessFunction ff, int n_p, int init_now = 1, int n_const = 1, int gmax = 1e5);
		~MPI_SRES();
};

#endif /* __OPTIM_H__ */


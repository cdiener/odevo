/*
 * optim.cxx
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

#include "optim.h"
#include <cmath>

void MPI_SRES::init()
{
	int ac = 0;
	char ** av = NULL;

	if(toFile)
	{ 
		out_file.open(logFile);
		out_file<<"Gen\t\tF";
		for(unsigned int i=0; i<n_p; i++) out_file<<"\t\tp"<<i+1;
		out_file<<std::endl;
	}

	ESInitial(&ac, &av, seed, &param, trsfm, f, es, n_const, n_p, ub, lb, miu, lambda, \
			gen_max, gamma, alpha, varphi, retry, &population, &stats);
}

double* MPI_SRES::strat_stat()
{
	double* stat = new double[2];
	stat[0] = stat[1] = 0.0;
	
	for(unsigned int i; i<n_p; i++) stat[0] += stats->bestindvdl->sp[i];
	for(unsigned int i; i<n_p; i++) stat[1] += pow(stats->bestindvdl->sp[i] - stat[0], 2);

	stat[0] = stat[0]/(double)n_p;
	stat[1] = sqrt( stat[1]/(n_p - 1.0) );

	return stat;
}

MPI_SRES::MPI_SRES(FitnessFunction ff, int n_p, int init_now, int n_const, int gmax)
{
	this->n_p = n_p;
	this->gen_max = gmax;
	this->n_const = n_const;
	this->f = ff;
	trsfm = NULL;
	
	this->lb = new double[n_p];
	for(unsigned int i=0;i<n_p;i++) lb[i] = 0.0;

	this->ub = new double[n_p];
	for(unsigned int i=0;i<n_p;i++) ub[i] = 1.0e6;
		
	
	// Choose standard parameters	
	seed = shareDefSeed;
	gamma = esDefGamma;
	alpha = esDefAlpha;
	varphi = esDefVarphi;
	retry = esDefRetry;
	miu = 50;
	lambda = 350;

	//Set log behaviour
	toFile = 1;
	logFreq = 100;
	logFile = (char*)"optim.log";

	//Set convergence tests
	rtol = 1.0e-6;
	atol = 1.0e-8;
	stol = 1.0e-6;
	gen_min = 10;
	
	if(init_now) init();
}

MPI_SRES::~MPI_SRES()
{
	delete[] lb;
	delete[] ub;
	
	if(out_file.is_open()) out_file.close();
	
	ESDeInitial(param, population, stats);
}

MPI_SRES& MPI_SRES::operator++()
{
	if(stats->curgen < param->gen) ESStep(population, param, stats, pf);

	return *this;
}

MPI_SRES& operator+=(MPI_SRES& opt, int steps)
{
	unsigned int j;

	int myid;

	if(opt.stats->curgen < opt.param->gen) 
	for(unsigned int i=0; i<steps; i++)
	{
		++opt;
		ESIndividual* best = opt.stats->bestindvdl;
		
		//Print logging information
		if(opt.stats->curgen % opt.logFreq == 0)
		{	
			MPI_Comm_rank(MPI_COMM_WORLD, &myid);
			if(opt.toFile && myid==0)
			{
				opt.out_file<<opt.stats->curgen<<"\t\t"<<best->f;
				ESfcnTrsfm* trsfm = opt.param->trsfm;
				
				
				if(trsfm == NULL) 
				for(j=0; j<opt.n_p; j++) opt.out_file<<"\t\t"<<best->op[j];
				else
				for(j=0; j<opt.n_p; j++)
				{
					if(trsfm[j] == NULL)
						opt.out_file<<"\t\t"<<best->op[j];
					else
						opt.out_file<<"\t\t"<<(trsfm[j])(best->op[j]);
				}
				opt.out_file<<std::endl;
			}
			else std::cout<<opt;			
		}
	}		
	
	return opt;
}

int MPI_SRES::converged()
{
	double* strat = strat_stat();
	
	if(strat[0] > stol)
	{
		delete[] strat;
		return 0;
	}
	else
	{
		delete[] strat;
		return 1;
	}
}

std::ostream& operator<<(std::ostream& out, MPI_SRES& sres)
{
	int myid;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if(myid != 0) return out;	

	double* strat = sres.strat_stat();

	out<<SEP<<std::endl;
	out<<"Generation: "<<sres.stats->curgen<<" | ";
	out<<"f = "<<sres.stats->bestindvdl->f<<" (Gen. "<<sres.stats->bestgen<<") | ";
	out<<"Mean strategy = "<<strat[0]<<" +- "<<strat[1]<<std::endl;
	out<<SEP<<std::endl;

	for(unsigned int i=0; i<sres.n_p; i++)
	{
		out<<"p["<<i+1<<"] = "<<sres.stats->bestindvdl->op[i]<<"\t";
		if((i+1) % MAXCOL == 0 && sres.n_p > i+1) out<<std::endl;
	}
	out<<std::endl;

	delete[] strat;

	return out;
}
 

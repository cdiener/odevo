/*
 * system.h
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

#ifndef __SYSTEM_H__
#define __SYSTEM_H__

#include <cmath>

using namespace std;

double ystart[] = {0.0, 1.0};
double p[] = {42.0};
const int use_jac = 0;
const unsigned int np = 1;
const unsigned int neq = 2;

void sys(int *neq, double *t, double *y, double *ydot)
{
	ydot[0] = y[1];  
	ydot[1] = p[0]*(1.0 - y[0]*y[0])*y[1] - y[0];
}

void jac(int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd){}

#endif /* __SYSTEM_H__ */


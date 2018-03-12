/*
 * extras.c
 * 
 * Copyright 2011 Fernando Pujaico Rivera <fernando.pujaico.rivera@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#include <stdio.h>
#include "extras.h"
#include <math.h>

#include <pds/pdsba.h>
#include <pds/pdsra.h>
#include <pds/pdsmath.h>

#include <pds/pdsit.h>


double fpar2(double p1,double p2)
{
	return p1+p2-2*p1*p2;
}

double fpar3(double p1,double p2,double p3)
{
	return fpar2(fpar2(p1,p2),p3);
}

double fpar4(double p1,double p2,double p3,double p4)
{
	return fpar2(fpar3(p1,p2,p3),p4);
}

double hu0Omega3_metodo1(PdsVector *Ps)
{
	double H=0;
	double p1,p2,p3;

	if(Ps->Nel<3)	return -1;

	p1=Ps->V[0];
	p2=Ps->V[1];
	p3=Ps->V[2];


	H=	  pds_hb(p1)+pds_hb(p2)+pds_hb(p3) 
		- pds_hb(fpar2(p1,p2))- pds_hb(fpar2(p1,p3))- pds_hb(fpar2(p3,p2))
		+ pds_hb(fpar3(p1,p2,p3));

	return H;
}

double hu0Omega4_metodo1(PdsVector *Ps)
{
	double H=0;
	double p1,p2,p3,p4;
	
	if(Ps->Nel<4)	return -1;

	p1=Ps->V[0];
	p2=Ps->V[1];
	p3=Ps->V[2];
	p4=Ps->V[3];


	H=	  pds_hb(p1)+pds_hb(p2)+pds_hb(p3) +pds_hb(p4) 
		- pds_hb(fpar2(p1,p2))- pds_hb(fpar2(p1,p3))- pds_hb(fpar2(p1,p4))- pds_hb(fpar2(p2,p3))- pds_hb(fpar2(p2,p4))- pds_hb(fpar2(p3,p4))
		+ pds_hb(fpar3(p1,p2,p3))+ pds_hb(fpar3(p2,p3,p4))+ pds_hb(fpar3(p1,p2,p4))+ pds_hb(fpar3(p1,p3,p4))
		- pds_hb(fpar4(p1,p2,p3,p4));

	return H;
}

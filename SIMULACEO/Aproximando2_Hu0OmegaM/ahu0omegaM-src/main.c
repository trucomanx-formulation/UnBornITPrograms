/*
 * main.c
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
#include <math.h>
#include <pds/pdsdatafunc.h>
#include <pds/pdsba.h>
#include <pds/pdsra.h>
#include <pds/pdsrv.h>
#include <pds/pdsmath.h>

#include "extras.h"

#include <pds/pdsit.h>

#define	DEFAULT_PS_INIT		0.3

#define P0 0.5

int main(int argc, char** argv)
{
	//PdsBaNatural M=DEFAULT_M;
	double rho=DEFAULT_PS_INIT;
	int i;

	int L=1000;

	double hcond,hcond2,ERROR;


	PdsVector *Ps=NULL;
	PdsUniform *RV=NULL;
	


	////////////////////////////////////////////////////////////////////////
	pds_get_double_param (argc, argv,"-p",&rho);
	
    float low_rho=0.01;
    float high_rho=rho;
	RV=pds_uniform_new(low_rho,high_rho);

    int Nrho=3;
	printf("\nrho[%d]:<%f,%f>\n",Nrho,low_rho,high_rho);

	Ps=pds_vector_new(Nrho);

	ERROR=0;
	for(i=0;i<L;i++)
	{
		pds_uniform_get_vector(RV,Ps);
		pds_entropy_u0_omega_bsc_model(Ps,P0, &hcond);

		hcond2=hu0Omega3_metodo1(Ps);
		ERROR=ERROR+fabs(hcond-hcond2)/hcond;

        printf("working %3d of %3d\r",i,L);
        fflush(stdout);
	}

	printf("ERROR percent: %f%c\n",100*ERROR/L,'%');

    printf("\n<<<< Last test >>>>\n");
    printf("Ps:\n");
    pds_vector_printf(Ps);
    printf("H(U0|Omega)[exact]:%f\n",hcond);
    printf("H(U0|Omega)[aprox]:%f\n",hcond2);
    printf("\n");
    
	return 0;
}





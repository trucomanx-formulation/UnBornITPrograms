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
#include <pds/pdsmath.h>
#include <pds/pdsit.h>
#include <pds/pdsrv.h>

#include "extras.h"

#define	DEFAULT_M			8
#define	DEFAULT_PS_INIT		0.1

#define PERCENT 			100

int main(int argc, char** argv)
{
	PdsBaNatural M=DEFAULT_M;
	int tmp,i;
	PdsVector *Ps=NULL;

	PdsUniform *RV=NULL;

	double rho=DEFAULT_PS_INIT,p0=0.5,teorho;

	double hjoint,fhjoint;
	PdsRaReal m;
	double Sumh,Pondh;
	double Psa,HPsa;


	////////////////////////////////////////////////////////////////////////
	tmp=M;
	pds_get_int_param (argc, argv,"-m",&tmp);
	M=(unsigned short)tmp;
	pds_get_double_param (argc, argv,"-p",&rho);
	

	printf("M: %d\n",M);
	printf("rho:\t%f\n",rho);

	////////////////////////////////////////////////////////////////////////
	Ps=pds_vector_new(M);
	//pds_vector_init_value (Ps,rho);

	RV=pds_uniform_new(0.001,0.002);
	pds_uniform_get_vector (RV,Ps);
	Sumh=0,Pondh=0;
	for(i=M-1;i>=0;i--)
	{
		Sumh=Sumh+pds_hb(Ps->V[i]);
		Pondh=Pondh+Ps->V[i]*pds_hb(Ps->V[i]);
	}

	HPsa=Sumh/M;
	Psa =pds_hbinv(HPsa);

	pds_vector_mean_vector(Ps, &m);
	printf("Ps:\t"); pds_vector_printf(Ps); 
	printf("E[Ps]           =\t%f\n",m);
	printf("E[hb(Ps)]  =HPsa=\t%f\n",HPsa);
	printf("hbinv(HPsa)=Psa =\t%f\n",Psa);
	printf("E[Psi,hb(Psi)]  =\t%f\n",Pondh/Sumh);
	printf("\n");

	////////////////////////////////////////////////////////////////////////
	pds_joint_entropy_bsc_model(Ps,p0, &hjoint);
	printf("Teor H(u1..uM)=\t%f\n",hjoint);
	pds_inv_symetric_joint_entropy_bsc_model (hjoint,M,&teorho);
	printf("Teor rho      =\t%e\n",teorho);
	printf("\n");

	////////////////////////////////////////////////////////////////////////
	pds_symetric_joint_entropy_bsc_model(Psa,M,&fhjoint);
	printf("     H(Psa,M) =\t%f\t%6.2f percent\n",fhjoint,(fhjoint-hjoint)*PERCENT/hjoint);


	////////////////////////////////////////////////////////////////////////
	printf("1+ Sum hb(pi) =\t%f\t%6.2f percent\n",Sumh+1,(Sumh+1-hjoint)*PERCENT/hjoint);
	////////////////////////////////////////////////////////////////////////
	pds_fast_joint_entropy_bsc_model(Ps,&fhjoint);
	printf("Fast H(u1..uM)=\t%f\t%6.2f percent\n",fhjoint,(fhjoint-hjoint)*PERCENT/hjoint);
	printf("========================================\n");

	////////////////////////////////////////////////////////////////////////
	//pds_fast2_joint_entropy_bsc_model(Ps,&fhjoint);
	for(i=0;i<M;i++)
	{
		//pds_ordenap2a(Ps);
		//pds_ordenap2(Ps);
		pds_ordenap_hu0rho2(Ps);
	}

	printf("PsX:\t"); pds_vector_printf(Ps); 
	pds_vector_mean_vector(Ps, &m);
	printf("E[PsX]= PsXm          =\t%f\n",m);
	////////////////////////////////////////////////////////////////////////
	pds_symetric_joint_entropy_bsc_model(m,M,&fhjoint);
	printf("H(PsXm,M)             =\t%f\t%6.2f percent\n",fhjoint,(fhjoint-hjoint)*PERCENT/hjoint);
	pds_symetric_entropy_u0_omega_bsc_model (m,M,&fhjoint);
	printf("1+Shb(pi)-H(u0|rho,M) =\t%f\t%6.2f percent\n",Sumh+1-fhjoint,(Sumh+1-fhjoint-hjoint)*PERCENT/hjoint);
	printf("\n");


	return 0;
}





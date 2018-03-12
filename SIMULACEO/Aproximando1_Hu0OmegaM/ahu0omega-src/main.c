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

#define	DEFAULT_M		6
#define	DEFAULT_PS_INIT		0.3

#define MILLON 			1000000

int main(int argc, char** argv)
{
	PdsBaNatural M=DEFAULT_M;
	int tmp,i;
	PdsVector *Ps=NULL;
	PdsVector *PsX=NULL;

	PdsUniform *RV;

	double rho=DEFAULT_PS_INIT,p0=0.5;

	double hcond,hcond2,hjoint;
	PdsRaReal ps_mean, ps_weight, psx_mean, psx_var;
	double ps_teori;
	double Sum_Hb_Psi,Sum_Psi_Hb_Psi,ps_invEHb;

	////////////////////////////////////////////////////////////////////////
	tmp=M;
	pds_get_int_param (argc, argv,"-m",&tmp);
	M=(unsigned short)tmp;
	pds_get_double_param (argc, argv,"-p",&rho);
	

	printf("M: %d\n",M);
	printf("rho:\t%e\n",rho);


	RV=pds_uniform_new(0.25-rho/2,0.25+rho/2);
	////////////////////////////////////////////////////////////////////////
	Ps=pds_vector_new(M);
	//pds_vector_init_value (Ps,rho);
	pds_uniform_get_vector (RV,Ps);

	Sum_Hb_Psi=0,
	Sum_Psi_Hb_Psi=0;
	for(i=M-1;i>=0;i--)
	{
		//Ps->V[i]=rho-i*rho/(2*M);
		Sum_Hb_Psi    =Sum_Hb_Psi     +          pds_hb(Ps->V[i]);
		Sum_Psi_Hb_Psi=Sum_Psi_Hb_Psi + Ps->V[i]*pds_hb(Ps->V[i]);
	}
	PsX=pds_vector_new_vector(Ps);

	printf("Calculando H(U0|OmegaM) ...\n");
	pds_entropy_u0_omega_bsc_model(Ps,p0, &hcond);
	//printf("H(U0|OmegaM) = %f\n",hcond);
	printf("Calculando rho of H(U0|rho,M) ...\n");
	pds_inv_symetric_entropy_u0_omega_bsc_model(hcond,Ps->Nel,&ps_teori);
	//printf("rho of H(U0|rho,M) = %f\n",ps_teori);

	pds_vector_mean_vector(Ps, &ps_mean);
	ps_invEHb = pds_hbinv(Sum_Hb_Psi/M);
	ps_weight = Sum_Psi_Hb_Psi/Sum_Hb_Psi;

	printf("Ps:\t"); pds_vector_printf(Ps); 
	printf("ps_teori  = \t%e <-- ## Teorico ##\n",ps_teori);
	printf("ps_invEHb = \t%e <-- hbinv(SUM{hb(Psi)}/M)\n",ps_invEHb);
	printf("ps_mean   = \t%e <-- SUM{Psi}/M \n",ps_mean);
	printf("ps_weight = \t%e <-- SUM{Psi*hb(Psi)}/SUM{hb(Psi)}\n",ps_weight);
	printf("\n");

	////////////////////////////////////////////////////////////////////////
	printf("Original H(u0|u1..uM) =\t%e\n",hcond);

	////////////////////////////////////////////////////////////////////////
	printf("Calculando H(OmegaM)...\n");
	pds_joint_entropy_bsc_model(Ps,p0,&hjoint);
	printf("Calculado....[OK]\n");
	hcond2=Sum_Hb_Psi+1-hjoint;
	printf("Sum{H(pi)}+1-H(Omega) =\t%e\t%6.2f %c\n",hcond2,(hcond2-hcond)*100/hcond,'%');


	////////////////////////////////////////////////////////////////////////
	pds_symetric_entropy_u0_omega_bsc_model(ps_teori,M,&hcond2);
	printf("Symetric H(u0|u1..uM) =\t%e\t%6.2f %c <-- ps_teori  \n",hcond2,(hcond2-hcond)*100/hcond,'%');

	////////////////////////////////////////////////////////////////////////
	pds_symetric_entropy_u0_omega_bsc_model(ps_invEHb,M,&hcond2);
	printf("Symetric H(u0|u1..uM) =\t%e\t%6.2f %c <-- ps_invEHb  \n",hcond2,(hcond2-hcond)*100/hcond,'%');

	////////////////////////////////////////////////////////////////////////
	pds_symetric_entropy_u0_omega_bsc_model(ps_mean,M,&hcond2);
	printf("Symetric H(u0|u1..uM) =\t%e\t%6.2f %c <-- ps_mean  \n",hcond2,(hcond2-hcond)*100/hcond,'%');

	////////////////////////////////////////////////////////////////////////
	pds_symetric_entropy_u0_omega_bsc_model(ps_weight,M,&hcond2);
	printf("Symetric H(u0|u1..uM) =\t%e\t%6.2f %c <-- ps_weight  \n",hcond2,(hcond2-hcond)*100/hcond,'%');


	////////////////////////////////////////////////////////////////////////
	printf("\n");
	for(i=0;i<M;i++)
	{
		//pds_ordenap_hu0rho2(PsX);
		pds_ordenap_hu0rhoN(PsX,1);
		//pds_ordenap_hu0rho4(PsX);
	}

	//printf("Ps :\t"); pds_vector_printf(Ps); 
	//printf("PsX:\t"); pds_vector_printf(PsX); 
	pds_vector_mean_vector(PsX, &psx_mean);
	pds_vector_var_vector (PsX, &psx_var);
	printf("    E[PsX]=\t%e\n",psx_mean);
	printf("SIGMA[PsX]=\t%e\n",sqrt(psx_var));
	printf("\n");

	pds_symetric_entropy_u0_omega_bsc_model(psx_mean,M,&hcond2);
	printf("Symetric H(u0|u1..uM) =\t%e\t%6.2f %c <-- E[PsX]  \n",hcond2,(hcond2-hcond)*100/hcond,'%');
	printf("\n");


	printf("\n");

	return 0;
}





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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/** \fn int pds_pomega(const PdsVector *Ps, double *pOmega)
 *  \brief Encuentra 
 *  \param[in] Ps Es el vector de probabilidades de error de las fuentes BSC.
 *  \param[out] pOmega Porbabilidade POmega.
 *  \return TRUE si todo fue bien, o FALSE sino. Ejem: Ps==NULL.
 *  \ingroup PdsBVectorGroup
 */
int pds_pomega(const PdsVector *Ps, double *pOmega)
{
	PdsBaNatural i;
	double p;

	if(Ps==NULL)		return FALSE;

	p=Ps->V[0];
	for(i=1;i<Ps->Nel;i++)
	{
		p=pds_hbinv(pds_hb(Ps->V[i])+pds_hb(p)-pds_hb(p + Ps->V[i] - 2*Ps->V[i]*p));
	}
	*pOmega=p;

	return TRUE;
}

/** \fn int pds_fast_joint_probability_bsc_model(const PdsVector *Ps, double *prob)
 *  \brief Encuentra la entropia conjunta para M fontes generadas
 *  pasando una fuente U0, con probabilidade Pr(U0=1)=0.5, atraves de M canales BSC
 *  con probabilidades de error Ps=[Ps0, Ps1, ..., Ps(M-1)].
 *  \param[in] Ps Es el vector de probabilidades de error de las fuentes BSC.
 *  \param[out] H Entropia conjuta.
 *  \return TRUE si todo fue bien, o FALSE sino. Ejem: Ps==NULL.
 *  \ingroup PdsBVectorGroup
 */
int pds_fast_joint_entropy_bsc_model(const PdsVector *Ps, double *H)
{
	PdsBaNatural i;
	double hjoint,pOmega;

	if(Ps==NULL)		return FALSE;

	hjoint=1;

	for(i=0;i<Ps->Nel;i++)
	{
		hjoint=hjoint+pds_hb(Ps->V[i]);
	}

	pOmega=0;

	pds_pomega(Ps,&pOmega);

	hjoint=hjoint-pds_hb(pOmega);

	*H=hjoint;
	return TRUE;
}


/** \fn int pds_fast2_joint_probability_bsc_model(const PdsVector *Ps, double *prob)
 *  \brief Encuentra la entropia conjunta para M fontes generadas
 *  pasando una fuente U0, con probabilidade Pr(U0=1)=0.5, atraves de M canales BSC
 *  con probabilidades de error Ps=[Ps0, Ps1, ..., Ps(M-1)].
 *  \param[in] Ps Es el vector de probabilidades de error de las fuentes BSC.
 *  \param[out] H Entropia conjuta.
 *  \return TRUE si todo fue bien, o FALSE sino. Ejem: Ps==NULL.
 *  \ingroup PdsBVectorGroup
 */
int pds_fast2_joint_entropy_bsc_model(const PdsVector *Ps, double *H)
{
	PdsBaNatural i,m;
	double hm,S,pm;

	if(Ps==NULL)		return FALSE;

	pm=Ps->V[0];
	hm=1;

	for(m=1;m<Ps->Nel;m++)
	{
		for(i=0,S=0;i<=m;i++)
		{
			S=S+pds_hb(Ps->V[i]);
		}
		hm= 1 + S-( pds_hb(Ps->V[m])+pds_hb(pm)-pds_hb(pm + Ps->V[m] - 2*Ps->V[m]*pm) );
		pm=pds_hbinv(1+S-hm);
	}
	*H=hm;
	return TRUE;
}

int pds_ordenap(PdsVector *Ps)
{
	PdsRaNatural i,M;
	PdsRaReal p,par;

	if(Ps==NULL)	return FALSE;

	M=Ps->Nel;

	pds_vector_asc_order (Ps);
	
	for(i=0;i<M/2;i++)
	{
		par=Ps->V[i]+Ps->V[M-1-i]-2.0*Ps->V[i]*Ps->V[M-1-i];
		if(par>0.5) par=1-par;

		p=0.5-0.5*sqrt(1-2*par);
		
		Ps->V[i]=p;
		Ps->V[M-1-i]=p;
	}
	
	return TRUE;
}

int pds_ordenap2(PdsVector *Ps)
{
	PdsRaNatural i,M;
	PdsRaReal p;

	if(Ps==NULL)	return FALSE;

	M=Ps->Nel;

	pds_vector_asc_order (Ps);

	for(i=0;i<M/2;i++)
	{
		p=pds_hbinv(0.5*pds_hb(Ps->V[i])+0.5*pds_hb(Ps->V[M-1-i]));
			
		Ps->V[i]=p;
		Ps->V[M-1-i]=p;
	}
	
	return TRUE;
}



int pds_ordenap2a(PdsVector *Ps)
{
	PdsRaNatural i,M;
	double p,h;

	if(Ps==NULL)	return FALSE;

	M=Ps->Nel;

	pds_vector_asc_order (Ps);

	for(i=0;i<M/2;i++)
	{
		h=pds_hb( pds_hb(Ps->V[i])+pds_hb(Ps->V[M-1-i])-2*pds_hb(Ps->V[i])*pds_hb(Ps->V[M-1-i]) )+1;
		pds_inv_symetric_joint_entropy_bsc_model (h,2,&p);
			
		Ps->V[i]=p;
		Ps->V[M-1-i]=p;
	}
	
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

double func_hu0u1u2(double p)
{
	return 2*pds_hb(p)-pds_hb(2*p*(1-p));
}

double func_hu0u1u2_inv(double h)
{
	double p_min,p_max,p0,h0,dh;
	if(h>=1.0)	return 0.5;
	if(h<=0.0)	return 0.0;

	p_min=0.0;	p_max=0.5;	

	dh=0.01;
	
	do{
		p0=(p_min+p_max)/2.0;

		h0=func_hu0u1u2(p0);

		if(h>h0)	p_min=p0;
		else		p_max=p0;
	}while( (fabs(h-h0)/h) > (dh/100.0) );

	return p0;
}

int pds_ordenap_hu0rho2(PdsVector *Ps)
{
	PdsRaNatural i,M;
	PdsRaReal p,k,par,p1,p2;

	if(Ps==NULL)	return FALSE;

	M=Ps->Nel;

	//printf("Ps:"); pds_vector_printf (Ps);
	pds_vector_asc_order (Ps);
	//printf("Ps:"); pds_vector_printf (Ps);

	for(i=0;i<M/2;i++)
	{
		p1=Ps->V[i];
		p2=Ps->V[M-1-i];
		par=p1+p2-2*p1*p2;

		k=pds_hb(p1)+pds_hb(p2)-pds_hb(par);

		p=func_hu0u1u2_inv(k);
			
		Ps->V[i]=p;
		Ps->V[M-1-i]=p;
	}
	
	return TRUE;
}

int pds_ordenap_hu0rho4(PdsVector *Ps)
{
	PdsRaNatural i,M;
	double p;
	double hcond;

	PdsVector *P=NULL;

	if(Ps==NULL)	return FALSE;

	M=Ps->Nel;

	P=pds_vector_new(4);

	pds_vector_asc_order (Ps);

	for(i=0;i<(M/2-1);i=i+2)
	{

		P->V[0] = Ps->V[i];
		P->V[1] = Ps->V[i+1];
		P->V[2] = Ps->V[M-1-i-1];
		P->V[3] = Ps->V[M-1-i];


		pds_entropy_u0_omega_bsc_model(P,0.5, &hcond);

		pds_inv_symetric_entropy_u0_omega_bsc_model(hcond,P->Nel,&p);

		Ps->V[i]       = p;
		Ps->V[i+1]     = p;
		Ps->V[M-1-i-1] = p;
		Ps->V[M-1-i]   = p;
	}

	pds_vector_free(P);
	
	return TRUE;
}

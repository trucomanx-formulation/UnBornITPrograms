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
#include "extrasbvector.h"
#include "extrascoin.h"
#include "extras.h"

#include <pds/pdsrv.h>
#include <pds/pdsra.h>
#include <pds/pdsba.h>

PdsBVector **pds_bvectors_new(PdsBaNatural m,PdsBaNatural k)
{
	PdsBVector **V=NULL;
	PdsBaNatural i,j;
	
	V=(PdsBVector **)calloc(m,sizeof(PdsBVector *));
	if(V==NULL)	return NULL;

	for(i=0;i<m;i++)
	{
		V[i]=pds_bvector_new(k);
		if(V[i]==NULL)
		{
			for(j=0;j<i;j++) pds_bvector_free(V[j]);
			return NULL;
		}
	}

	return V;
}


int pds_bvectors_load_coins(PdsBVector **u,PdsCoin **c,PdsBaNatural m)
{
	int d;
	PdsRaNatural i;
	
	for(i=0;i<m;i++)
	{
		d=pds_coin_get_bvector(c[i],u[i]);
		if(d==FALSE)	return FALSE;
	}

	return TRUE;
}



int pds_bvectors_xor(PdsBVector **U,PdsBVector **E,PdsBaNatural m,PdsBVector *u)
{
	PdsBaNatural i;
	int d;

	if(U==NULL)			return FALSE;
	if(E==NULL)			return FALSE;
	if(u==NULL)			return FALSE;
	if(U[0]->Nel!=E[0]->Nel)	return FALSE;
	if(U[0]->Nel!=u->Nel)		return FALSE;

	for(i=0;i<m;i++)
	{
		d=pds_bvector_xor(U[i],E[i],u);
		if(d==FALSE)	return FALSE;
	}

	return TRUE;
}

int pds_bvectors_media(PdsBVector *_u_,PdsBVector **U,PdsBaNatural m)
{
	PdsBaNatural i,j;
	PdsBaNatural S;
	PdsBaBit b,x;
	int d;

	if(U==NULL)			return FALSE;
	if(_u_==NULL)			return FALSE;
	if(U[0]->Nel!=_u_->Nel)		return FALSE;

	
	for(i=0;i<U[0]->Nel;i++)
	{
		S=0;
		for(j=0;j<m;j++)
		{
			d=pds_bvector_get_bit(U[j],i,&b);
			if(d==FALSE)	return FALSE;
			S=S+b;
		}

		if(S>(m/2.0))			pds_bvector_set_bit(_u_,i,1);
		else if(S<(m/2.0))		pds_bvector_set_bit(_u_,i,0);
		else if((S*1.0)==(m/2.0))	
		{
			x=rand()%2;
			pds_bvector_set_bit(_u_,i,x);
		}
	}

	return TRUE;
}

//probabilidad media de 1 en los vectores
int pds_bvectors_prob(PdsBaReal *p,PdsBVector **U,PdsBaNatural m)
{
	PdsBaNatural i;
	PdsBaNatural S,n;

	if(U==NULL)			return FALSE;
	if(p==NULL)			return FALSE;
	
	for(i=0,S=0;i<m;i++)
	{
		pds_bvector_get_ones(U[i],&n);
		S=S+n;
	}

	*p=(1.0*S)/(m*U[0]->Nel);

	return TRUE;
}

int pds_bvectors_printf (PdsBVector **U,PdsBaNatural m)
{
	PdsBaNatural i;

	if(U==NULL)	return FALSE;

	for(i=0;i<m;i++)
	{
		printf("u[%3d]: ",i);	pds_bvector_printf(U[i]);
	}
		
	return TRUE;
}


void pds_bvectors_free(PdsBVector **X,PdsRaNatural m)
{
	PdsRaNatural i;
	
	if(X!=NULL)	for(i=0;i<m;i++)	pds_bvector_free(X[i]);
}



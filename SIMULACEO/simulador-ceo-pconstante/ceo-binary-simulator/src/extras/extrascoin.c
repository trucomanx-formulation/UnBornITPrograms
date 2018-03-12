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
#include "extrascoin.h"

#include <pds/pdsrv.h>
#include <pds/pdsra.h>
#include <pds/pdsba.h>


int pds_coin_get_bvector(PdsCoin *X, PdsBVector *V)
{
	PdsRvNaturalD xn;
	PdsRaNatural i;

	if(X==NULL)		return FALSE;
	if(V==NULL)		return FALSE;

	for(i=0;i<V->Nel;i++)
	{
		pds_congruential_get_value(X->X1,&xn);
	
		if( (xn*1.0/PDS_RAND_MAX)<X->p )	
		{
			X->x=+1;
			pds_bvector_set_bit(V,i,1);
		}
		else
		{
			X->x=-1;
			pds_bvector_set_bit(V,i,0);
		}

	}

	return TRUE;
}


int pds_bvector_load_coin(PdsBVector *u,PdsCoin *c)
{
	int d;
	
	d=pds_coin_get_bvector(c,u);
	if(d==FALSE)	return FALSE;

	return TRUE;
}

PdsCoin **pds_coins_new(PdsVector *p)
{
	PdsCoin **C=NULL;
	PdsRaNatural i,j;
	
	C=(PdsCoin **)calloc(p->Nel,sizeof(PdsCoin *));
	if(C==NULL)	return NULL;

	for(i=0;i<p->Nel;i++)
	{
		C[i]=pds_coin_new(p->V[i]);
		if(C[i]==NULL)
		{
			for(j=0;j<i;j++) pds_coin_free(C[j]);
			return NULL;
		}
	}

	return C;
}

int pds_coins_set_p(PdsCoin **X, PdsVector *P)
{

	PdsRaNatural i;

	if(X==NULL)	return FALSE;
	if(P==NULL)	return FALSE;

	for(i=0;i<P->Nel;i++)	X[i]->p=P->V[i];

	return TRUE;
}

void pds_coins_free(PdsCoin **X,PdsRaNatural m)
{
	PdsRaNatural i;
	
	if(X!=NULL)	for(i=0;i<m;i++)	pds_coin_free(X[i]);
}


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
#include "extras.h"

#include <pds/pdsrv.h>
#include <pds/pdsra.h>
#include <pds/pdsba.h>



#include <math.h>

#include <pds/pdsdatafunc.h>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/** \fn int pds_bvector_get_ones(const PdsBVector *BVector,PdsBaNatural *n)
 *  \brief Encuentra el número de unos en el vector BVector.
 *  \param[in] BVector Es el vector binarioa a consultar.
 *  \param[out] n Es el número unos en el vector BVector.
 *  \return TRUE si todo fue bien, o FALSE sino. Ejem: BVector==NULL.
 *  \ingroup PdsBVectorGroup
 */
int pds_bvector_get_ones(const PdsBVector *BVector,PdsBaNatural *n)
{
	PdsBaNatural i,j,k;
	PdsBaByte Byte;

	*n=0;

	if(BVector==NULL)			return FALSE;

	for(i=0,k=0;i<BVector->Nbytes;i++)
	{
		Byte=BVector->V[i];
		for(j=0;(j<8)&&(k<BVector->Nel);j++,k++)
		{
			*n=*n+((Byte&(1<<j))>>j);
		}
	}
	return TRUE;
}

int pds_bvector_xor(PdsBVector *BVector1,const PdsBVector *BVector2,const PdsBVector *BVector3)
{
	PdsBaNatural i;

	if(BVector1==NULL)			return FALSE;
	if(BVector2==NULL)			return FALSE;
	if(BVector3==NULL)			return FALSE;
	if(BVector1->Nel!=BVector2->Nel)	return FALSE;	
	if(BVector2->Nel!=BVector3->Nel)	return FALSE;	
	if(BVector1->Nbytes!=BVector2->Nbytes)	return FALSE;	
	if(BVector2->Nbytes!=BVector3->Nbytes)	return FALSE;	

	for(i=0;i<BVector1->Nbytes;i++)
	{
		BVector1->V[i]=(BVector2->V[i])^(BVector3->V[i]);
	}
	return TRUE;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int pds_bvector_mean_bvector (const PdsBVector *VectorX,PdsBaReal *mx)
{
	PdsBaBit m;
	PdsBaReal s;
	PdsBaNatural i;

	if(VectorX==NULL)	return FALSE;

	for(i=0,s=0;i<VectorX->Nel;i++)
	{
		pds_bvector_get_bit (VectorX, i, &m);
		s=s+m;
	}
	*mx=s/VectorX->Nel;

	return TRUE;
}

/** \fn int pds_bvector_cor_bvector (const PdsBVector *VectorX,const PdsBVector *VectorY, PdsBaReal *c)
 *  \brief Devuelve el valor del coeficiente de correlación muestral de los vectores VectorX y VectorY.
 *
 *  \f[ \eta_x=\frac{\sum_{i=0}^{Nel-1} {X_i}}{Nel} \f]
 *  \f[ \eta_y=\frac{\sum_{i=0}^{Nel-1} {Y_i}}{Nel} \f]
 *  \f[ cor(X,Y)=\frac{\sum_{i=0}^{Nel-1} {(X_i -{\eta}_x)(Y_i -{\eta}_y)}}{\sqrt{\sum_{i=0}^{Nel-1} {(X_i -{\eta}_x)^2}}\sqrt{\sum_{i=0}^{Nel-1} {(Y_i -{\eta}_y)^2}}} \f]
 *  \param[in] VectorX El vector en consulta.
 *  \param[in] VectorY El vector en consulta.
 *  \param[out] c El valor de la correlación de los vectores VectorX y VectorY.
 *  \return TRUE si todo fue bien o FALSE si no (ej: VectorX==NULL, VectorY==NULL o
 *  longitudes distintas). 
 *  \ingroup PdsVectorGroup
 */
int pds_bvector_cor_bvector (const PdsBVector *VectorX,const PdsBVector *VectorY, PdsBaReal *c)
{
	PdsBaReal mx,my,txy,tx,ty;
	PdsBaBit a,b;
	PdsBaNatural i,id;

	*c=0;

	if(VectorX==NULL)		return FALSE;
	if(VectorY==NULL)		return FALSE;
	if(VectorX->Nel!=VectorY->Nel)	return FALSE;

	id=pds_bvector_mean_bvector (VectorX,&mx);
	if(id==FALSE)			return FALSE;
	id=pds_bvector_mean_bvector (VectorY,&my);
	if(id==FALSE)			return FALSE;

	for(txy=0,tx=0,ty=0,i=0;i<VectorX->Nel;i++)
	{
		pds_bvector_get_bit (VectorX, i, &a);
		pds_bvector_get_bit (VectorY, i, &b);

		txy=txy+(a-mx)*(b-my);
		tx =tx +(a-mx)*(a-mx);
		ty =ty +(b-my)*(b-my);
	}

	(*c)=txy/sqrt(tx*ty);
	return TRUE;
}

/*
 * extras.h
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

#ifndef _EXTRAS_H_
#define _EXTRAS_H_

#include <pds/pdsrv.h>
#include <pds/pdsra.h>
#include <pds/pdsba.h>


////////////////////////////////////////////////////////////////////////////////
// pds/pdsba.h
////////////////////////////////////////////////////////////////////////////////
// BVector1 = BVector2 XOR BVector3
int pds_bvector_xor(PdsBVector *BVector1,const PdsBVector *BVector2,const PdsBVector *BVector3);
/** \fn int pds_bvector_get_ones(const PdsBVector *BVector,PdsBaNatural *n)
 *  \brief Encuentra el número de unos en el vector BVector.
 *  \param[in] BVector Es el vector binarioa a consultar.
 *  \param[out] n Es el número unos en el vector BVector.
 *  \return TRUE si todo fue bien, o FALSE sino. Ejem: BVector==NULL.
 *  \ingroup PdsBVectorGroup
 */
int pds_bvector_get_ones(const PdsBVector *BVector,PdsBaNatural *n);


////////////////////////////////////////////////////////////////////////////////
// pds/pdsbm.h
////////////////////////////////////////////////////////////////////////////////

int pds_bvector_mean_bvector (const PdsBVector *VectorX,PdsBaReal *mx);
int pds_bvector_cor_bvector (const PdsBVector *VectorX,const PdsBVector *VectorY, PdsBaReal *c);
#endif /* _EXTRAS_H_ */

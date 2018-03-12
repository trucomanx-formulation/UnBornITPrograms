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

#include <pds/pdsba.h>
#include <pds/pdsra.h>


/** \fn int pds_fast_joint_probability_bsc_model(const PdsVector *Ps, double *prob)
 *  \brief Encuentra la entropia conjunta para M fontes generadas
 *  pasando una fuente U0, con probabilidade Pr(U0=1)=0.5, atraves de M canales BSC
 *  con probabilidades de error Ps=[Ps0, Ps1, ..., Ps(M-1)].
 *  \param[in] Ps Es el vector de probabilidades de error de las fuentes BSC.
 *  \param[out] H Entropia conjuta.
 *  \return TRUE si todo fue bien, o FALSE sino. Ejem: Ps==NULL.
 *  \ingroup PdsBVectorGroup
 */
int pds_fast_joint_entropy_bsc_model(const PdsVector *Ps, double *H);


int pds_fast2_joint_entropy_bsc_model(const PdsVector *Ps, double *H);

int pds_ordenap(PdsVector *Ps);
int pds_ordenap2(PdsVector *Ps);
int pds_ordenap2a(PdsVector *Ps);
int pds_ordenap_hu0rho2(PdsVector *Ps);
int pds_ordenap_hu0rho4(PdsVector *Ps);

#endif /* _EXTRAS_H_ */

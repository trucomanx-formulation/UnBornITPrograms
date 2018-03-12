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

#ifndef __EXTRASBVECTOR_H__
#define __EXTRASBVECTOR_H__

#include <pds/pdsrv.h>
#include <pds/pdsra.h>
#include <pds/pdsba.h>

PdsBVector **pds_bvectors_new(PdsBaNatural m,PdsBaNatural k);

int pds_bvectors_printf (PdsBVector **U,PdsBaNatural m);

int pds_bvectors_xor(PdsBVector **U,PdsBVector **E,PdsBaNatural m,PdsBVector *u);

int pds_bvectors_load_coins(PdsBVector **u,PdsCoin **c,PdsBaNatural m);

void pds_bvectors_free(PdsBVector **X,PdsRaNatural m);

int pds_bvectors_media(PdsBVector *u1,PdsBVector **U,PdsBaNatural m);

//probabilidad media de 1 en los vectores
int pds_bvectors_prob(PdsBaReal *p,PdsBVector **U,PdsBaNatural m);

#endif /* __EXTRASBVECTOR_H__ */

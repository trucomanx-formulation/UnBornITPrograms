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

#ifndef __EXTRASCOIN_H__
#define __EXTRASCOIN_H__

#include <pds/pdsrv.h>
#include <pds/pdsra.h>
#include <pds/pdsba.h>


int pds_coin_get_bvector(PdsCoin *X, PdsBVector *V);

// Carga resultado de iterar una moneda en cada elemento del vector.
int pds_bvector_load_coin(PdsBVector *u,PdsCoin *c);

PdsCoin **pds_coins_new(PdsVector *p);

int pds_coins_set_p(PdsCoin **X, PdsVector *P);

void pds_coins_free(PdsCoin **X,PdsRaNatural m);


#endif /* __EXTRASCOIN_H__ */

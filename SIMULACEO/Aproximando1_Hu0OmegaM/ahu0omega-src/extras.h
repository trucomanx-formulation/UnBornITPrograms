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


int pds_ordenap2(PdsVector *Ps);

int pds_ordenap_hu0rho2(PdsVector *Ps);

int pds_ordenap_hu0rho4(PdsVector *Ps);

int pds_ordenap_hu0rhoN(PdsVector *Ps,short int halfN);

#endif /* _EXTRAS_H_ */

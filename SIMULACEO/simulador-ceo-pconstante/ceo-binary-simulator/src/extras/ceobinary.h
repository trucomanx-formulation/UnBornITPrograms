/*
 * ceobinary.h
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

#ifndef __CEOBINARY_H__
#define __CEOBINARY_H__

#include <stdlib.h>
#include <stdio.h>
#include <pds/pdsra.h>
#include <pds/pdsba.h>

void ceobinary_help(void);

void ceobinary_load_data(int argc,char *argv[],char **file,char **filelog,int *m,float *beta,float *PeMINIMO);

void ceobinary_data_fprintf(FILE *fd,char *file,char *filelog,int m,float beta,float PeMINIMO);

void ceobinary_itera_fprintf(FILE *fd,PdsRaReal Ps,PdsRaReal Pe,PdsRaNatural k,PdsRaNatural ITER);

void ceobinary_next_prob(PdsRaReal p_now,PdsRaReal *p_next,PdsRaReal p_min,PdsRaReal beta);

void carga_correlacion(PdsMatrix *Cor,PdsBVector **U,int M);

#endif /* __CEOBINARY_H__ */

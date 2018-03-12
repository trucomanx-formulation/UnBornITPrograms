/*
 * main.c
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
#include <pds/pdsrv.h>
#include <pds/pdsra.h>
#include <pds/pdsba.h>

#include <math.h>

#define	DEFAULT_K		1000
#define	DEFAULT_ERRORESMAX	1000
#define	DEFAULT_M		10
#define	DEFAULT_BETA		1.0
#define	DEFAULT_PEMINIMO	1.0e-6
#define	DEFAULT_PS_INIT		0.5

#include "extras/extrascoin.h"
#include "extras/extrasbvector.h"
#include "extras/extras.h"

#include "extras/ceobinary.h"

int main(int argc, char** argv)
{
	char  *file=NULL;	FILE *fd=NULL;
	char  *filelog=NULL;	FILE *fdlog=NULL;
	int   m=DEFAULT_M;
	float beta=DEFAULT_BETA;
	float PeMINIMO=DEFAULT_PEMINIMO;
	PdsRaReal ps=DEFAULT_PS_INIT;
	PdsBaReal pstmp;
	double psmedio;

	PdsRaNatural k=DEFAULT_K;
	PdsBaNatural n;
	unsigned long int n_unos_u;
	unsigned long int Error;
	unsigned long int ERRORESMAX;

	PdsRaNatural ITER,i;
	PdsRaReal    Pe;

	PdsBVector **U=NULL;

	PdsBVector **E=NULL;
	PdsVector  *Ps=NULL;
	PdsCoin   **C=NULL;

	PdsBVector *_u_=NULL;
	PdsBVector *u=NULL;
	PdsCoin   *c=NULL;

	PdsMatrix *Cor=NULL;
	PdsMatrix *CorS=NULL;

	ceobinary_load_data(argc,argv,&file,&filelog,&m,&beta,&PeMINIMO);

	Cor =pds_matrix_new(m,m);
	CorS=pds_matrix_new(m,m);
	
	ERRORESMAX=(unsigned long int)(DEFAULT_ERRORESMAX*beta);

	fdlog=fopen(filelog,"w");
	ceobinary_data_fprintf(stdout,file,filelog,m,beta,PeMINIMO);
	ceobinary_data_fprintf(fdlog ,file,filelog,m,beta,PeMINIMO);
	fclose(fdlog);

	remove(file);	remove(filelog);
	
	// fontes ui
	U=pds_bvectors_new(m,k);

	// fontes ei
	E=pds_bvectors_new(m,k);		
	Ps=pds_vector_new(m);
	C=pds_coins_new(Ps);
	
	// fonte u
	u=pds_bvector_new(k);
	c=pds_coin_new(0.5);
	
	// vector fonte estimado
	_u_=pds_bvector_new(k);

	Pe=1.0;

	while( (ps>0)&&(Pe>PeMINIMO) )
	{
		//Creando las monedas de los [E1, E2, ...,Em]
		pds_vector_init_value(Ps,ps);
		pds_coins_set_p(C,Ps);

		//Erro: Bits errados al enviar los ITER*k bits
		Error=0;
		n_unos_u=0;
		psmedio=0;
		pds_matrix_init_value (Cor,0);
		for(i=0;Error<ERRORESMAX;i++)
		{
			//Juego una moneda c para cargar el vector u
			pds_bvector_load_coin(u,c);
			pds_bvector_get_ones(u,&n);
			n_unos_u=n_unos_u+n;

			//Juego C monedas para cargar los m vectores Ej
			pds_bvectors_load_coins(E,C,m);

			//Probabilidad media de 1 en los vectores Ej
			pds_bvectors_prob(&pstmp,E,m);
			psmedio=psmedio+pstmp;

			//Uj=Ej XOR u
			pds_bvectors_xor(U,E,m,u);

			///////////////////////////////
			carga_correlacion(CorS,U,m);
			pds_matrix_add_matrix(Cor,CorS);
			///////////////////////////////

			// Obtengo el valor estimado de u a partir de los vectores Uj
			pds_bvectors_media(_u_,U,m);

			// Compara el vector u con el vector estimado _u_
			pds_bvector_cmp(u,_u_,&n);
	
			Error=Error+n;	
		}
		ITER=i;

		printf("(1-2 rho)^2 : %f\n",(1.0-2.0*ps)*(1.0-2.0*ps));
		pds_matrix_mul_value (Cor,1.0/ITER);
		printf("Cor:\n");pds_matrix_printf(Cor);

		psmedio=psmedio/ITER;
		printf("P(u=1)=%f\t",(1.0*n_unos_u)/(k*ITER));

		// Calula Pe 
		Pe=(1.0*Error)/(k*ITER);


		// Escribe registros de iteracion en el stdout
		ceobinary_itera_fprintf(stdout,psmedio,Pe,k,ITER);

		// Escribe registros de iteracion en el log
		fdlog=fopen(filelog,"a+");
		ceobinary_itera_fprintf(fdlog ,psmedio,Pe,k,ITER);
		fclose(fdlog);

		// Escribe los datos de la simulación.
		fd=fopen(file,"a+");
		fprintf(fd,"%e\t%e\n",psmedio,Pe);
		fclose(fd);

		ps=ps-0.0500000000;
	}

	pds_bvectors_free(U,m);
	pds_bvector_free(_u_);
	
	pds_bvectors_free(E,m);
	pds_vector_free(Ps);
	pds_coins_free(C,m);

	pds_bvector_free(u);
	pds_coin_free(c);
	return 0;
}





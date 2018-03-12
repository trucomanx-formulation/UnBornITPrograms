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
#include <pds/pdsdatafunc.h>

#include "extras/extrascoin.h"
#include "extras/extrasbvector.h"
#include "extras/extras.h"
#include "extras/ceobinary.h"

#include <math.h>

#define	DEFAULT_M		2
#define	DEFAULT_PS_INIT		0.5
#define	DEFAULT_PU0_INIT	0.5

float InfMutuaExperimental(float Px0,float Px1,float Py0,float Py1,float P00,float P01,float P10,float P11);
void mutua_load_data(int argc,char *argv[],char **file,float *pu0,int *m);
void mutua_help(void);

int main(int argc, char** argv)
{
	char  *file=NULL;	FILE *fd=NULL;
	int   m=DEFAULT_M;


	PdsRaNatural k=4096;
	PdsBaNatural n,ntotal;

	PdsRaNatural ITER,i,j;

	unsigned long int  Lx0,Lx1,Ly0,Ly1,L00,L01,L10,L11;
	float Px0,Px1,Py0,Py1,P00,P01,P10,P11;
	float Inf,Ps_mix,Inf_teorico,pu0=DEFAULT_PU0_INIT,Pu0;
	PdsBaBit dx,dy;

	PdsBVector **U=NULL;

	PdsBVector **E=NULL;
	PdsVector  *Ps=NULL;
	PdsCoin   **C=NULL;

	PdsBVector *_u_=NULL;
	PdsBVector *u=NULL;
	PdsCoin   *c=NULL;

	PdsUniform *RV_Ps=NULL;

	mutua_load_data(argc,argv,&file,&pu0,&m);


	remove(file);	
	
	// fontes ui
	U=pds_bvectors_new(m,k);

	// fontes ei
	RV_Ps=pds_uniform_new (0.05, 0.95);
	E=pds_bvectors_new(m,k);		
	Ps=pds_vector_new(m);
	C=pds_coins_new(Ps);
	
	// fonte u
	u=pds_bvector_new(k);
	c=pds_coin_new(pu0);
	
	// vector fonte estimado
	_u_=pds_bvector_new(k);

	//Creando las monedas de los [E1, E2, ...,Em]
	pds_uniform_get_vector (RV_Ps, Ps);
	pds_coins_set_p(C,Ps);

	// ITER: Numero de paquetes enviados de tamaño k
	ITER=4096;
	if(ITER<=8)	ITER=8;

	//Erro: Bits errados al enviar los ITER*k bits
	ntotal=0;
	Lx0=0;	Lx1=0;	Ly0=0;	Ly1=0;
	L00=0;	L01=0;	L10=0;	L11=0;
	for(i=0;i<ITER;i++)
	{
		//Juego una moneda c para cargar el vector u
		pds_bvector_load_coin(u,c);
		pds_bvector_get_ones(u,&n);
		ntotal=ntotal+n;//ntotal:numero de unos en el vector u.

		//Juego C monedas para cargar los m vectores Ej
		pds_bvectors_load_coins(E,C,m);

		//Uj=Ej XOR u
		pds_bvectors_xor(U,E,m,u);

		// cargo la cantidad de 1 y 0 en la 1ra fuente.
		pds_bvector_get_ones(U[0],&n);
		Lx1=Lx1+n;
		Lx0=Lx0+(k-n);

		// cargo la cantidad de 1 y 0 en la 2da fuente.
		pds_bvector_get_ones(U[1],&n);
		Ly1=Ly1+n;
		Ly0=Ly0+(k-n);

		// cargo la cantidad de 00,01,10,11 en las fuentes.
		for(j=0;j<k;j++)
		{
			pds_bvector_get_bit (U[0],j, &dx);
			pds_bvector_get_bit (U[1],j, &dy);

			if((dx==0)&&(dy==0))	//0,0
			{	L00=L00+1;	}

			if((dx==0)&&(dy==1))	//0,1
			{	L01=L01+1;	}

			if((dx==1)&&(dy==0))	//1,0
			{	L10=L10+1;	}

			if((dx==1)&&(dy==1))	//1,1
			{	L11=L11+1;	}
		}

	}

	Pu0=(ntotal*1.0)/(k*ITER);

	Px0=(Lx0*1.0)/(k*ITER);	Px1=(Lx1*1.0)/(k*ITER);
	Py0=(Ly0*1.0)/(k*ITER);	Py1=(Ly1*1.0)/(k*ITER);

	P00=(L00*1.0)/(k*ITER);
	P01=(L01*1.0)/(k*ITER);
	P10=(L10*1.0)/(k*ITER);
	P11=(L11*1.0)/(k*ITER);

	Inf=InfMutuaExperimental(Px0,Px1,Py0,Py1,P00,P01,P10,P11);

	printf("\n");
	printf("m=%d\n",m);
	printf("Pu0[1]=%f\n",Pu0);
	printf("k=%d\n",k);
	printf("Bloques=%d\n",ITER);
	printf("L=Bloques*k=%d\t=Bits enviados por fuente\n",ITER*k);
	printf("\n");
	printf("Ps1=%f\n",Ps->V[0]);
	printf("Ps2=%f\n",Ps->V[1]);
	printf("\n");
	printf("Pu1[0]=%f\tPu1[1]=%f\n",Px0,Px1);
	printf("Pu2[0]=%f\tPu2[1]=%f\n",Py0,Py1);		
	printf("\n");
	printf("P00=%f\tP01=%f\n",P00,P01);
	printf("P10=%f\tP11=%f\n",P10,P11);
	printf("\n");
	Ps_mix=Ps->V[0]+Ps->V[1]-2.0*Ps->V[0]*Ps->V[1];
	Inf_teorico=1+(Ps_mix*log2(Ps_mix)+(1.0-Ps_mix)*log2(1.0-Ps_mix));
	printf("Inf. Mutua Exp.=%f\n",Inf);
	printf("Inf. Mutua Teo.=%f\t=1-h(Ps1+Ps2-2*Ps1*Ps2)\n",Inf_teorico);
	printf("\n");



	// Escribe los datos de la simulación.
	fd=fopen(file,"a+");
	fprintf(fd,"Inf. Mutua Exp.=%f\n",Inf);
	fprintf(fd,"Inf. Mutua Teorica.=%f\n",Inf_teorico);
	fclose(fd);


	pds_bvectors_free(U,m);
	pds_bvector_free(_u_);
	
	pds_bvectors_free(E,m);
	pds_vector_free(Ps);
	pds_coins_free(C,m);

	pds_bvector_free(u);
	pds_coin_free(c);
	return 0;
}

float InfMutuaExperimental(float Px0,float Px1,float Py0,float Py1,float P00,float P01,float P10,float P11)
{
	float I00,I01,I10,I11;
	float I;

	I00=log2(P00/(Px0*Py0));
	I01=log2(P01/(Px0*Py1));
	I10=log2(P10/(Px1*Py0));
	I11=log2(P11/(Px1*Py1));

	I= P00*I00+P01*I01+P10*I10+P11*I11;
	return I;
}

void mutua_load_data(int argc,char *argv[],char **file,float *pu0,int *m)
{
	int id;

	if(argc==1)	ceobinary_help();

	//Muestra la ayuda
	id=pds_exist_param(argc,argv,"-h");	if(id==TRUE)	mutua_help();
	id=pds_exist_param(argc,argv,"--help");	if(id==TRUE)	mutua_help();
	//Numero de fuentes
	id=pds_get_int_param(argc,argv,"-m",m);

	//Probabilidade de la fuente original u0
	id=pds_get_float_param(argc,argv,"-p",pu0);
	id=pds_get_float_param(argc,argv,"--pu0",pu0);

	if(*m!=2)
	{
		printf("\"m\" solo puede ser 2. Este programa es de testeo de inf. mutua.\n");
		*m=2;
	}

	//Nombre del archivo de salida de la simulación.
	id=pds_get_chars_param(argc,argv,"-f",file);
	id=pds_get_chars_param(argc,argv,"--file",file);

	//Defino el nombre del archivo de salida
	if(*file==NULL)	
	{	
		*file=(char *)calloc(128,sizeof(char));
		sprintf(*file,"data%d.txt",*m);
	}

}

void mutua_help(void)
{
	printf("\n");
	printf("Ui=U XOR Ei\t\ti ={1, 2}\n\n");
	printf("\t\t\tP(u0=1) =p\n");
	printf("\t\t\tP(ei=1)=pi\n");
	printf("\t\t\tP(ui=1)=p+pi-2*p*pi\n\n");
	printf("El programa muestra la informacion mutua ente u1 y u2.\n");
	printf("\n");
	printf("mutua [Opciones]\n");
	printf("\n");
	printf("Opciones:\n");
	printf("\t--help\n");
	printf("\t-h            \tMuestra esta ayuda y termina el programa.\n");
	printf("\n");
	printf("\t-m      numero\tNúmero de fuentes, siempre es m=2\n");
	printf("\n");
	printf("\t--pu0   prob  \n");
	printf("\t-p      prob  \tProbabilidad de Pr[u0=1]. \n");
	printf("\t              \tPor defecto Pr[u0=1]=0.5\n");
	printf("\n");
	printf("\t--file  nombre\n");
	printf("\t-f      nombre\tNombre del archivo de salida de probabilidades. \n");
	printf("\t              \tEstá ordenado en dos columnas la primera de probabilidades\n");
	printf("\t              \tde P(ei=1) y la segunda de probabilidades de error al\n");
	printf("\t              \testimar U. Por defecto nombre=data{m}.txt\n");
	printf("\n");
	exit(0);
}


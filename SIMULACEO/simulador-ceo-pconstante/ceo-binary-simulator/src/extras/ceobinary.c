/*
 * ceobinary.c
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
#include "ceobinary.h"
#include <pds/pdsdatafunc.h>
#include <pds/pdsra.h>
#include <pds/pdsba.h>
#include "extras.h"

void ceobinary_help(void)
{
	printf("\n");
	printf("U_i=U XOR E_i\t\ti ={1, 2, ..., m}\n\n");
	printf("\t\t\tP(u=1)  =p   = 0.5\n");
	printf("\t\t\tP(e_i=1)=p_i = rho\n");
	printf("\t\t\tP(u_i=1)=p+p_i-2*p*p_i\n\n");
	printf("El programa estima U a partir de U_i, i={0, 1, ..., M}.\n");
	printf("Retorna un archivo de texto con dos columnas, la primer con\n");
	printf("valores de probabilidad rho y la segunda com la probabilidade\n");
	printf("de error de la estimacion.\n");
	printf("\n");
	printf("ceo-binary-simulator [Opciones]\n");
	printf("\n");
	printf("Opciones:\n");
	printf("\t--help\n");
	printf("\t-h            \tMuestra esta ayuda y termina el programa.\n");
	printf("\n");
	printf("\t-m      numero\tNúmero de fuentes. Por defecto m=10\n");
	printf("\n");
	printf("\t--beta  factor\n");
	printf("\t-b      factor\tFactor de incremento del número de errores.\n");
	printf("\t              \tPor defecto factor=1. Por defecto \n");
	printf("\t              \tERRORESMAX=factor*1000.\n");//,DEFAULT_ERRORESMAX);
	printf("\t              \tLa simulacion se detiene despues de ERRORESMAX\n");
	printf("\t              \terrores encontrados.\n");
	printf("\n");
	printf("\t--pminimo pmin\n");
	printf("\t-p        pmin\tProbabilidad mínima confiable del programa, mientras\n");
	printf("\t              \tmas pequeño mas confiable pero mas lento. Por defecto\n");
	printf("\t              \tpmin=1.0e-6\n");
	printf("\n");
	printf("\t--file  nombre\n");
	printf("\t-f      nombre\tNombre del archivo de salida de probabilidades. \n");
	printf("\t              \tEstá ordenado en dos columnas la primera de probabilidades\n");
	printf("\t              \tde P(ei=1) y la segunda de probabilidades de error al\n");
	printf("\t              \testimar U. Por defecto nombre=data{m}.txt\n");
	printf("\n");
	exit(0);
}


void ceobinary_load_data(int argc,char *argv[],char **file,char **filelog,int *m,float *beta,float *PeMINIMO)
{
	int id;

	if(argc==1)	ceobinary_help();

	//Muestra la ayuda
	id=pds_exist_param(argc,argv,"-h");	if(id==TRUE)	ceobinary_help();
	id=pds_exist_param(argc,argv,"--help");	if(id==TRUE)	ceobinary_help();
	//Numero de fuentes
	id=pds_get_int_param(argc,argv,"-m",m);
	//Factor de prediccion de paquetes enviados en el simulador
	id=pds_get_float_param(argc,argv,"-b",beta);	
	id=pds_get_float_param(argc,argv,"--beta",beta);	
	//Probabilidad minima simulada, despues se detiene.
	id=pds_get_float_param(argc,argv,"-p",PeMINIMO);	
	id=pds_get_float_param(argc,argv,"--pminimo",PeMINIMO);	
	//Nombre del archivo de salida de la simulación.
	id=pds_get_chars_param(argc,argv,"-f",file);
	id=pds_get_chars_param(argc,argv,"--file",file);

	//Defino el nombre del archivo de salida
	if(*file==NULL)	
	{	
		*file=(char *)calloc(128,sizeof(char));
		sprintf(*file,"data%d.txt",*m);
	}

	//Defino el nombre del archivo de registro de actividad.
	*filelog=(char *)calloc(strlen(*file)+8,sizeof(char));
	sprintf(*filelog,"%s.log",*file);

}


void ceobinary_data_fprintf(FILE *fd,char *file,char *filelog,int m,float beta,float PeMINIMO)
{
	if(fd!=NULL)
	{
		fprintf(fd,"m=%d\n",m);
		fprintf(fd,"beta=%f\n",beta);
		fprintf(fd,"Pe(minimo)=%e\n",PeMINIMO);
		fprintf(fd,"file=%s\n",file);
		fprintf(fd,"filelog=%s\n",filelog);
	}
}

void ceobinary_itera_fprintf(FILE *fd,PdsRaReal Ps,PdsRaReal Pe,PdsRaNatural k,PdsRaNatural ITER)
{
	if(fd!=NULL)
	{
		fprintf(fd,"Ps:%f\tPe:%e\tITERrecom:%d\tITER:%d\n"
							,Ps
							,Pe
							,(int)(1000.0/((Pe+1.0e-12)*k))
							,(int)ITER);
	}
}

void ceobinary_next_prob(PdsRaReal p_now,PdsRaReal *p_next,PdsRaReal p_min,PdsRaReal beta)
{
	static PdsRaReal p[3]={0.0,0.0,0.0};
	int N=3;
	int i;

	//verifico que no tenga datos cero
	if(p_now<=p_min)	p_now=p_min;

	// coloco el nuevo dato p_now en p[N-1] y disloco los datos anteriores.
	for(i=0;i<(N-1);i++)	p[i]=p[i+1];
	p[N-1]=p_now;

	// Predigo el p siguiente.
	if(p[0]<=0.0)	*p_next=p_now/10;
	else
	{
		*p_next=(p[0]/beta)*(p[2]/p[1])*(p[2]/p[1])*(p[2]/p[1]);
	}
}


void carga_correlacion(PdsMatrix *Cor,PdsBVector **U,int M)
{
	int i,j;
	PdsBaReal c;


	for(i=0;i<M;i++)
	for(j=i;j<M;j++)
	{
		c=0;
		pds_bvector_cor_bvector (U[i],U[j], &c);
		Cor->M[i][j]=c;
		Cor->M[j][i]=c;
	}
}

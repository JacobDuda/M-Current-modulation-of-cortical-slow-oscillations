#ifndef Network_H
#define Network_H

#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>

// #include "Randomnumber.h"
// #include "Variables.h"

using namespace std;

long sed=1234567;

double ran1(long *idum) {

	int j;
	long k;
	static long iy=0;
	static long iv[NTAB1];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum = 1;
		else *idum = -(*idum);
		for (j=NTAB1+7;j>=0;j--) {
			k = (*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-IR1*k;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB1) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-IR1*k;
	if (*idum < 0) *idum += IM1;
	j = iy/NDIV1;
	iy = iv[j];
	iv[j] = *idum;
	if ((temp=AM1*iy) > RNMX1) return RNMX1;
	else return temp;
}

double gaussian(long *idum) {

        static int iset=0;
        static double gset;
        double fac, rsq, v1, v2;

        if (iset == 0) {
                do {
                        v1 = 2.0*ran1(idum)-1.0;
                        v2 = 2.0*ran1(idum)-1.0;
                        rsq=v1*v1+v2*v2;
                } while (rsq>=1.0 || rsq == 0);

                fac = sqrt(-2.0*log(rsq)/rsq);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset = 0;
                return gset;
        }
}


/* routine that implements cyclic boundary conditions */
int check(int i, int n){
	int j = i;
	while (j < 0)
		j += n;
	while (j >= n)
		j -= n;
	return j;
}


void setConnectivity(int (*mat)[100], int numPre, int numPost, int conn,
			    double sigma, double jplus, int type_of_conn){
	int nConnx, postN, j;
	double frac = (double) numPost/numPre;

	for (int i=0;i<numPre;i++){
		nConnx = (int) rint(conn*(0.25*gaussian(&sed)+1.));
		// nConnx = (int) rint(conn*(0.25*1.000523+1.));

		if (nConnx>98) nConnx=98;
		if (nConnx<0) nConnx=0;

		int actual_j=0;

		for (j=0; j<nConnx; j++)
		{
			postN = (int) rint(double(i*frac + numPost*sigma*gaussian(&sed)));
			// postN = (int) rint(double(i*frac + numPost*1.000523));
			if ((numPre==numPost)&&(postN==i))
				j=j-1; //do not allow for autapses
			else if ((postN>=0)&&(postN<numPost)) //remove the if-clause if ring model
			{
				mat[i][actual_j] = check(postN,numPost);	//ToDo: hace los chequeos la guarda, con lo cual la llamada a check es al dope, deberia ser postN
				actual_j++;
			}
		}
		mat[i][actual_j] = numPost;	//el numOfXCell indica el fin de las conexiones para la fila i
	}
}



void construct_network(){
  FILE *outfile;

	const int numPre_EE=1024;
	const int numPre_EI=1024;
	const int numPre_IE=256;
	const int numPre_II=256;

	const int numPost_EE=1024;
	const int numPost_EI=256;
	const int numPost_IE=1024;
	const int numPost_II=256;

	double frac_EE= (double) numPost_EE/numPre_EE;
	double frac_EI= (double) numPost_EI/numPre_EE;
	double frac_II= (double) numPost_II/numPre_II;
	double frac_IE= (double) numPost_IE/numPre_IE;

	const double sigmaEE=0.05;
	const double sigmaEI=0.05;
	const double sigmaII=0.025;
	const double sigmaIE=0.025;

	const double conn=20.0; /*Number of Connec*/



	Wee = new int[numEcells][100];
	Wei = new int[numEcells][100];
	Wie = new int[numIcells][100];
	Wii = new int[numIcells][100];

	setConnectivity(Wee,numEcells,numEcells,conn,sigmaEE,10,1);
	setConnectivity(Wei,numEcells,numIcells,conn,sigmaEI,10,1);
	setConnectivity(Wie,numIcells,numEcells,conn,sigmaIE,1.62,1);
	setConnectivity(Wii,numIcells,numIcells,conn,sigmaII,1.62,1);


}






#endif

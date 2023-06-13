#ifndef Channels_ICells_H
#define Channels_ICells_H

#include<cmath>
#include "Variables.h"

/**********************************************************/
/*                Channels for Interneurons               */
/**********************************************************/
double Na_steadystate_m_I(double memV){
  double alpham=5.*0.1*(memV+35.)/(1.-exp(-0.1*(memV+35.)));
  double betam=5.*4.*exp(-0.055555555556*(memV+60.));
  return alpham/(alpham+betam);
}
double Na_inactgate_h_I(double variable,double memV){
  double alphah=5.*0.07*exp(-0.05*(memV+58.));
  double betah=5.*1./(1.+exp(-0.1*(memV+28.)));
  return (alphah*(1.0-variable))-(betah*variable);
}

double K_actgate_n_I(double variable, double memV){
  double betan=5.*0.125*exp(-0.0125*(memV+44.));
  double alphan=0.0;
  if(fabs(memV+34.)<0.5)
		alphan=5.*0.01/(0.1-0.1*0.1*0.5*(memV+34.));
	else
		alphan=5.*0.01*(memV+34.)/(1.-exp(-0.1*(memV+34.)));

  return (alphan*(1.0-variable))-(betan*variable);
}


/**********************************************************/
/*                The Dynamics start here                 */
/**********************************************************/
double syn_sigmoidI(double memV){
  return 1./(1.+exp(-0.5*(memV)));
}
void allderivsI(double *x, double *dx, double Iext, int index){
  double NaChannel=0.0;
  NaChannel=Na_steadystate_m_I(x[0])*Na_steadystate_m_I(x[0])*Na_steadystate_m_I(x[0])*x[1];
	NaChannel=gNa_I*0.2*NaChannel*(x[0]-55.0);

	double KChannel=0.0;
	KChannel=x[2]*x[2]*x[2]*x[2];
	KChannel=gK_I*0.2*KChannel*(x[0]+90.0);

	double LeakChannel=0.0;
	LeakChannel=glI[index]*(x[0]-VlI[index]);

  double inp = -LeakChannel - NaChannel - KChannel + Iext;
  dx[0] = inp/0.2;
  dx[1] = Na_inactgate_h_I(x[1],x[0]);
  dx[2] = K_actgate_n_I(x[2],x[0]);
  dx[3] = 1.*1.*syn_sigmoidI(x[0])-x[3]/10.0; /*CHECK HERE WITH COMPTE 1.*y[4]*syn_sigm*/

}
/**/
void rk4I(double curr, int ind){
	int i;
	double hh=dt*0.5;
	double h6=dt*0.166666666666666667;
	for (i=0;i<nVarsI;i++)
		xt[i]=varsI[i]+hh*dVarsI[i];
	allderivsI(xt,dxt,curr,ind);
	for (i=0;i<nVarsI;i++)
		xt[i]=varsI[i]+hh*dxt[i];
	allderivsI(xt,dxm,curr,ind);
	for (i=0;i<nVarsI;i++) {
		xt[i]=varsI[i]+dt*dxm[i];
		dxm[i] += dxt[i];
	}
	allderivsI(xt,dxt,curr,ind);
	for (i=0;i<nVarsI;i++)
		varsI[i] += h6*(dVarsI[i]+dxt[i]+2.0*dxm[i]);
}
/**/
void rungeKutta4I(double timesim, int index, double currExt){
	varsI[0]=x[0];
	varsI[1]=x[1];
	varsI[2]=x[2];
  varsI[3]=x[3];
	allderivsI(varsI,dVarsI,currExt,index);
  rk4I(currExt,index);
	x[0]=varsI[0];
	x[1]=varsI[1];
	x[2]=varsI[2];
  x[3]=varsI[3];
}
/**/





#endif

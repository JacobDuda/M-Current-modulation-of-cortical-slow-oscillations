#ifndef Channels_ECells_H
#define Channels_ECells_H

#include<cmath>
#include "Variables.h"

/***********************************************************/
/*         FAST SODIUM CHANNEL FOR E-CELLS                 */
/***********************************************************/
double Na_steadystate_m_E(double memV){
  double alpham=4.*0.1*(memV+33.)/(1.-exp(-0.1*(memV+33.)));;
  double betam=4.*4.*exp(-0.083333333333*(memV+53.7));
  return alpham/(alpham+betam);
}
double Na_inactgate_h_E(double variable,double memV){
  double alphah=4.*0.07*exp(-0.1*(memV+50.));;
  double betah=4.*1./(1.+exp(-0.1*(memV+20.)));
  return (alphah*(1.0-variable))-(betah*variable);
}
/**/
/***************************************************************/
/*          DELAYED RECTIFIER - K CHANNEL- FOR E-CELLS         */
/***************************************************************/
double K_actgate_n_E(double variable, double memV){
  double betan=4.*0.125*exp(-0.04*(memV+44.));
  double alphan=4*0.01*(memV+34)/(1.-exp(-0.1*(memV+34.)));
  return (alphan*(1.0-variable))-(betan*variable);
}
/**/
/**********************************************************/
/*                A-TYPE K CHANNEL FOR E-CELLS            */
/**********************************************************/
double A_steadystate_m_E(double memV){
  return 1./(1.+exp(-(memV+50.)/20.));
}
double A_inactgate_h_E(double variable, double memV){
  double inactivOpenRate=1./(1.+exp((memV+80.)/6.))/15.;
  double inactivCloseRate=1./15.-1./(1.+exp((memV+80.)/6.))/15.;

  return (inactivOpenRate+inactivCloseRate)*((inactivOpenRate/(inactivOpenRate+inactivCloseRate))-variable);
}
/**/
/**********************************************************/
/* SLOWLY INACTIVATING A-TYPE K CHANNEL FOR E-CELLS - KS  */
/**********************************************************/
double KS_actgate_m_E(double variable, double memV){
  double minf=1./(1.+exp(-(memV+34.)/6.5));
  double tau=8./(exp(-(memV+55.)/30.)+exp((memV+55.)/30.));
  return (minf-variable)/tau;
}
/**/
/**********************************************************/
/*   HIGH THRESHOLD CALCIUM CHANNEL FOR E-CELLS  - Ca     */
/**********************************************************/
double Ca_act_m_E(double memV){
  return 1./(1.+exp(-0.11111111*(memV+20.)));
}
/**/
/**********************************************************/
/*          ANOMALOUS RECTIFIER K CHANNEL FOR E-CELLS     */
/**********************************************************/
double AR_act_h_E(double memV){
  return 1./(1.+exp((memV+75.0)*0.25));
}
/**/
/**********************************************************/
/*          PERSISTENT NA CHANNEL FOR E-CELLS - NaP       */
/**********************************************************/
double NaP_act_m_E(double memV){
  return 1./(1.+exp(-(memV+55.7)/7.7));
}
/**/
/**********************************************************/
/*         NA DEPENDENT K  CHANNEL FOR E-CELLS * KNa       */
/**********************************************************/
double KNa_winf_E(double variable){
  double fac=pow(38.7/variable,3.5);
  return 0.37/(1.+fac);
}
/**/




double mcurrentMcCormick(double variable, double memV){
  double minf=1./(1.+exp(-0.1*(memV+35.)));
  double taum=1000.0/(3.3*(exp(0.05*(memV+35.))+exp(-0.05*(memV+35.))));
  return (minf-variable)/taum;
}

double mcurrentWang(double variable, double memV){
  double minf=1./(1.+exp(-0.1666*(memV+44.)));
  double taum=100.0/((exp(-0.083*(memV+44.))+exp(0.083*(memV+44.))));
  return (minf-variable)/taum;
}

double mcurrentHay(double variable, double memV){
  double alpham=0.0033*exp(0.1*(memV+35.));
  double betam=0.0033*exp(-0.1*(memV+35.));
  double minf=alpham/(alpham+betam);
  double tadj=pow(2.3,(34.-21.)*0.1);
  double taum=1./(tadj*(alpham+betam));

  return (minf-variable)/taum;

}

/*H-Channels Affect Frequency, Power and Amplitude Fluctuations of Neuronal Network Oscillations*/
double hchannel(double variable, double memV){
  double linf=1./(1.+exp((memV+81.0)/7.0));
  double taul=exp(0.033*(memV+75.0))/(0.02*(1.0+exp(0.083*(memV+75.0))));

  return (linf-variable)/taul;
}

double hchannelHay(double variable, double memV){
  double fact = 1./(exp((memV+154.9)/11.9)-1.);
  double alpham=0.0063*(memV+154.9)*fact;
  double betam=0.00193*exp(memV/33.1);
  double minf=alpham/(alpham+betam);
  double taum=1./(alpham+betam);

  return (minf-variable)/taum;
}


double hchannelHill(double variable, double memV){
  double mh=1./(1.+exp((memV+75.0)/5.5));
  double aux_taum=exp(-14.59-0.086*memV) + exp(-1.87+0.0701*memV);
  double taum=1./aux_taum;

  return (mh-variable)/taum;
}


/**********************************************************/
/*                The Dynamics start here                 */
/**********************************************************/
double syn_sigmoid(double memV){
  return 1./(1.+exp(-0.5*(memV-20.)));
}
// this function is called multiple more times...
void allderivsE(double *y, double *dy, double IsynSoma, double IsynDend, int index, double gKNa_cond, double gKCa_cond,  double factorM, double factorH, int key_h){


	const double phi_eq = Rpump*naEq*naEq*naEq/(naEq*naEq*naEq+kp3);

  // the pow function raises the first argument to the power of the second
  double NaChannel=0.0;
	NaChannel=pow(Na_steadystate_m_E(y[0]),3)*y[1]*gNa_E*0.5*0.3*(y[0]-55.0);
  //you are raising some voltage dependent gating variable m to the third power
  //then multiplying by h_Na_E (not clear what this is)
  //then you are multiplying by the driving voltage (or whatever it is called, the voltage diff)

	double KChannel=0.0;
	KChannel=y[2]*y[2]*y[2]*y[2]*gK_E*0.5*0.3*(y[0]+100.0);

	double LeakChannelS=0.0;
	LeakChannelS=glE[index]*0.3*(y[0]-VlE[index]);

  double LeakChannelD=0.0;
	LeakChannelD=glE[index]*0.7*(y[5]-VlE[index]);

  double AChannel=0.0;
  AChannel=A_steadystate_m_E(y[0])*A_steadystate_m_E(y[0])*A_steadystate_m_E(y[0]);
  AChannel*=gA*0.5*0.3*y[3]*(y[0]+100.0);

  double KSChannel=0.0;
  KSChannel=gKS*0.5*0.3*y[4]*(y[0]+100.0);

  /*Dendrite below - Soma above*/

  double CaChannel=0.0;
  CaChannel=gCa*0.5*0.7*Ca_act_m_E(y[5])*Ca_act_m_E(y[5])*(y[5]-120.0);

  double NaPChannel=0.0;
  NaPChannel=NaP_act_m_E(y[5])*NaP_act_m_E(y[5])*NaP_act_m_E(y[5]);
  NaPChannel*=gNaP*0.5*0.7*(y[5]-55.0);

  double ARChannel=0.0;
  ARChannel=gAR*0.5*0.7*AR_act_h_E(y[5])*(y[5]+100.0);

  double KCaChannel=0.0;
  KCaChannel=gKCa_cond*0.5*0.7*y[6]/(y[6]+30.)*(y[5]+100.0);
  KCaChannel=KCaChannel;

  double KNaChannel=0.0;
  KNaChannel=gKNa_cond*0.5*0.3*KNa_winf_E(y[7])*(y[0]+100.0);
  KNaChannel=KNaChannel;

  // if(keyKCa==1) cout << "WITH " << KCaChannel << endl;
  // if(keyKCa==0) cout << "WITHOUT " << KCaChannel << endl;

  double MChannel=0.0;
  MChannel=gM*0.5*0.3*y[11]*(y[0]+100);

  // we are going to have to track down every mention of gH maybe?
  //depends on what happens with this HChannel
  double Hchannel=0.0;
  Hchannel=key_h*gH*0.5*0.7*y[12]*(y[5]+45.0); //here we are adding our factorH blocking variable

  const double iCsoma = 1./0.5/0.3; //wtf is this
  const double iCdend = 1./0.5/0.7;
  //looks like we are summing all of the various currents together to get the current for the soma and dendrite
  //we are also adding the current that is coming from synaptic inputs (see last entries), I believe these are calculated before this function is called
  double iSoma = -LeakChannelS-NaChannel-KChannel-0.9*AChannel-0.7*KSChannel-0.56*KNaChannel-factorM*0.083*MChannel-gsd[index]*(y[0]-y[5])+IsynSoma;

  double iDend = -LeakChannelD-CaChannel-1.*NaPChannel-1.*ARChannel-1.*KCaChannel-factorH*0.0115*Hchannel-gsd[index]*(y[5]-y[0])+IsynDend;

/*    double iSoma = -LeakChannelS-NaChannel-KChannel-0.9*AChannel-0.7*KSChannel-0.56*KNaChannel-0.083*MChannel-gsd[index]*(y[0]-y[5])+IsynSoma;
  double iDend = -LeakChannelD-CaChannel-1.*NaPChannel-1.*ARChannel-1.*KCaChannel-0.012*Hchannel-gsd[index]*(y[5]-y[0])+IsynDend;
*/

  dy[0] = iSoma*iCsoma;
  dy[1] = Na_inactgate_h_E(y[1],y[0]);
  dy[2] = K_actgate_n_E(y[2],y[0]);
  dy[3] = A_inactgate_h_E(y[3],y[0]);
  dy[4] = KS_actgate_m_E(y[4],y[0]);

  dy[5] = iDend*iCdend;

  dy[6] = -alphaCa*CaChannel-y[6]/tauCa;
  double na3=y[7]*y[7]*y[7];
  dy[7] = -naAlpha*(NaChannel+NaPChannel)-Rpump*na3/(na3+kp3) + phi_eq;

  dy[8] = 3.48*syn_sigmoid(y[0]) - y[8]/2.0;
  dy[9] = 0.5*y[10]*(1.-y[9])-y[9]/100.;
  dy[10] = 3.48*syn_sigmoid(y[0]) - y[10]/2.;

  dy[11] = mcurrentMcCormick(y[11],y[0]);
  dy[12] = hchannelHill(y[12],y[5]);
}

/**/
//idk understand how this method works... I think it is like eulers method but maybe taking multiple smaller steps?
// it is a fourth order Runge-Kutta method, you take the derivitive at the beggining, end, and first and third quarters of the interval, and average them (look it up)
void rk4E(double curr, double cur, int ind, double gKNa_cond, double gKCa_cond,  double keyKNa, double keyKCa, double key_h){
	int i;
	double hh=dt*0.5;
	double h6=dt*0.166666666666666667;
	for (i=0;i<nVarsE;i++)
		yt[i]=varsE[i]+hh*dVarsE[i];
	allderivsE(yt,dyt,curr,cur,ind, gKNa_cond, gKCa_cond, keyKNa, keyKCa,key_h);
	for (i=0;i<nVarsE;i++)
		yt[i]=varsE[i]+hh*dyt[i];
	allderivsE(yt,dyt,curr,cur,ind, gKNa_cond, gKCa_cond, keyKNa, keyKCa,key_h);
	for (i=0;i<nVarsE;i++) {
		yt[i]=varsE[i]+dt*dym[i];
		dym[i] += dyt[i];
	}
	allderivsE(yt,dyt,curr,cur,ind, gKNa_cond, gKCa_cond, keyKNa, keyKCa,key_h);
	for (i=0;i<nVarsE;i++)
		varsE[i] += h6*(dVarsE[i]+dyt[i]+2.0*dym[i]); //this is the actual averaging of the values
}
/**/
// this is the function originally called in main
void rungeKutta4E(double timesim, int index, double currExt1, double currExt2, double gKNa_cond, double gKCa_cond, double keyKNa, double keyKCa, double key_h){
  // this y vector is from the main function, idk how it got here...
	varsE[0]=y[0]; //memV_SomaE
	varsE[1]=y[1]; //h_Na_E, i think these have to be activation gates
	varsE[2]=y[2]; //n_K_E
  varsE[3]=y[3]; //h_A_E
  varsE[4]=y[4]; //m_KS_E
  varsE[5]=y[5]; //memV_DendE
  varsE[6]=y[6]; //conc_Ca_E
  varsE[7]=y[7]; //conc_Na_E
  varsE[8]=y[8]; //sAmpaE
  varsE[9]=y[9]; //sNmdaE
  varsE[10]=y[10]; //xNmdaE
  varsE[11]=y[11]; //mcurrent
  varsE[12]=y[12]; //hcurrent
  // note, I think dVarsE holds the derivitive of the varibles in VarsE, or how much they change each step
	allderivsE(varsE,dVarsE,currExt1,currExt2,index, gKNa_cond, gKCa_cond, keyKNa, keyKCa, key_h); //I don't understand how to feed variables into functions in C, 
  rk4E(currExt1,currExt2,index, gKNa_cond, gKCa_cond, keyKNa, keyKCa,  key_h);
	y[0]=varsE[0];
	y[1]=varsE[1];
	y[2]=varsE[2];
  y[3]=varsE[3];
  y[4]=varsE[4];
  y[5]=varsE[5];
  y[6]=varsE[6];
  y[7]=varsE[7];
  y[8]=varsE[8];
  y[9]=varsE[9];
  y[10]=varsE[10];
  y[11]=varsE[11];
  y[12]=varsE[12];
}
/**/





#endif

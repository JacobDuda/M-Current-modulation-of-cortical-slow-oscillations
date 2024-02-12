#ifndef Variables_H
#define Variables_H


// long seed=124579; /*For Randomnumber.h*/
const int NumOfCells=1280; /*For Network.h*/
const int numcells=1280;
const int numIcells=256;
const int numEcells=1024;
#define vth -20


long seed;


#define dt 0.06
#define timesimulation 80000
//
#include "Randomnumber.h"


// #define gEE_A 0.0054
// #define gEE_N 0.0009
//
// #define gEI_A 0.00225
// #define gEI_N 0.0005
//
// #define gIE 0.00415
// #define gII 0.000165

#define gEE_A 0.0030
#define gEE_N 0.002

#define gEI_A 0.0005625
#define gEI_N 0.000125

#define gIE 0.01366
#define gII 0.00066



/*Below, parameters form Compte code:*/
/* #define EE_CONDUCTANCE_AMPA 0.108 */
/* #define EE_CONDUCTANCE_NMDA 0.018*/
/* #define EI_CONDUCTANCE_AMPA 0.045*/
/* #define EI_CONDUCTANCE_NMDA 0.01*/
/* #define IE_CONDUCTANCE 0.083*/
/* #define II_CONDUCTANCE 0.0033*/


/*Ecells
eCondA =  pn->EEcondAMPA*EESynInputAMPA[i]/pn->EsparseConn;
eCondN = condNMDA(ECellRing[i].Vdendrite,pn->EEcondNMDA)*EESynInputNMDA[i]/pn->EsparseConn; //input to dendrite
iCond = pn->IEcond*IESynInput[i]/pn->IsparseConn*frac;
*/
/*Icells
eCondA =  pn->EIcondAMPA*EISynInputAMPA[i]/pn->EsparseConn/frac;
eCondN = condNMDA(ICellRing[i].potential,pn->EIcondNMDA)*EISynInputNMDA[i]/pn->EsparseConn/frac;
iCond = pn->IIcond*IISynInput[i]/pn->IsparseConn*frac;
*/

#define vSynGABA -70.0
#define vSynGlu 0.0

/*Pointers for the network*/
int (*Wee)[100];
int (*Wei)[100];
int (*Wie)[100];
int (*Wii)[100];

double *EESynInputAMPA;
double *EESynInputNMDA;
double *IESynInput;
double *EISynInputAMPA;
double *EISynInputNMDA;
double *IISynInput;

/*Variables for Pyramidal Neurons*/
#define nVarsE 13
double *y;double *dydt;double *varsE;double *dVarsE;
double *dym;double *dyt;double *yt;
/**/

/*Variables for Interneurons Neurons*/
#define nVarsI 4
double *x;double *dxdt;double *varsI;double *dVarsI;
double *dxm;double *dxt;double *xt;
/**/

/*parameters for concentration of Na*/
#define kp3 3375.0
#define Rpump 0.018
#define naEq 9.5
#define naAlpha 0.01
/**/
/*parameters for concentration of Ca*/
#define alphaCa 0.005
#define tauCa 150.0
/**/


#define gM 1.0
#define gH 1.0

/*Conductances for PY-Cells*/
#define gNa_E 50.0
#define gK_E 10.5
#define gA 1.0
#define gKS 0.576
#define gCa 0.43
#define gNaP 0.0686
#define gAR 0.0257
#define gKCa 0.57
#define gKNa 1.33
/**/

/*Conductances for Interneurons*/
#define gNa_I 35.0
#define gK_I 9.0
/**/



/*Vector for dynamic variables E-Cells*/
double *memV_SomaE; double *h_Na_E;
double *n_K_E; double *h_A_E;
double *m_KS_E; double *memV_DendE;
double *conc_Ca_E; double *conc_Na_E ;
double *mcurrent; double *hcurrent;
/**/
/*Vector for dynamic variables I-Cells*/
double *memV_SomaI;
double *h_Na_I;
double *n_K_I;
/**/

/*Vector for synapses*/
double *sAmpaE; double *xNmdaE;
double *sNmdaE; double *sGabaI;
double *auxGabaI;
/**/

/*Aux for RK4E*/
double **aux_dvarsE;
double **aux_dxm;
double **aux_dxt;
double **aux_xt;
/**/
/*Aux for RK4I*/
double **aux_dvarsI;
double **aux_dym;
double **aux_dyt;
double **aux_yt;
/**/
/*Random Variables*/
/*Random variables*/
double *glI;double *VlI;
double *glE;double *VlE;double *gsd;
/**/

int *spikesindexes;
double *timeindexes;
double **monitoringneurons;
double *LFP;

/*Start dynamic variables*/
void startpointers(){
  /**/
  y= new double [nVarsE];dydt= new double[nVarsE];varsE = new double[nVarsE];
  dVarsE = new double[nVarsE];dym = new double[nVarsE];dyt = new double[nVarsE];
  yt = new double[nVarsE];
  /**/
  /**/
  x= new double [nVarsI];dxdt= new double[nVarsI];varsI = new double[nVarsI];
  dVarsI = new double[nVarsI];dxm = new double[nVarsI];dxt = new double[nVarsI];
  xt = new double[nVarsI];
  /**/

  /**/
  /**/
  memV_SomaE= new double[NumOfCells]; h_Na_E= new double[NumOfCells];
  n_K_E= new double[NumOfCells]; h_A_E= new double[NumOfCells];
  m_KS_E= new double[NumOfCells]; memV_DendE= new double[NumOfCells];
  conc_Ca_E= new double[NumOfCells]; conc_Na_E = new double[NumOfCells];
  mcurrent=new double[NumOfCells]; hcurrent=new double[NumOfCells];
  /**/
  /**/
  memV_SomaI= new double[NumOfCells];
  h_Na_I= new double[NumOfCells];
  n_K_I= new double[NumOfCells];
  /**/
  /**/
  /*SYNAPSES*/
  sAmpaE = new double[NumOfCells]; xNmdaE = new double[NumOfCells];
  sNmdaE = new double[NumOfCells]; sGabaI = new double[NumOfCells];
  auxGabaI = new double[NumOfCells];
  /**/

  // /*Aux for RK4E*/
  aux_dvarsE = new double*[NumOfCells];
  for(int i=0;i<NumOfCells;i++) aux_dvarsE[i]= new double[nVarsE];
  aux_dym = new double*[NumOfCells];
  for(int i=0;i<NumOfCells;i++) aux_dym[i]= new double[nVarsE];
  aux_dyt = new double*[NumOfCells];
  for(int i=0;i<NumOfCells;i++) aux_dyt[i]= new double[nVarsE];
  aux_yt = new double*[NumOfCells];
  for(int i=0;i<NumOfCells;i++) aux_yt[i]= new double[nVarsE];
  /**/
  /*Aux for RK4I*/
  aux_dvarsI = new double*[NumOfCells];
  for(int i=0;i<NumOfCells;i++) aux_dvarsI[i]= new double[nVarsI];
  aux_dxm = new double*[NumOfCells];
  for(int i=0;i<NumOfCells;i++) aux_dxm[i]= new double[nVarsI];
  aux_dxt = new double*[NumOfCells];
  for(int i=0;i<NumOfCells;i++) aux_dxt[i]= new double[nVarsI];
  aux_xt = new double*[NumOfCells];
  for(int i=0;i<NumOfCells;i++) aux_xt[i]= new double[nVarsI];
  /**/
  /**/
  for(int i=0;i<NumOfCells;i++){
    memV_SomaE[i]=-65.0; h_Na_E[i]=0.0;
    n_K_E[i]=0.0; h_A_E[i]=0.0;
    m_KS_E[i]=0.0; memV_DendE[i]=memV_SomaE[i];
    conc_Ca_E[i]=0.0; conc_Na_E[i]=10.0;
    /**/
    memV_SomaI[i]=-62.35;
    h_Na_I[i]=0.0;
    n_K_I[i]=0.0;
    /**/
    sAmpaE[i]=0.0; xNmdaE[i]=0.0;
    sNmdaE[i]=0.0; sGabaI[i]=0.0;
    auxGabaI[i]=0.0;
  }
  /**/
	y[0]=-65.0; /*V Soma*/ y[1]=0.0; /*h for Na*/
	y[2]=0.0; /*n for K*/ y[3]=0.0; /*h for A*/
  y[4]=0.;  /*m for KS*/ y[5]=-65.0; /*V dendrite*/
  y[6]=0.0; /*Ca conc*/  y[7]=0.0; /*Na Conc*/
  y[8]=0.0; /*Ca conc*/  y[9]=0.0; /*Na Conc*/
  y[10]=0.0; /*Ca conc*/
  /**/
  x[0]=-62.35; /*V Soma Interneurons*/
  x[1]=0.0; /*h for Na*/
  x[2]=0.0; /*n for K*/
  x[3]=0.0;
  /**/
  glI=new double[numIcells];
  VlI=new double[numIcells];
  glE=new double[numEcells];
  VlE=new double[numEcells];
  gsd=new double[numEcells];

  EESynInputAMPA=new double[numEcells]; /*E receives from E*/
  EESynInputNMDA=new double[numEcells]; /*E receives from E*/
  IISynInput=new double[numIcells]; /*I receives from I*/
  IESynInput=new double[numEcells]; /*E receives from I*/
  EISynInputAMPA=new double[numIcells]; /*I receives from E*/
  EISynInputNMDA=new double[numIcells]; /*I receives from E*/


  spikesindexes=new int[timesimulation];
  timeindexes=new double[timesimulation];
  for(int i=0;i<timesimulation;i++){
    spikesindexes[i]=-1;
    timeindexes[i]=-1.0;
  }

  monitoringneurons=new double*[numbermonitoringneurons*2];
  for(int i=0;i<numbermonitoringneurons*2;i++) monitoringneurons[i]= new double[timesimulation];

  LFP=new double[timesimulation];


  /**/
  for(int i=0;i<numEcells;i++){
    glE[i]=0.01+0.001*gaussian_generator(0,1,seed);
    VlE[i]=-60.95+0.3*gaussian_generator(0,1,seed);
    gsd[i]=1.75+0.1*gaussian_generator(0,1,seed);
  }
  for(int i=0;i<numIcells;i++){
    glI[i]=0.0205+0.0005*gaussian_generator(0,1,seed);
    VlI[i]=-63.8+0.15*gaussian_generator(0,1,seed);
  }
}

void ComputeSynpases(){//i think we are just resetting any old values
  for(int i=0;i<numEcells;i++){EESynInputAMPA[i]=0.0;
    EESynInputNMDA[i]=0.0;IESynInput[i]=0.0;}
  for(int i=0;i<numIcells;i++){EISynInputAMPA[i]=0.0;
    EISynInputNMDA[i]=0.0;IISynInput[i]=0.0;}
  int jj=0;
  /*excitatory*/
  for(int i=0;i<numEcells;i++){ //for all excitatory cells
    jj=0;
    while(Wee[i][jj]<numEcells){ //we are adding all the activations from all jj excitatory cells onto cell i?
      EESynInputAMPA[Wee[i][jj]]+=sAmpaE[i];//sEAMPA[i];
      EESynInputNMDA[Wee[i][jj]]+=sNmdaE[i];//sENMDA[i];
      jj++;
    }
    jj=0;
    while(Wei[i][jj]<numIcells){ //same thing adding all excitatory input onto inhibitory cells
      EISynInputAMPA[Wei[i][jj]]+=sAmpaE[i];//sEAMPA[i];
      EISynInputNMDA[Wei[i][jj]]+=sNmdaE[i];//sENMDA[i];
      jj++;
    }
  }/*end lop exc*/
  /*Inhibitory*/
  for(int i=0;i<numIcells;i++){
    jj=0;
    while(Wie[i][jj]<numEcells){
      IESynInput[Wie[i][jj]]+=sGabaI[i];//sI[i];
      jj++;
    }
    jj=0;
    while(Wii[i][jj]<numIcells){
      IISynInput[Wii[i][jj]]+=sGabaI[i];//sI[i];
      jj++;
    }
  }/*end loop inh*/
}

void deletepointers(){
  delete [] y;delete [] dydt;delete [] varsE;
  delete [] dVarsE;delete [] dym;delete [] dyt;
  delete [] yt;
  /**/
  delete [] x;delete [] dxdt;delete [] varsI;
  delete [] dVarsI;delete [] dxm;delete [] dxt;
  delete [] xt;
  /**/
  for(int i=0;i<NumOfCells;i++) aux_dvarsI[i];delete [] aux_dvarsI;
  for(int i=0;i<NumOfCells;i++) aux_dym[i];delete [] aux_dym;
  for(int i=0;i<NumOfCells;i++) aux_dyt[i];delete [] aux_dyt;
  for(int i=0;i<NumOfCells;i++) aux_yt[i];delete [] aux_yt;
  for(int i=0;i<NumOfCells;i++) delete aux_dvarsE[i];delete [] aux_dvarsE;
  for(int i=0;i<NumOfCells;i++) delete aux_dxm[i];delete [] aux_dxm;
  for(int i=0;i<NumOfCells;i++) delete aux_dxt[i];delete [] aux_dxt;
  for(int i=0;i<NumOfCells;i++) delete aux_xt[i];delete [] aux_xt;

  for(int i=0;i<numbermonitoringneurons*2;i++) delete monitoringneurons[i]; delete  [] monitoringneurons;
  /**/
  /**/
  delete [] memV_SomaE; delete [] h_Na_E;
  delete [] n_K_E; delete [] h_A_E;
  delete [] m_KS_E; delete [] memV_DendE;
  delete [] conc_Ca_E; delete [] conc_Na_E;
  delete [] mcurrent;
  delete [] hcurrent;
  /**/
  /**/
  delete [] memV_SomaI;
  delete [] h_Na_I;
  delete [] n_K_I;
  /**/
  delete [] sAmpaE; delete [] xNmdaE;
  delete [] sNmdaE; delete [] sGabaI;
  delete [] auxGabaI;

  delete [] Wee;
  delete [] Wei;
  delete [] Wie;
  delete [] Wii;

  delete [] EESynInputAMPA;
  delete [] EESynInputNMDA;
  delete [] IESynInput;
  delete [] EISynInputAMPA;
  delete [] EISynInputNMDA;
  delete [] IISynInput;

  delete [] spikesindexes;
  delete [] timeindexes;
  delete [] LFP;
}
/**/














#endif

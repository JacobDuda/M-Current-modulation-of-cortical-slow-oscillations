#include<iostream>
#include "Randomnumber.h"
#include "Channels_ECells.h"
#include "Channels_ICells.h"
#include "Variables.h"
#include "Network.h"
#include "NeuronsStructure.h"

using namespace std;




int main (int argc, char *argv[]){



  double transient_samples=10000.0; 
  double factorM=0.;

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /*&&*/if(argc != 3){
	/*&&*/	cout << "#########-ERROR-########" << endl;
	/*&&*/	cout << "# Give 2 paramenters as input:" << endl;
	/*&&*/	cout << "# 1) Percentage factor for M current (betw 0-1) " << endl;
  /*&&*/	cout << "# 2) seed " << endl;
	/*&&*/	cout << "# Example: ./exct 0.1 1221314" << endl;
	/*&&*/	cout << "########################" << endl;
	/*&&*/	return 1;
	/*&&*/}
  factorM=(double) atof(argv[1]);
  seed=(long) atoi(argv[2]);

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  const double KNa_percentage= 1.0; /*1st argument - Conductance for IKCa */
  const double gKNa_conductance= 1.37;/*2nd argument - Conductance for IKNa */
  const double KCa_percentage= 1.0; /*3rd argument - Conductance for IKCa */
  const double gKCa_conductance= 0.57;/*4th argument - Conductance for IKNa */
  const int keychannel=1;
  long seed1=seed;
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/


  startpointers();


  /*Only 30% of py neurons have the H channel*/  // maybe this is going to be important later, looks like we are just adding an H current to 30% of py neurons
  int *neuronsH=new int[300];
  for(int i=0;i<300;i++){ // this is just a for loop, for i in range(0,300)
    int auxkey=0;
  	if(i==0) neuronsH[i]=int(ran(&seed)*1024);
    if(i!=0){
      do{
        neuronsH[i]=int(ran(&seed)*1024);
        for(int k=0;k<i;k++){
          if(neuronsH[i]==neuronsH[k]){
            auxkey=1;
            break;
          }
          else auxkey=0;
        }
      }
      while(auxkey==1);
    }
  }
  /**/


	construct_network();

  int *electrode=new int [2];

  electrode[0]=400+int(ran(&seed)*400);
  electrode[1]=600+int(ran(&seed)*400);

  int countLPF=0;
  int countSingleNeuron=0;
  int auxindex=0;


  FILE *outf, *outf2, *outf3, *outf4,*outf5;
	char *Outfile;
  char *Outfile2;
  char *Outfile3;
  char *Outfile4;
  char *Outfile5;

	Outfile=new char[500];
  Outfile2=new char[500];
  Outfile3=new char[500];
  Outfile4=new char[500];
  Outfile5=new char[500];

  sprintf(Outfile,"Raster_factorM_%.3lf_seed%ld.txt",factorM,seed1);
  sprintf(Outfile2,"LFP_factorM_%.3lf_seed%ld.txt",factorM,seed1);
  sprintf(Outfile3,"NEURONS_factorM_%.3lf_seed%ld.txt",factorM,seed1);
  sprintf(Outfile4,"FR_GLOBAL_factorM_%.3lf_seed%ld.txt",factorM,seed1);
  sprintf(Outfile5,"Synchrony__factorM_%.3lf_seed%ld.txt",factorM,seed1);


  outf = fopen(Outfile,"w");
  outf2 = fopen(Outfile2,"w");
  outf3 = fopen(Outfile3,"w");



  const int number_sync_cells=512;

  double **monitoring_cell;
  double meanV=0.0;
  double meanV_sqr=0.0;
  double aux_meanV=0.0;
  monitoring_cell = new double*[number_sync_cells];
  for(int i=0;i<number_sync_cells;i++) monitoring_cell[i]= new double[2];


  for(int t=1;t<timesimulation;t++){ //for t in range(0,timesteps)
    double auxsumLFP=0.0;
    int aux_FR=0;
    double currtimesimu=t*dt*0.001; /*current time simulation*/

    //this is another function stored in the variables.h file, i think it just adds the activations from all cells onto eachother... but check this
    ComputeSynpases(); // i think we are just computing what is happening at all the synapses when neurotransmitters land there

    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
    /*Excitatory looping over time - EXC DYNAMICS*/
    /**/for(int j=0;j<numEcells;j++){ //for i in number of cells 
    /**/  double eCondAmpa=gEE_A*EESynInputAMPA[j]; // AMBA and GABA conductances are not voltage dependent, NMDA conductance is
    /**/  double eCondNMDA=gEE_N*EESynInputNMDA[j]/(1.0 + 0.2801*exp(-0.062*memV_DendE[j]));
    /**/  double iCondGABA=gIE*IESynInput[j];
    /**/  double currtodend=(eCondAmpa+eCondNMDA)*(vSynGlu-memV_DendE[j]); //the current to every dendrite is the conductances of the AMPA and NMDA receptors multiplied by their reversal potential compared to current Vmem
    /**/  double currtosoma=iCondGABA*(vSynGABA-memV_SomaE[j]); //same thing as above, but the GABA cells are landing on the soma
    /**/
    /**/  y[0]=memV_SomaE[j];y[1]=h_Na_E[j]; // all of these weird entries are related to "dynamic variables" for e-cells, whatever that means
    /**/  y[2]=n_K_E[j];y[3]=h_A_E[j];
    /**/  y[4]=m_KS_E[j];y[5]=memV_DendE[j];
    /**/  y[6]=conc_Ca_E[j];y[7]=conc_Na_E[j];
    /**/  y[8]= sAmpaE[j];y[9]= sNmdaE[j];
    /**/  y[10]= xNmdaE[j]; y[11]=mcurrent[j]; y[12]=hcurrent[j];
    /**/  for(int z=0;z<nVarsE;z++){ // nVarsE is how many variables are associated with an excitatory cell
    /**/    dVarsE[z]=aux_dvarsE[j][z]; //idk what these are, storing variables ig 
    /**/    dym[z]=aux_dym[j][z];
    /**/    dyt[z]=aux_dyt[j][z];
    /**/    yt[z]=aux_yt[j][z];
    /**/  }
    /**/
          int key_h=0;
          for(int ii=0;ii<300;ii++){
            if(j==neuronsH[ii]) key_h=1; // this is checking if neuron[j] was randomly given H_channels earlier I think
          }
    /**/  rungeKutta4E(currtimesimu,j,currtosoma,currtodend, gKNa_conductance, gKCa_conductance, factorM, 1.0, key_h);

          /*Implementing the local synchrony*/
          // wth is this doing
          if(t>transient_samples){ //transient samples are the number of steps the systems simulates before data collection (warm up)
            if(j>=256 and j<768){
              for(int ii=0;ii<number_sync_cells;ii++){
                if((ii+256)==j){
                  monitoring_cell[j-256][0]+=y[0]*y[0]; //monitoring cell is only called in writing to output, so maybe this is just for some stat... idk
                  monitoring_cell[j-256][1]+=y[0]/(timesimulation-transient_samples-1.0);
                  aux_meanV+=y[0]*(1./number_sync_cells);
                  break;
                }
              }
            }
          }

          if(t>transient_samples){
            if(j<=electrode[0]+25 and j>=electrode[0]-25) auxsumLFP+=abs(currtosoma)+abs(currtodend);
          }

          if(t>transient_samples){

            if(((memV_SomaE[j]-vth)*(y[0]-vth)<0.0) && (y[0]>vth)){
              fprintf(outf,"%lf\t%d\t1\n", (t-transient_samples)*dt*0.001, j);
            }
            
            if(j==500)   fprintf(outf3,"%lf\t%lf\t", (t-transient_samples)*dt*0.001, y[0]);
            if(j==501)   fprintf(outf3,"%lf\t", y[0]);
            if(j==502)   fprintf(outf3,"%lf\t", y[0]);
            if(j==503)   fprintf(outf3,"%lf\t", y[0]);
            if(j==504)   fprintf(outf3,"%lf\t", y[0]);
            if(j==505)   fprintf(outf3,"%lf\t", y[0]);
            if(j==506)   fprintf(outf3,"%lf\t", y[0]);
            if(j==507)   fprintf(outf3,"%lf\t", y[0]);
            if(j==508)   fprintf(outf3,"%lf\t", y[0]);
            if(j==509)   fprintf(outf3,"%lf\n", y[0]);
            
          }

    /**/ // ok it looks like we are updating these dynamic values
    /**/  memV_SomaE[j]=y[0];h_Na_E[j]=y[1];
    /**/  n_K_E[j]=y[2];h_A_E[j]=y[3];
    /**/  m_KS_E[j]=y[4];memV_DendE[j]=y[5];
    /**/  conc_Ca_E[j]=y[6];conc_Na_E[j]=y[7];
    /**/  sAmpaE[j]=y[8];sNmdaE[j]=y[9];
    /**/  xNmdaE[j]=y[10]; mcurrent[j]=y[11]; hcurrent[j]=y[12];
    /**/  for(int z=0;z<nVarsE;z++){
    /**/    aux_dvarsE[j][z]=dVarsE[z];
    /**/    aux_dym[j][z]=dym[z];
    /**/    aux_dyt[j][z]=dyt[z];
    /**/    aux_yt[j][z]=yt[z];
    /**/  }
    /**/}/*End looping EXC*/
    if(t>transient_samples){
      meanV+=(aux_meanV)/(timesimulation-transient_samples-1.0);
      meanV_sqr+=aux_meanV*aux_meanV;
      aux_meanV=0.0;
    }

    /**//*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
    if(t>transient_samples){
      fprintf(outf2,"%lf\t%lf\n", (t-transient_samples)*dt*0.001, auxsumLFP);
    }
    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
    /*Inhibitory looping over time - EXC DYNAMICS*/
    /**/for(int j=0;j<numIcells;j++){
    /**/  double eCondAmpa=gEI_A*EISynInputAMPA[j];
    /**/  double eCondNMDA=gEI_N*EISynInputNMDA[j]/(1.0 + 0.2801*exp(-0.062*memV_SomaI[j]));
    /**/  double iCondGABA=gII*IISynInput[j];
    /**/  double currtodend=(eCondAmpa+eCondNMDA)*(vSynGlu-memV_SomaI[j]);
    /**/  double currtosoma=iCondGABA*(vSynGABA-memV_SomaI[j]);
    /**/
    /**/  x[0]=memV_SomaI[j];
    /**/  x[1]=h_Na_I[j];
    /**/  x[2]=n_K_I[j];
    /**/  x[3]=sGabaI[j];
    /**/  for(int z=0;z<nVarsI;z++){
    /**/    aux_dvarsI[j][z]=dVarsI[z];
    /**/    aux_dxm[j][z]=dxm[z];
    /**/    aux_dxt[j][z]=dxt[z];
    /**/    aux_xt[j][z]=xt[z];
    /**/  }
    /**/
    /**/  rungeKutta4I(currtimesimu,j, currtodend+currtosoma);
          if(t>transient_samples){
            if(((memV_SomaI[j]-vth)*(x[0]-vth)<0.0) && (x[0]>vth))   fprintf(outf,"%lf\t%d\t0\n", (t-transient_samples)*dt*0.001, j);
          }
    /**/  //if(((memV_SomaI[j]-vth)*(x[0]-vth)<0.0))  cout << t*dt*0.001 << "\t" << j << "\t" << 0 << endl;
    /**/
    /**/  memV_SomaI[j]=x[0];h_Na_I[j]=x[1];
    /**/  n_K_I[j]=x[2];sGabaI[j]=x[3];
    /**/  for(int z=0;z<nVarsI;z++){
    /**/    aux_dvarsI[j][z]=dVarsI[z];
    /**/    aux_dxm[j][z]=dxm[z];
    /**/    aux_dxt[j][z]=dxt[z];
    /**/    aux_xt[j][z]=xt[z];
    /**/  }
    /**/}/*End Looping inh*/
    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

	}/*END LOOP TIME*/

  fclose(outf);
  fclose(outf2);
  fclose(outf3);

  // fclose(outf4);



  /*Synchorny Parameter*/
  double variance_glogal = (meanV_sqr/(timesimulation-transient_samples-1.0))-(meanV*meanV);

  double variance_local=0.0;
  for(int i=0;i<number_sync_cells;i++){
    variance_local+=((monitoring_cell[i][0]/(timesimulation-transient_samples-1.0))-(monitoring_cell[i][1]*monitoring_cell[i][1]))/number_sync_cells;
  }

  double synchrony=sqrt(variance_glogal/(variance_local));
  outf5 = fopen(Outfile5,"w");
  fprintf(outf5,"#Global_variance Local_variance Synchronization Parameter\n");
  fprintf(outf5,"%lf\t%lf\t%lf\n",variance_glogal,variance_local,synchrony);
  fclose(outf5);


  for(int i=0;i<number_sync_cells;i++) monitoring_cell[i];
  delete [] monitoring_cell;



  deletepointersMonitoringChannels();
  deletepointers();
  delete [] electrode;

	return 0;
}

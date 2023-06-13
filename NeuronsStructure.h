#ifndef NeuronsStructure_H
#define NeuronsStructure_H


#include "Randomnumber.h"



int* neurons_WITH_KNaChannel;
int* neurons_WITH_KCaChannel;
int* fewNeuronsWITH;
int* fewNeuronsWITHOUT;
int* auxiliar;
int numb_neurons_withKNa;
int numb_neurons_withKCa;


void neurons_with_KNA(double percentageKNa, double percentageKCa){

  /*WarmingUp*/
  for(int i=0;i<1000;i++) double lixo=ran(&seed);
  /**/

  fewNeuronsWITH=new int[numbermonitoringneurons];
  fewNeuronsWITHOUT=new int[numbermonitoringneurons];

  numb_neurons_withKNa= int(percentageKNa*1024.0);
  numb_neurons_withKCa= int(percentageKCa*1024.0);

  neurons_WITH_KNaChannel=new int[numb_neurons_withKNa]; /*Neurons that will have KNA*/
  neurons_WITH_KCaChannel=new int[numb_neurons_withKCa]; /*Neurons that will have KCA*/

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /*Neurons with and without the KNa Channels*/
  for(int u=0;u<numb_neurons_withKNa;u++){
    int auxkey=0;
  	if(u==0) neurons_WITH_KNaChannel[u]=int(ran(&seed)*1024);
  	if(u!=0){
  		do{
  			neurons_WITH_KNaChannel[u]=int(ran(&seed)*1024);
  			for(int k=0;k<u;k++){
  				if(neurons_WITH_KNaChannel[u]==neurons_WITH_KNaChannel[k]){
  					auxkey=1;
  					break;
  				}
  				else auxkey=0;
  			}
  		}
  		while(auxkey==1);
  	}
  }
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /*Neurons with and without the KCa Channels*/
  for(int u=0;u<numb_neurons_withKCa;u++){
    int auxkey=0;
  	if(u==0) neurons_WITH_KCaChannel[u]=int(ran(&seed)*1024);
  	if(u!=0){
  		do{
  			neurons_WITH_KCaChannel[u]=int(ran(&seed)*1024);
  			for(int k=0;k<u;k++){
  				if(neurons_WITH_KCaChannel[u]==neurons_WITH_KCaChannel[k]){
  					auxkey=1;
  					break;
  				}
  				else auxkey=0;
  			}
  		}
  		while(auxkey==1);
  	}
  }
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /*5 neurons WITH channel to be monitored*/
  for(int k=0;k<numbermonitoringneurons;k++){
    int keyaux=0;
    int aux_index=0;
    if(k==0){
      aux_index=int(ran(&seed)*1024);
      fewNeuronsWITH[k]=aux_index;
    }
    else{
      do{
        aux_index=int(ran(&seed)*1024);
        for(int j=0;j<k;j++){
          if(aux_index==fewNeuronsWITH[j]){
            keyaux=1;
            break;
          }
          else{
            keyaux=0;
            fewNeuronsWITH[k]=aux_index;
          }
        }
      }
      while(keyaux==1);
    }
  }
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/


  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /*5 neurons WITHOUT channel to be monitored*/
  auxiliar=new int[numbermonitoringneurons*2];
  for(int k=0;k<(numbermonitoringneurons*2);k++){
    int keyauxx=0;
    if(k==0) auxiliar[k]=int(ran(&seed)*1024);
    else{
      do{
        auxiliar[k]=int(ran(&seed)*1024);
        for(int u=0;u<k;u++){
          if(auxiliar[u]==auxiliar[k]){
            keyauxx=1;
            break;
          }
          else keyauxx=0;
        }
      }
      while(keyauxx==1);
    }
  }

  for(int k=0;k<(numbermonitoringneurons*2);k++){
    if(k<numbermonitoringneurons) fewNeuronsWITH[k]=auxiliar[k];
    else{
      fewNeuronsWITHOUT[k-numbermonitoringneurons]=auxiliar[k];
    }
  }
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/


  // for(int i=0;i<numb_neurons_withKNa;i++){
  //   cout << i << "\t" << neurons_WITH_KNaChannel[i] << endl;
  // }
  // cout << "####################" << endl;
  // for(int i=0;i<numb_neurons_withKCa;i++){
  //   cout << i << "\t" << neurons_WITH_KCaChannel[i] << endl;
  // }
  // cout << "####################" << endl;
  // for(int i=0;i<numb_neurons_withKCa;i++){
  //   for(int j=0;j<numb_neurons_withKCa;j++){
  //     if(neurons_WITH_KCaChannel[i]==neurons_WITH_KCaChannel[j] and i!=j){
  //       cout << "ERROR! Neurons for monitoring with and without are the same!" << endl;
  //       cout << "With: " << i << " " << neurons_WITH_KCaChannel[i] << "\tWithout: " << j << " " << neurons_WITH_KCaChannel[j] << endl;
  //     }
  //   }
  // }


  // /*CONSISTENCY*/
  // for(int i=0;i<numbermonitoringneurons;i++){
  //   for(int j=0;j<numbermonitoringneurons;j++){
  //     if(fewNeuronsWITH[i]==fewNeuronsWITHOUT[j]){
  //       cout << "ERROR! Neurons for monitoring with and without are the same!" << endl;
  //       cout << "With: " << i << " " << fewNeuronsWITH[i] << "\tWithout: " << j << " " << fewNeuronsWITHOUT[j] << endl;
  //     }
  //   }
  // }
  // /**/
  // cout << "####################" << endl;
  // for(int i=0;i<numbermonitoringneurons;i++){
  //   cout << i << "\t WITH: " << fewNeuronsWITH[i] << "\t WITHOUT: " << fewNeuronsWITHOUT[i] << endl;
  // }

}/*END VOID neurons_with_KNA*/


void deletepointersMonitoringChannels(){
  delete [] neurons_WITH_KNaChannel;
  delete [] neurons_WITH_KCaChannel;
  delete [] fewNeuronsWITH;
  delete [] auxiliar;
  delete [] fewNeuronsWITHOUT;
}









#endif

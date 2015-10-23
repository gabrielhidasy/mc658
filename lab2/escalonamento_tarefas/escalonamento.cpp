/* Copyright 2014, Gurobi Optimization, Inc. */

#include <math.h>
#include <iostream> 
#include <fstream> 
#include <cstdlib>
#include <cstring>
#include<lemon/math.h>
#include "gurobi_c++.h"
using namespace std;

int main(int argc,char *argv[]) {
  srand48(1);
  if (argc!=2) {cout<<endl<<"Usage: "<< argv[0]<<"  <filename>"<< endl << endl;
    cout << "Example:      " << argv[0] << " arq1.in" << endl << endl; exit(0);}
	
  ifstream ifile;
  ifile.open(argv[1]);  if (!ifile) return(false);  

  int n, m;
  ifile >> n; ifile >> m;
  vector<double> jobsize(n);
  vector<double> machinespeed(m);
  for(int i=0; i<n; i++) ifile >> jobsize[i];
  for(int i=0; i<m; i++) ifile >> machinespeed[i];
  ifile.close();

  cout << endl;
  cout << "Numero de tarefas: " << n << endl;
  cout << "Tamanho das tarefas" << endl;
  for(int i=0; i<n; i++)
    cout << "T_" << i << ": " << jobsize[i] << endl;
  cout << "Velocidade das maquinas" << endl;
  for(int i=0; i<m; i++)
    cout << "M_" << i << ": " << machinespeed[i] << endl;
  cout << endl << endl;
  
  try {
/*--------------------------------------------------------------------------------- */	
/*--------------------------- ALTERE DAQUI PARA ABAIXO ---------------------------- */
/*--------------------------------------------------------------------------------- */	

    // Voce deve atualizar esta variavel com valores 0 e 1 de tal forma que
    // tarefa_maquina[i][j] = 1 SE_E_SOMENTE_SE a tarefa i foi escalonada na maq. j
    vector<vector<int> > tarefa_maquina(n, vector<int>(m));

    // mais abaixo eh feito um escalonamento aleatorio.
    // Nao se esqueca de alterar a geracao do escalonamento aleatorio para o
    // escalonamento obtido pelo seu programa linear inteiro


    // cabecalho comum para os programas lineares (descomente)
     int seed=0;
     GRBEnv env = GRBEnv();
     GRBModel model = GRBModel(env);
     model.getEnv().set(GRB_IntParam_Seed, seed);
     model.set(GRB_StringAttr_ModelName, "Escalonamento de Tarefas"); // prob. name
     model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // is a minimization problem  

	GRBVar min;
	
	GRBVar x[m][n];

	min = model.addVar(0, 30000, 1 , GRB_CONTINUOUS, "min");

	
	int i,j;
	for(i=0;i<m;i++){
		for(j=0;j<n;j++) {
			char name[100];
			sprintf(name,"I%dM%d",i,j);
			x[i][j] = model.addVar(0, 1, 0 , GRB_BINARY, name);	
		}
	}
		
	model.update();

	//cada tarefa em 1 máquina
	
	for (int i=0; i<n; i++){
		GRBLinExpr sum2;
      		for (int j=0; j<m; j++){
			sum2= sum2 + x[j][i];
		}
		model.addConstr(sum2 == 1);
	}	

	//min eh pior caso
	for(int j=0; j<m; j++){
		GRBLinExpr sum3;
		for (i=0; i<n ; i++) {
			sum3 = sum3 + jobsize[i]*x[j][i]/ machinespeed[j];
		}
		model.addConstr(min >= sum3);
	}

	model.update();
	model.write("model.lp");
	system("cat model.lp");

	model.optimize();
	// Verifica se obteve solucao otima
	if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
       cout << "Erro, sistema impossivel" << endl;
       exit(1);
	}
    for (int i=0; i<n; i++){
      for (int j=0; j<m; j++){
	if(x[j][i].get(GRB_DoubleAttr_X)>0.999)
		tarefa_maquina[i][j] = 1; 
	else
		tarefa_maquina[i][j] = 0; 
	}
	}
 

	  
/*--------------------------------------------------------------------------------- */	
/*--------------------------- ALTERE DAQUI PARA CIMA ------------------------------ */
/*--------------------------------------------------------------------------------- */	

    cout << "\n\n";
      double makespan=0;
      cout << "Escalonamento obtido das tarefas nas maquinas\n";
      cout << "Notacao para tarefas : T_id(tamanho_original)\n";
      cout << "Notacao para maquinas: M_id(velocidade)\n\n";
      
      for (int j=0; j<m; j++) {
	double tmaq=0.0;
	cout << "M_"<<j<<"("<< machinespeed[j] << ") : ";
	for(int i=0; i<n; i++) {
	  if (tarefa_maquina[i][j] == 1) {
	    cout << "T_" << i << "(" << jobsize[i] << ") ";
	    tmaq += jobsize[i] / machinespeed[j];
	  }
	}
	cout << endl << "Tempo gasto pela maquina M_" <<j << ": " << tmaq << endl<<endl;
	if (tmaq>makespan) makespan = tmaq;
      }
      cout << "\nTempo para completar todas as tarefas: " << makespan << "s\n\n";
      //  } else  cout << "No solution" << "\n";
  }
  // catch (GRBException e) {
  //   cout << "Error code = " << e.getErrorCode() << endl;
  //   cout << e.getMessage() << endl;
  // }
  catch (...) {
    cout << "Exception during optimization" << endl;
  }
  return 0;
}


// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
// Problems with connectivity: Minimum Cost k-paths (edge disjoint)
#include <gurobi_c++.h>
#include <iostream>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <gurobi_c++.h>
#include <iostream>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <queue>
#include <lemon/list_graph.h>
#define BC_EPS 0.001
#define BC_INF 1000000000000.0
#include "mygraphlib.h"
#include <string>
#include "myutils.h"
#include <lemon/concepts/digraph.h>
#include <lemon/preflow.h>
using namespace lemon;
using namespace std;
#define EPS 0.00001


int cutcount = 0;

class Problem_Data {
public:
  Problem_Data(ListGraph &g,
	       NodeName &nodename,
	       NodePos &posx,
	       NodePos &posy,
	       vector<Node> &s,
	       vector<Node> &t,
	       vector<double> &tmax,
	       vector<double> &q,
	       EdgeWeight &custo,
	       EdgeWeight &capacidade,
	       EdgeWeight &latencia);
  ListGraph &g;
  int NNodes,NEdges,NPairs;
  NodeName &nodename;
  NodePos &posx;
  NodePos &posy;
  vector<Node> &s;
  vector<Node> &t;
  vector<double> &Tmax;
  vector<double> &q;
  EdgeWeight &custo; 
  EdgeWeight &capacidade; 
  EdgeWeight &latencia;
  vector<vector<Node> > BestSol;
	
	
  double BestVal;
  double BestLB;
};

Problem_Data::Problem_Data(ListGraph &graph,
			   NodeName &vname,
			   NodePos &posicaox,
			   NodePos &posicaoy,
			   vector<Node> &source,
			   vector<Node> &tail,
			   vector<double> &tempomax,
			   vector<double> &qtd,
			   EdgeWeight &cust,
			   EdgeWeight &capac,
			   EdgeWeight &lat):
  g(graph), 
  nodename(vname), 
  posx(posicaox), 
  posy(posicaoy), 
  s(source), t(tail), 
  Tmax(tempomax), 
  q(qtd),
  custo(cust),
  capacidade(capac),
  latencia(lat)
{
  BestVal = 0;
  BestLB = 0;
}

void extract_solution(Problem_Data &G,
		      ListGraph::EdgeMap< vector <GRBVar> > &g,
		      ListGraph::EdgeMap< vector <GRBVar> > &b) {
  // for (int z = 0; z < G.NPairs; z++) {
  //   for (ListGraph::EdgeIt e(G.g); e != INVALID; ++e) {
  //     //Para cada aresta desse no
  //     if (g[e][z].get(GRB_DoubleAttr_X) > 0.9 ||
  // 	  b[e][z].get(GRB_DoubleAttr_X) > 0.9) {
  // 	G.BestVal += G.custo[e]*G.q[z];
  //     } 
  //   }
  // } //Obtem a best bound manualmente
  G.BestSol.resize(G.NPairs);
  for (int z = 0; z < G.NPairs; z++) {
    vector <Node> npath[2];
    for (ListGraph::EdgeIt e(G.g); e != INVALID; ++e) {
      //Para cada aresta desse no
      if (g[e][z].get(GRB_DoubleAttr_X) > 0.9) {
	npath[0].push_back(G.g.u(e));
	npath[1].push_back(G.g.v(e));
      }
      if (b[e][z].get(GRB_DoubleAttr_X) > 0.9) {
	npath[0].push_back(G.g.v(e));
	npath[1].push_back(G.g.u(e));
      }
    }
    Node ns = G.s[z];
    while (G.BestSol[z].size() != npath[0].size()) {
      G.BestSol[z].push_back(ns);
      int c = 0;
      while (npath[0][c] != ns) c++;
      ns = npath[1][c];
    }
    G.BestSol[z].push_back(G.t[z]);
  }
}

bool transmissoes(Problem_Data &G, long maxtime)
{
  int seed = 0;
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.getEnv().set(GRB_IntParam_Seed, seed);
  model.set(GRB_StringAttr_ModelName, "Transmissão");
  model.getEnv().set(GRB_DoubleParam_TimeLimit, maxtime);
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
  // is a minimization problem
  G.NPairs = G.s.size(); //Correcting strange problem in input

  //Two variables for each edge for each pair
  ListGraph::EdgeMap< vector <GRBVar> > g(G.g);
  ListGraph::EdgeMap< vector <GRBVar> > b(G.g);

  char nammeG[100];
  char nammeB[100];
  double tmp = 0;
  //Funcao de minimizacao
  for (ListGraph::EdgeIt e(G.g); e != INVALID; ++e) {
    tmp = G.custo[e];
    g[e].resize(G.NPairs);
    b[e].resize(G.NPairs);
    for (int z = 0; z < G.NPairs; z++) {
      sprintf(nammeG,"EG_%d_%d_%d",G.g.id(G.g.u(e))+1,G.g.id(G.g.v(e))+1,z);
      sprintf(nammeB,"EB_%d_%d_%d",G.g.id(G.g.v(e))+1,G.g.id(G.g.u(e))+1,z);
      g[e][z] = model.addVar(0,1,tmp*G.q[z],GRB_BINARY,nammeG);
      b[e][z] = model.addVar(0,1,tmp*G.q[z],GRB_BINARY,nammeB);
    }
  }
  model.update();
  // Restricao de latencia
  for (int z = 0; z < G.NPairs; z++) {
    GRBLinExpr expr;
    for (ListGraph::EdgeIt e(G.g); e != INVALID; ++e) {
      expr += g[e][z]*G.latencia[e];
      expr += b[e][z]*G.latencia[e];
    }
    model.addConstr(expr <= G.Tmax[z]);
  }
  model.update();
  // Restricao de capacidade
  for (ListGraph::EdgeIt e(G.g); e != INVALID; ++e) {
    GRBLinExpr expr;
    for (int z = 0; z < G.NPairs; z++) {
      expr += g[e][z]*G.q[z];
      expr += b[e][z]*G.q[z];
    }
    model.addConstr(expr <= G.capacidade[e]);
  }
  model.update();
  // Restricao de conservacao de fluxo
  for (int z = 0; z < G.NPairs; z++) { //Para cada par
    for (ListGraph::NodeIt v(G.g); v != INVALID; ++v) {
      //Para cada no
      GRBLinExpr expr;
      for (ListGraph::IncEdgeIt e(G.g, v); e != INVALID; ++e) {
	//Para cada aresta desse no
	if (G.g.u(e) == v) {
	  expr += g[e][z];
	  expr -= b[e][z];
	} else {
	  expr -= g[e][z];
	  expr += b[e][z];
	}
      }
      if (G.g.id(v) == G.g.id(G.s[z])) {
	//se v for membro de s
	model.addConstr(expr == 1);
      } else {
	if (G.g.id(v) == G.g.id(G.t[z])) {
	  //se v for membro de t
	  model.addConstr(expr == -1);
	} else {
	  model.addConstr(expr == 0);
	}
      }
    }
  }
  model.update();
  model.optimize();
  //model.write("model.lp"); system("cat model.lp");
  if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
    cout << "Erro, sistema impossivel" << endl;
    exit(1);
  }
  cout << "Time: " <<  model.get(GRB_DoubleAttr_Runtime) << endl;
  try {
    if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
      G.BestLB = model.get(GRB_DoubleAttr_ObjBound);
      G.BestVal = G.BestLB;
      return false;
    }
  } catch (...) {
    cout << "Erro, nenhuma solução encontrada no tempo" << endl;
  }
  
  if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
      extract_solution(G,g,b);
      G.BestVal = model.get(GRB_DoubleAttr_ObjVal);
      G.BestLB = G.BestVal;
      return true;
  }
  return true;
}



int main(int argc, char *argv[]) 
{
  ListGraph g;
  int NNodes,NEdges,NPairs;
  NodeName nodename(g);
  ListGraph::NodeMap<double> posx(g), posy(g);
  vector<Node> s;
  vector<Node> t;
  vector<double> Tmax;
  vector<double> q;
  EdgeWeight custo(g);
  EdgeWeight capacidade(g);
  EdgeWeight latencia(g);
  vector<vector<Node> > BestSol;

  //double BestVal;
  //double BestLB;
	
  int seed=0;
  srand48(seed);

  // uncomment one of these lines to change default pdf reader, or insert new one
  //set_pdfreader("open");    // pdf reader for Mac OS X
  //set_pdfreader("xpdf");    // pdf reader for Linux
  //set_pdfreader("poppler");  // pdf reader for Linux
  //set_pdfreader("open -a Skim.app");
  // double cutoff;   // used to prune non promissing branches (of the B&B tree)
  if (argc!=2) {
    cout << endl << "Usage: " << argv[0] << "  <filename>" << endl << endl;
    cout << "Example:      " << argv[0] << " arq1.in" << endl << endl;
    exit(0);
  }

  string filename = argv[1];
  
  ReadListGraph3(filename, g, NNodes, NEdges, NPairs, nodename, posx, posy, latencia, capacidade, custo, s, t, Tmax, q);
  Problem_Data dt(g, nodename, posx, posy, s, t, Tmax, q, custo, capacidade, latencia);
  try {
    bool res = transmissoes(dt, 100);
    if (res) {
      cout << "Melhor resultado: " << dt.BestVal << endl;
    } else {
      cout << "Resultaso obtido: " << dt.BestVal << endl;
    }
    
  } catch (...) {
    cout << "Error during callback..." << endl;
  }
}



// Project and Analysis of Algorithms
// Fl�vio Keidi Miyazawa
// Problems with connectivity: Minimum Cost k-paths (edge disjoint)
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

bool transmissoes(Problem_Data &G, long maxtime)
{
  int seed = 0;
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.getEnv().set(GRB_IntParam_Seed, seed);
  model.set(GRB_StringAttr_ModelName, "Transmiss�o"); 
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
  // is a minimization problem
  G.NPairs = G.s.size(); //NPairs ta vindo errado!!!
  cout << "NPairs " << G.NPairs << endl;
  vector <vector <GRBVar> > x(G.NEdges, vector <GRBVar> (G.NPairs));
  // ListGraph::NodeMap<GRBVar> x(g); // variables for each node
  // ListGraph::EdgeMap<GRBVar> y(g); // variables for each edge
  // for (unsigned int i = 0; i < G.s.size(); i++) {
  //   (x[i])(G.g); // variables for each edge
  // }
  char namme[100];
  int i = 0;
  for(ListGraph::EdgeIt e(G.g); e != INVALID; ++e) {
    for (int j = 0; j < G.NPairs; j++) {
      sprintf(namme,"V_%d_%d",i,j);
      //printf("V_%d_%d ",i,j);
      double tmp  = G.custo[e]*G.q[j];
      
      //cout << "G.Q > " << G.q[j] << endl;
      x[i][j] = model.addVar(0.0, 1.0, tmp,GRB_BINARY,namme);
    }
    i++;
  }
  model.update();
  try {
    model.optimize();
    //model.write("model.lp"); system("cat model.lp");
    return 1;
  } catch (...)
    {
      cout << "Error during callback 2..." << endl;
      return(0);
    }
  return 0;
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
  vector<vector<Node>> BestSol;

  double BestVal;
  double BestLB;
	
  int seed=0;
  srand48(1);

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

    bool res = transmissoes(dt, 10000);
    cout << "Melhor resultado: " << dt.BestVal << endl;
    /*
      double soma=0.0;
      for (ListGraph::NodeIt v(g);v!=INVALID;++v) vcolor[v]=BLUE; // all nodes BLUE

      for (ListGraph::EdgeIt e(g); e!=INVALID; ++e) ecolor[e]=BLUE;

      for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
      lpvar[v] = x[v].get(GRB_DoubleAttr_X);
      if (lpvar[v] > 1.0 - EPS){ soma += weight[v]; vcolor[v]=RED;}
      }
      cout << "vCover Tree Value = " << soma << endl;

      //-----------------------------------------------------------------

      ViewListGraph(g,vname,ename,px,py,vcolor,ecolor,
      "vCover cost in graph with "+IntToString(T.nnodes)+
      " nodes: "+DoubleToString(soma));
    */
  } catch (...) {cout << "Error during callback..." << endl; }
}


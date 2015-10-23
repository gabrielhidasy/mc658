// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
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

int monitoramento_em_grafo_bipartido( ListGraph &g, NodeName &vname, ListGraph::NodeMap<double> &custo, ListGraph::NodeMap<int> &solucao)
{
  int seed=0;
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.getEnv().set(GRB_IntParam_Seed, seed);
  model.set(GRB_StringAttr_ModelName, "Monitoramento em Grafo Bipartido"); // prob. name
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // is a minimization problem
  // ------------------------------------------------------
  // Construa o modelo daqui para baixo
  // ------------------------------------------------------

  // Exemplos de como voce pode declarar variaveis indexadas nos vertices ou nas arestas.
  // Nao necessariamente voce precisa dos dois tipos
  // ListGraph::NodeMap<GRBVar> x(g); // variables for each node
  // ListGraph::EdgeMap<GRBVar> y(g); // variables for each edge
  ListGraph::NodeMap<GRBVar> x(g); // variables for each node
  ListGraph::EdgeMap<GRBVar> y(g); // variables for each edge
  int name = 0;
  char namme[100];
  for(ListGraph::NodeIt v(g); v != INVALID; ++v) {
    sprintf(namme,"PC_%s",vname[v].c_str());
    x[v] = model.addVar(0.0, 1.0, custo[v],GRB_CONTINUOUS,namme); }
  model.update();
  try {
    for(ListGraph::EdgeIt e(g); e != INVALID; ++e) {
      //Para cada aresta, um dos lados e 1
      GRBLinExpr expr;
      expr += x[g.u(e)];
      expr += x[g.v(e)];
      model.addConstr(expr >= 1);
    }
    model.update();
    // ------------------------------------------------------
    // Construa o modelo daqui para cima
    // ------------------------------------------------------
    //model.write("model.lp"); system("cat model.lp");
    model.optimize();
    for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
      if (x[v].get(GRB_DoubleAttr_X)>1-EPS) solucao[v] = 1;
      else solucao[v] = 0;
      //solucao[v] = 1;
    }
    return(1);

  } catch (...) {cout << "Error during callback..." << endl; return(0);}
}

// vCover_Instance put all relevant information in one class.
class vCover_Instance {
public:
  vCover_Instance(ListGraph &graph,NodeName &vvname,
		  NodePos &posx,NodePos &posy,NodeWeight &vweight);
  ListGraph &g;
  NodeName &vname;
  NodePos &px;
  NodePos &py;
  NodeWeight &weight;
  int nnodes;
};

vCover_Instance::vCover_Instance(ListGraph &graph, NodeName &vvname,
				 NodePos &posx,NodePos &posy,
				 NodeWeight &vweight):
  g(graph), vname(vvname), px(posx), py(posy), weight(vweight)
{
  nnodes = countNodes(g);
}

int main(int argc, char *argv[]) 
{
  ListGraph g;  // graph declaration
  string digraph_vcover_filename;
  NodeName vname(g);  // name of graph nodes
  EdgeName ename(g);  // name of graph edges
  ListGraph::NodeMap<double> px(g),py(g);  // xy-coodinates for each node
  ListGraph::NodeMap<int> vcolor(g);// color of nodes
  ListGraph::EdgeMap<int> ecolor(g); // color of edges
  ListGraph::NodeMap<int> solucao(g);  // vetor 0/1 que identifica vertices da solucao
  NodeWeight lpvar(g);    // used to obtain the contents of the LP variables
  NodeWeight weight(g);   // node weights
  vector <Node> V;
  srand48(1);

  // uncomment one of these lines to change default pdf reader, or insert new one
  //set_pdfreader("open");    // pdf reader for Mac OS X
  //set_pdfreader("xpdf");    // pdf reader for Linux
  set_pdfreader("evince");  // pdf reader for Linux
  //set_pdfreader("open -a Skim.app");
  SetGraphLabelFontSize(50);
  
  if (argc!=2) {cout<<endl<<"Usage: "<< argv[0]<<"  <digraph_vcover_filename>"<< endl << endl;
    cout << "Example:      " << argv[0] << " gr_bipartido_1.in" << endl << endl; exit(0);}

  digraph_vcover_filename = argv[1];
  ReadListGraph3(digraph_vcover_filename,g,vname,weight,px,py);
  vCover_Instance T(g,vname,px,py,weight);
  for (EdgeIt e(g); e != INVALID; ++e) ename[e] = vname[g.u(e)]+" , "+vname[g.v(e)];

  cout << "Rede com " << T.nnodes << endl << endl;
  cout << "Computadores e seus custos:" << endl << endl;
  for (NodeIt v(g);v!=INVALID;++v)
    cout << "Computador[ " << vname[v] << " ] com custo " << weight[v] << "." << endl;
  cout << endl << "Conexoes entre computadores:" << endl << endl;
  for (EdgeIt e(g);e!=INVALID;++e)
    cout << "Conexao " << ename[e] << "." << endl;
  cout << endl << endl;

  

  // for (NodeIt v(g);v!=INVALID;++v) vcolor[v]=BLACK;
  // for (EdgeIt e(g); e != INVALID; ++e) ecolor[e] = BLUE;
  // for (EdgeIt e(g); e != INVALID; ++e) ename[e] = "";
  // ViewListGraph(g,vname,ename,px,py,vcolor,ecolor,digraph_vcover_filename);
  // cout << "Pause\n";  getchar();
  
  // Generate the binary variables and the objective function
  // Add one binary variable for each edge and set its cost in the objective function

  if (monitoramento_em_grafo_bipartido(g, vname, weight, solucao)) {
    // verificacao e apresentacao da solucao obtida

    for (ListGraph::EdgeIt e(g); e!=INVALID; ++e) 
      if (!(solucao[g.u(e)] || solucao[g.v(e)])) {
	cout << "Nao foi obtida uma solucao viavel.\nConexao {"
	     << vname[g.u(e)] << "---" << vname[g.v(e)] << "} nao monitorada."
	     << endl << endl; exit(0);}

    double soma=0.0;
    for (ListGraph::NodeIt v(g);v!=INVALID;++v) vcolor[v]=BLUE; // all nodes BLUE
    for (ListGraph::EdgeIt e(g); e!=INVALID; ++e) ecolor[e]=BLUE;
    cout << "Computadores Selecionados" << endl ;
    for (ListGraph::NodeIt v(g); v!=INVALID; ++v) 
      if (solucao[v]){soma += weight[v]; vcolor[v]=RED; cout << vname[v] << endl;}

    cout << endl << "Valor da Solucao = " << soma << "." << endl;

    ViewListGraph(g,vname,ename,px,py,vcolor,ecolor,
		  "Solucao do Monitoramento com custo "+DoubleToString(soma));
    return(0);
  }else{cout << "Programa linear gerado eh inviavel." << endl;return(1);}
}


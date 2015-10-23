// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
// Problems with connectivity: Maximum matching
#define BC_EPS 0.001
#define BC_INF 1000000000000.0
#include <stdio.h>
#include <string>
#include "mygraphlib.h"
#include "myutils.h"
#include <lemon/lp.h>
#include <lemon/list_graph.h>
#include <gurobi_c++.h>
using namespace lemon;
using namespace std;
#define EPS 0.00001

int main(int argc, char *argv[]) 
{
  ListGraph g;  // graph declaration
  string graph_filename;
  NodeName vname(g);  // name of graph nodes
  EdgeName ename(g);  // name of graph nodes
  ListGraph::NodeMap<double> px(g),py(g);  // xy-coodinates for each node
  ListGraph::NodeMap<int> vcolor(g);// color of nodes
  ListGraph::EdgeMap<int> ecolor(g); // color of edges
  EdgeWeight lpvar(g);    // used to obtain the contents of the LP variables
  EdgeWeight weight(g);   // edge weights
  srand48(1);

  // uncomment one of these lines to change default pdf reader, or insert new one
  //set_pdfreader("open");    // pdf reader for Mac OS X
  //set_pdfreader("xpdf");    // pdf reader for Linux
  //set_pdfreader("evince");  // pdf reader for Linux

  // double cutoff;   // used to prune non promissing branches (of the B&B tree)
  if (argc!=2) {cout<<endl << "Usage: "<< argv[0]<<" <graphfilename>"<< endl << endl;
    cout << "Example:      " << argv[0] << " gr_att48" << endl << endl;

	exit(0);}

  graph_filename = argv[1];
  ReadListGraph(graph_filename,g,vname,weight,px,py);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE); // is a maximization problem

  /* LPI variables */
  ListGraph::EdgeMap<GRBVar> x(g); // variable for connections, 1=connected, 0=not connected
  
  GRBLinExpr expressao;
  for (ListGraph::EdgeIt e(g); e != INVALID; ++e)
    x[e] = model.addVar(0.0, 1.0, weight[e], GRB_BINARY);  
  model.update();
  
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    GRBLinExpr expr;
    for (ListGraph::IncEdgeIt e(g,v); e != INVALID; ++e) expr += x[e];
    model.addConstr(expr <= 1 );
    vcolor[v] = BLUE;
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  try {
    model.optimize();
    double soma=0.0;
    for (ListGraph::EdgeIt e(g); e!=INVALID; ++e) {
      lpvar[e] = x[e].get(GRB_DoubleAttr_X);
      ename[e] = DoubleToString(weight[e]);
      if (lpvar[e] > 1.0 - EPS) { soma += weight[e]; ecolor[e] = BLUE; }
      else ecolor[e] = NOCOLOR; }
    cout << "Maximum weighted matching = " << soma << endl;

    // Esta rotina precisa do programa neato/dot do ListGraphviz 
    ViewListGraph(g,vname,ename,px,py,vcolor,ecolor,
    "Maximum weighted matching in graph with "+IntToString(countNodes(g))+
    	    " nodes:"+DoubleToString(soma));
  } catch(GRBException e) {
    cerr << "It was not possible to solve the ILP." << endl;
    cerr << "Error code = " << e.getErrorCode() << endl;
    cerr << e.getMessage();
  }
  return 0;
}


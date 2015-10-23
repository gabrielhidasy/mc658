// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
// Problems with connectivity: Generate Triangulated Digraph 
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lemon/list_graph.h>
#define BC_EPS 0.001
#define BC_INF 1000000000000.0
#include "mygraphlib.h"
#include <string>
#include "myutils.h"
#include <lemon/concepts/digraph.h>
#include <lemon/preflow.h>
using namespace lemon;
#define EPS 0.00001


// Para compilar:
//  g++ -lemon geompack.cpp myutils.cpp  mygraphlib.cpp generate_triangulated_digraph.cpp -o out


int main(int argc, char *argv[]) 
{
  int n;
  double box_width,box_height;
  Digraph g;  // graph declaration
  DiNodeName vname(g);  // name of graph nodes
  Digraph::NodeMap<double> px(g),py(g);  // xy-coodinates for each node
  Digraph::NodeMap<int> vcolor(g);// color of nodes
  Digraph::ArcMap<int> ecolor(g); // color of edges
  ArcWeight lpvar(g);    // used to obtain the contents of the LP variables
  ArcWeight weight(g);   // edge weights
  vector <DiNode> V;
  srand48(1);

  // double cutoff;   // used to prune non promissing branches (of the B&B tree)
  if (argc!=4) {cout<<"Usage: "<< argv[0]<<"<number_of_nodes_in_graph> <box_width> <box_height>"<< endl;exit(0);}

  n = atoi(argv[1]);
  box_width = atof(argv[2]);
  box_height = atof(argv[3]);

  GenerateTriangulatedListDigraph(g,vname,px,py,weight,n,box_width,box_height);

  cout << countNodes(g) << " " << countArcs(g) << endl;
  for (DiNodeIt v(g);v!=INVALID;++v) 
    cout << vname[v] << " " << px[v] << " " << py[v] << endl;
  for (ArcIt e(g);e!=INVALID;++e) 
    cout << vname[g.source(e)] << " " << vname[g.target(e)] << " " << weight[e] << endl;
  return 0;
}


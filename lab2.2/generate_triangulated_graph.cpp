// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
// Problems with connectivity: Generate Triangulated Graph 
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
using namespace lemon;
#define EPS 0.00001


// Para compilar:
//  g++ -lemon geompack.cpp myutils.cpp  mygraphlib.cpp generate_triangulated_graph.cpp -o out


int main(int argc, char *argv[]) 
{
  int n;
  double box_width,box_height;
  ListGraph g;  // graph declaration
  NodeName vname(g);  // name of graph nodes
  ListGraph::NodeMap<double> px(g),py(g);  // xy-coodinates for each node
  ListGraph::NodeMap<int> vcolor(g);// color of nodes
  ListGraph::EdgeMap<int> ecolor(g); // color of edges
  EdgeWeight lpvar(g);    // used to obtain the contents of the LP variables
  EdgeWeight weight(g);   // edge weights
  vector <Node> V;
  srand48(1);

  // double cutoff;   // used to prune non promissing branches (of the B&B tree)
  if (argc!=4) {cout<<"Usage: "<< argv[0]<<"<number_of_nodes_in_graph> <box_width> <box_height>"<< endl;exit(0);}

  n = atoi(argv[1]);
  box_width = atof(argv[2]);
  box_height = atof(argv[3]);

  GenerateTriangulatedListGraph(g,vname,px,py,weight,n,box_width,box_height);

  cout << countNodes(g) << " " << countEdges(g) << endl;
  for (NodeIt v(g);v!=INVALID;++v) 
    cout << vname[v] << " " << px[v] << " " << py[v] << endl;
  for (EdgeIt e(g);e!=INVALID;++e) 
    cout << vname[g.u(e)] << " " << vname[g.v(e)] << " " << weight[e] << endl;
  return 0;
}


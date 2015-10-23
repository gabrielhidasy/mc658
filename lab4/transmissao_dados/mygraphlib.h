// =============================================================
#ifndef MYGRAPHLIB_DEFINE
#define MYGRAPHLIB_DEFINE

#include<float.h>
#include<fstream>
#include<iomanip>
#include<iostream>
#include<lemon/concepts/digraph.h>
#include<lemon/concepts/graph.h>
#include<lemon/concepts/maps.h>
#include<lemon/list_graph.h>
#include<lemon/math.h>
#include<lemon/preflow.h>
#include<string>
#include "myutils.h"
#include "geompack.hpp"

using namespace std;
using namespace lemon;

/* read a graph from a text file (using edge format / ListGraph or ListDigraph)*/

#define MAXLABELNAME 200
#define MAXLINE 1000

#define VIEW_DOT 0
//#define VIEW_NEATO 1

typedef ListGraph::Edge Edge;
typedef ListGraph::EdgeIt EdgeIt;
typedef ListGraph::EdgeMap<double> EdgeWeight;
typedef ListGraph::EdgeMap<int> EdgeColor;
typedef ListGraph::EdgeMap<int> EdgeIndex;
typedef ListGraph::EdgeMap<string> EdgeName;
typedef ListGraph::Node Node;
typedef ListGraph::NodeIt NodeIt;
typedef ListGraph::NodeMap<bool> CutMap;
typedef ListGraph::NodeMap<double> NodeWeight;
typedef ListGraph::NodeMap<int> NodeColor;
typedef ListGraph::NodeMap<int> NodeIndex;
typedef ListGraph::NodeMap<string> NodeName;
typedef ListGraph::NodeMap<double> NodePos;
typedef Preflow<ListGraph, EdgeWeight> PFType;


typedef ListDigraph Digraph;
typedef Digraph::Node DiNode;
typedef Digraph::NodeIt DiNodeIt;
typedef Digraph::Arc Arc;
typedef Digraph::ArcIt ArcIt;
typedef Digraph::ArcMap<double> ArcWeight;
typedef Digraph::NodeMap<string> DiNodeName;
typedef Digraph::ArcMap<string> ArcName;
typedef Digraph::NodeMap<int> DiNodeColor;
typedef Digraph::ArcMap<int> ArcColor;
typedef Digraph::NodeMap<bool> DiCutMap;
typedef Digraph::NodeMap<double> DiNodePos;



// read a list digraph. If go_and_back is true, for a line [u,v,cost] the 
// routine insert the arc (u,v) and (v,u), both with cost custo. Otherwise,
// it insert only the arc (u,v).

bool ReadListDigraph(string filename,
		     ListDigraph &g,
		     ListDigraph::NodeMap<string>& nodename,
		     ListDigraph::ArcMap<double>& custo,
		     const bool go_and_back);

bool ReadListDigraph(string filename,
		     ListDigraph &g,
		     ListDigraph::NodeMap<string>& vname,
		     ListDigraph::ArcMap<double>& weight,
		     ListDigraph::NodeMap<double>& posx,
		     ListDigraph::NodeMap<double>& posy,
		     const bool go_and_back);

// Read a geometric graph (points in the euclidean plane) or a list graph
// If the graph is geometric, the positions (posx and posy) are the given points. 
// If the graph is a list graph, the positions are computed by the neato program.
bool ReadListGraph(string filename,
		   ListGraph &g,
		   ListGraph::NodeMap<string>& nodename,
		   ListGraph::EdgeMap<double>& weight,
		   ListGraph::NodeMap<double>& posx,
		   ListGraph::NodeMap<double>& posy);

bool ReadListGraph3(string filename,
		   ListGraph &g,
		   int& NNodes,
		   int& NEdges,
		   int& NPairs,
		   ListGraph::NodeMap<string>& nodename,
		   ListGraph::NodeMap<double>& posx,
		   ListGraph::NodeMap<double>& posy,
		   ListGraph::EdgeMap<double>& latencia,
		   ListGraph::EdgeMap<double>& capacidade,
		   ListGraph::EdgeMap<double>& custo,
		   vector<Node>& s,
		   vector<Node>& t,
		   vector<double>& Tmax,
		   vector<double>& q);

// ==============================================================



bool GenerateVertexPositions(ListGraph &g,
			     ListGraph::EdgeMap<double>& custo,
			     ListGraph::NodeMap<double>& posx,
			     ListGraph::NodeMap<double>& posy);

void PrintListGraph(ListGraph &g, ListGraph::NodeMap<string> &vname, EdgeWeight &graphweight);


// This routine visualize a graph using a pdf viewer. It uses neato (from
// graphviz.org) to generate a pdf file and a program to view the pdf file. The
// pdf viewer name is given in the viewername parameter.
int ViewListGraph(ListGraph &g,
	      ListGraph::NodeMap<string> &vname, // name of the nodes
	      ListGraph::EdgeMap<string> &ename,  // name of edges
	      ListGraph::NodeMap<double>& px, // x-position of the nodes
	      ListGraph::NodeMap<double>& py, // y-position of the nodes
	      ListGraph::NodeMap<int>&    vcolor, // color of node (see myutils.h)
	      ListGraph::EdgeMap<int>&    ecolor, // color of edge
		  //string viewername, // name of the program that opens a pdf file
		  string text); // text displayed below the figure

int ViewEuclideanGraph(ListGraph &g,
		       ListGraph::NodeMap<string> &vname, // nome dos vertices
		       ListGraph::NodeMap<double> &posx, // coord. x dos vertices
		       ListGraph::NodeMap<double> &posy, // coord. y dos vertices
		       ListGraph::NodeMap<int> &vcolor,  // cor dos vertices
		       ListGraph::EdgeMap<int> &acolor);  // cor das arestas


// This routine visualize a digraph using a pdf viewer. It uses neato (from
// graphviz.org) to generate a pdf file and a program to view the pdf file. The
// pdf viewer name is given in the viewername parameter.
int ViewListDigraph(ListDigraph &g,
      ListDigraph::NodeMap<string> &vname, // node names
      ListDigraph::NodeMap<double>& px, // x-position of the nodes
      ListDigraph::NodeMap<double>& py, // y-position of the nodes
      ListDigraph::NodeMap<int>&    vcolor, // color of node (see myutils.h)
      ListDigraph::ArcMap<int>&    ecolor, // color of edge 
      string text); // text displayed below the figure


// Generate a triangulated ListDigraph, building the Delaunay
// triangulation of random points. Each edge of the Delaunay triangulation
// leads to two arcs (in both senses)  
// Uses the geompack program, available in
// http://people.sc.fsu.edu/~jburkardt/cpp_src/geompack/geompack.html
bool GenerateTriangulatedListDigraph(ListDigraph &g,
		 DiNodeName &vname, // name of the nodes
		 ListDigraph::NodeMap<double>& px, // x-position of the nodes
		 ListDigraph::NodeMap<double>& py, // y-position of the nodes
		 ListDigraph::ArcMap<double>& weight, // weight of edges
		 int n, // number of nodes
		 double SizeX, // coordinate x is a random number in [0,SizeX)
                 double SizeY); // coordinate y is a random number in [0,SizeY)

// the same as above, but for non-oriented edges
bool GenerateTriangulatedListGraph(ListGraph &g, // return with generated graph
			  ListGraph::NodeMap<string> &vname, // return with name of the nodes
			  ListGraph::NodeMap<double>& px, // return with x-position of the nodes
			  ListGraph::NodeMap<double>& py, // return with y-position of the nodes
			  ListGraph::EdgeMap<double>& weight, // return with weight of edges
			  int n, // number of nodes 
			  double SizeX, // coordinate x is a random number in [0,SizeX)
			  double SizeY); // coordinate y is a random number in [0,SizeY)

/* // Save a Digraph in a text file */
/* bool SaveListDigraph(ListDigraph &g, */
/* 		     DiNodeName &vname, // name of the nodes */
/* 		     ListDigraph::NodeMap<double>& px, // x-position of the nodes */
/* 		     ListDigraph::NodeMap<double>& py, // y-position of the nodes */
/* 		     ListDigraph::ArcMap<double>& weight, // weight of edges */
/* 		     string outputfilename); // name of the destination text file */

// =============================================================

class AdjacencyMatrix {
public:
  AdjacencyMatrix(ListGraph &graph,EdgeWeight &graphweight,double NonEdgeValue);
  ~AdjacencyMatrix();
  double *AdjMatrix;
  ListGraph *g;
  EdgeWeight *weight;
  int Nnodes,Nedges,Nmatrix;
  double NonEdgeValue;
  Node *Index2Node;
  Edge *Index2Edge;
  double Cost(Node,Node);
  double Cost(Edge);
  NodeIndex Node2Index;
  EdgeIndex Edge2Index;
};

//Generate a random complete euclidean ListGraph
bool GenerateRandomEuclideanListGraph(ListGraph &g,
		  ListGraph::NodeMap<string> &vname, // node names
		  ListGraph::NodeMap<double>& px, // x-position of the nodes
		  ListGraph::NodeMap<double>& py, // y-position of the nodes
		  ListGraph::EdgeMap<double>& weight, // weight of edges
		  int n, // number of nodes
		  double SizeX, // coordinate x is a random number in [0,SizeX)
  	          double SizeY); // coordinate y is a random number in [0,SizeY)

//Generate a random complete euclidean ListDigraph
bool GenerateRandomEuclideanListDigraph(ListDigraph &g,
			  ListDigraph::NodeMap<string> &vname, // node name
			  ListDigraph::NodeMap<double>& px, // x-position of the node
			  ListDigraph::NodeMap<double>& py, // y-position of the node
			  ListDigraph::ArcMap<double>& weight, // weight of edges
			  int n, // number of nodes
			  double SizeX, // coordinate x is a random number in [0,SizeX)
					double SizeY); // coordinate y is a random number in [0,SizeY)


// Given a color code, return its name
//std::string ColorName(int cor);

// Obtain a mininum cut for directed graphs from s to t.
// The returned cut is given by the vector of nodes 'cut' (boolean
// vector: nodes v in the same side of s have cut[v]=true, otherwise cut[v]=false.
double DiMinCut(ListDigraph &g, ArcWeight &weight, DiNode &s,DiNode &t, DiCutMap &cut);


#endif

  



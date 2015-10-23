//----------------------------------------------------------------------
// Example of an exact program to solve the Minimum Traveling
// Salesman Problem, using LEMON Library and Integer Linear Programming
// with GUROBI Solver.
//
// This program finds a minimum Traveling Salesman Tour (TSP) via a branch
// and cut approach. The formulation used is given in the page 100 of the
// slides in the link (in portuguese)
//
// http://www.ic.unicamp.br/~fkm/lectures/proglin.pdf
//
// To obtain the cuts that are violated, it uses the Gomory-Hu subroutine, 
// available from LEMON package. For a short explanation why it is interesting
// to use Gomory-Hu tree, see pages 141-142 of the above slides (in portuguese).
//
// OBS.: The edges cost in the graphs in the same directory does not have the
// same costs computed by TSPLIB
//
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------
#include <gurobi_c++.h>
#include <float.h>
#include <math.h>
#include <lemon/list_graph.h>
#include <lemon/gomory_hu.h>
#define BC_EPS 0.001
#define BC_INF 1000000000000.0
#include "mygraphlib.h"
#include "myutils.h"

// there are some differences from Gurobi Version 5.5, see below
#define GUROBI_NEWVERSION 1  

// This is the type used to obtain the pointer to the problem data. This pointer
// is stored in the branch and cut tree. And when we define separation routines,
// we can recover the pointer and access the problem data again.
class TSP_Data {
public:
  TSP_Data(ListGraph &graph,NodeName &nodename,NodePos &posicaox,
	   NodePos &posy,EdgeWeight &eweight);
  ListGraph &g;
  int NNodes,NEdges;
  int max_perturb2opt_it; // maximum number of iterations for heuristic Perturb2OPT
  NodeName &vname;
  EdgeName ename;
  NodeColor vcolor;
  EdgeColor ecolor;
  EdgeWeight &weight;
  NodePos &posx;
  NodePos &posy;
  AdjacencyMatrix AdjMat; // adjacency matrix
  vector<Node> BestCircuit; // vector containing the best circuit found
  double BestCircuitValue;
};

TSP_Data::TSP_Data(ListGraph &graph,NodeName &nodename,NodePos &posicaox,
		   NodePos &posicaoy,EdgeWeight &eweight):
  g(graph),
  vname(nodename),
  ename(graph),
  vcolor(graph),
  ecolor(graph),
  weight(eweight),
  posx(posicaox),
  posy(posicaoy),
  AdjMat(graph,eweight,BC_INF), //
  BestCircuit(countEdges(graph)) 
{
  NNodes=countNodes(this->g);
  NEdges=countEdges(this->g);
  BestCircuitValue = DBL_MAX;
  max_perturb2opt_it = 3000; // default value
}


bool Heuristic_2_OPT(AdjacencyMatrix &A,vector<Node> &Circuit,double &BestCircuitValue, int &NNodesCircuit)
{
  double CurrentWeight=0.0,Remove,Insert;
  bool globalimproved,improved;
  vector<Node> CircuitAux(NNodesCircuit);
  int i,j,k,l;
  for (int i=0;i<NNodesCircuit-1;i++) 
    CurrentWeight += A.Cost(Circuit[i],Circuit[i+1]);
  CurrentWeight += A.Cost(Circuit[NNodesCircuit-1],Circuit[0]);

  i=0;
  globalimproved = false;
  do {
    int i1,i2,j1,j2;
    improved=false;
    i1 = i;
    do {
      // one edge is (i1 , i1+1)
      i2 = (i1+1)%NNodesCircuit;
      for (j=2;j<NNodesCircuit-1;j++) {
	j1 = (i1+j)%NNodesCircuit;
	j2 = (j1+1)%NNodesCircuit;
	Remove = A.Cost(Circuit[i1],Circuit[i2])+A.Cost(Circuit[j1],Circuit[j2]);
	Insert = A.Cost(Circuit[i1],Circuit[j1])+A.Cost(Circuit[i2],Circuit[j2]);
	if (Remove-Insert > 0) {
	  k = 0;
	  CircuitAux[k++] = Circuit[i1];
	  for (l=j1;l!=i1;l=(l-1+NNodesCircuit)%NNodesCircuit)
	    CircuitAux[k++] = Circuit[l];
	  for (l=j2;l!=i1;l=(l+1)%NNodesCircuit)
	    CircuitAux[k++] = Circuit[l];

	  for (k=0;k<NNodesCircuit;k++) Circuit[k] = CircuitAux[k];
	  CurrentWeight = CurrentWeight - Remove + Insert;
	  if (CurrentWeight < BestCircuitValue-BC_EPS) {
	    BestCircuitValue = CurrentWeight;
	    globalimproved = true;
	  }
	  improved = true; 
	}
      }
      i1=(i1+1)%NNodesCircuit;
    } while (i1!=i);
    i = (i+1)%NNodesCircuit;
  }while (improved);
  if (globalimproved) cout << "[Heuristic: 2OPT] New Solution of value " << BestCircuitValue << "\n";
  return(globalimproved);
}

// This routine must be called when the vector x (indexed on the edges) is integer.
// The contained circuit is transformed into a circuit represented by a sequence of nodes.
// The vector x must represent a circuit (must be conected)
bool Update_Circuit(TSP_Data &tsp, ListGraph::EdgeMap<GRBVar>& x)
{
  ListGraph::NodeMap<Node> Adj1(tsp.g), Adj2(tsp.g);
  double CircuitValue,f;
  Node FirstNode=INVALID,u,ant;
  int n,i,NNodesCircuit;
  
  NNodesCircuit=0; // used to count nodes, only to verify correctness
  // Adj1[v] and Adj2[v] have the two adjacent nodes of v in the circuit
  for (NodeIt v(tsp.g); v != INVALID; ++v) { Adj1[v]=INVALID; Adj2[v]=INVALID; }
  CircuitValue = 0.0;
  for (EdgeIt e(tsp.g); e != INVALID; ++e) { 
    f = x[e].get(GRB_DoubleAttr_X);
    assert(!(IsFrac(f))); //f cannot be fractional;}
    if (f > 1.0-BC_EPS) {  // if the edge is in the solution
      Node u,v;
      NNodesCircuit++;
      u = (tsp.g.u(e)); v = (tsp.g.v(e));   // then, obtain the edge nodes u and v
      CircuitValue += tsp.weight[e];
      assert((Adj1[u]==INVALID)||(Adj2[u]==INVALID));
      assert((Adj1[v]==INVALID)||(Adj2[v]==INVALID));  // and say that
      if (Adj1[u]==INVALID) Adj1[u] = v; else Adj2[u] = v; // v is adjacent to u
      if (Adj1[v]==INVALID) Adj1[v] = u; else Adj2[v] = u; // and u is adj. to v
    }
  }
  //cout << "CircuitValue = " << CircuitValue << " BestPrevious = "<< tsp.BestCircuitValue << endl;
  if (CircuitValue > tsp.BestCircuitValue-BC_EPS) return(false);
  // First verify if the circuit is valid
  assert (NNodesCircuit==tsp.NNodes);
  for (NodeIt v(tsp.g); v != INVALID; ++v) {
    assert((Adj1[v]!=INVALID) && (Adj2[v]!=INVALID));
    // cout << "V_" << tsp.vname[v];
    // if (Adj1[v]!=INVALID) cout << " Adj1[ V_" << tsp.vname[v] << " ] = "<< tsp.vname[Adj1[v]];
    // else cout << " Adj1[ V_" << tsp.vname[v] << " ] = INV";
    // if (Adj2[v]!=INVALID) cout << " Adj2[ V_" << tsp.vname[v] << " ] = "<< tsp.vname[Adj2[v]];
    // else cout << " Adj2[ V_" << tsp.vname[v] << " ] = INV";
    // cout << endl;
  }

  // find one node of the circuit
  for (NodeIt v(tsp.g); v != INVALID; ++v) {FirstNode = v; break;}
  
  n = 1;     ant = FirstNode;     u = Adj1[FirstNode];
  // verify if each node have two adjacent nodes
  while (u!=FirstNode) {
    assert(u!=INVALID);
    if (Adj1[u]==ant) {ant = u; u=Adj2[u];} else {assert(Adj2[u]==ant); ant=u; u=Adj1[u]; }
    n++;
  }
  if (n!=tsp.NNodes) cout << "Solution is not connected" << endl;
  assert(n==tsp.NNodes); 

  // the circuit is ok. So, update the best circuit (remember that at this point
  tsp.BestCircuitValue = CircuitValue; // the new circuit is better than the 
  i = 1;                               // previous best circuit
  tsp.BestCircuit[0] = FirstNode;
  u = Adj2[FirstNode];
  ant = FirstNode;
  while (u!=FirstNode) {
    tsp.BestCircuit[i] = u;      i++;
    if (Adj1[u]==ant) {ant=u; u=Adj2[u];} else {ant=u; u=Adj1[u];}
  }
  Heuristic_2_OPT(tsp.AdjMat,tsp.BestCircuit,tsp.BestCircuitValue,tsp.NNodes);
  return(true);
}

class subtourelim: public GRBCallback
{
  TSP_Data &tsp;
  ListGraph::EdgeMap<GRBVar>& x;

  double (GRBCallback::*solution_value)(GRBVar);
public:
  subtourelim(TSP_Data &tsp, ListGraph::EdgeMap<GRBVar>& x) : tsp(tsp),x(x)
  {    }
protected:
  void callback()
  {
    if  (where==GRB_CB_MIPSOL) // if mip solution (in this case, integer solution, ***but can be unconnected***)
      {solution_value = &subtourelim::getSolution;}  // pointer to function that obtain solution
    else // it is not a candidate solution
      // Exercise: try to comment the next two lines and test the execution time
      if ((where==GRB_CB_MIPNODE) && (getIntInfo(GRB_CB_MIPNODE_STATUS)==GRB_OPTIMAL))// if node with optimal fractional solution
	{solution_value = &subtourelim::getNodeRel;} else // pointer to function that obtain node (fractional) solution
	return;
    try {
      typedef ListGraph::EdgeMap<double> capacityType;
      capacityType capacity(tsp.g);
      for (ListGraph::EdgeIt e(tsp.g); e != INVALID; ++e) {
	capacity[e] = (this->*solution_value)(x[e]);
      }
      GomoryHu<ListGraph, capacityType> ght(tsp.g, capacity);
      ght.run();

      // The Gomory-Hu tree is given as a rooted directed tree. Each node has
      // an arc that points to its father. The root node has father -1.
      // Remember that each arc in this tree represents a cut and the value of
      // the arc is the weight of the corresponding cut. So, if an arc has weight
      // less than 2 than we found a violated cut and in this case, we insert the
      // corresponding constraint.
      // This may be a good approach in the beggining of the execution, but as the 
      // current fractional solution becomes "more connected", few cuts are violated
      // and GH tree becomes a heavy approach.
      for (NodeIt u(tsp.g); u != INVALID; ++u) {
	GRBLinExpr expr = 0;
	double vcut;
	if ((tsp.g).id(ght.predNode(u))==-1) continue;
	vcut =  ght.predValue(u);
	if (vcut > 2.0 - BC_EPS) continue;
	for(GomoryHu<ListGraph,EdgeWeight>::MinCutEdgeIt a(ght,u,ght.predNode(u));a!=INVALID;++a) 
	  expr += x[a];
	addLazy( expr >= 2 );
      }
    } catch (...) {
      cout << "Error during callback..." << endl;
    }
  }
};


void ChangeNode(Node &a,Node &b)
{ Node c;  c = a; a = b; b = c; }
  
// This routine starts with some solution (current best or a random generated
// solution) and iteratively perturb changing some nodes and applying 2OPT.
bool TSP_Perturb2OPT(TSP_Data &tsp)
{
  ListGraph *g;
  int nchanges,i,j;
  g = &tsp.g;
  vector<Node> Circuit(tsp.NNodes), BestCircuit(tsp.NNodes);
  double BestCircuitValue;
  //nchanges = 2 (tests with some instances indicate that nchanges=1 is better
  nchanges = 1;

  // Start with a initial solution (if there is no solution, generate any sequence)
  if (tsp.BestCircuitValue < DBL_MAX) {
    BestCircuitValue = tsp.BestCircuitValue;
    for (int i=0; i<tsp.NNodes; i++) BestCircuit[i] = tsp.BestCircuit[i];
  } else {
    int i=0;
    BestCircuitValue = 0.0;
    for (ListGraph::NodeIt v(*g); v!=INVALID; ++v) BestCircuit[i++]=v;
    for (i=0;i<tsp.NNodes;i++) 
      BestCircuitValue += 
	tsp.AdjMat.Cost(BestCircuit[i] , BestCircuit[(i+1)%tsp.NNodes]);
  }

  for (int it=0;it<tsp.max_perturb2opt_it;it++) {
    if (!(it%100)) printf("[Heuristic: Perturbation+2OPT] it = %d (of %d)\n",it+1,tsp.max_perturb2opt_it);
    for (int k=0;k<tsp.NNodes;k++) Circuit[k] = BestCircuit[k];
    for (int nc=0;nc<nchanges;nc++) {
      i = (int) (drand48()*tsp.NNodes);  // get two random nodes and exchange
      j = (int) (drand48()*tsp.NNodes);  // their positions
      if (i!=j) ChangeNode(Circuit[i],Circuit[j]);
    }
    Heuristic_2_OPT(tsp.AdjMat,Circuit,BestCircuitValue,tsp.NNodes);
    if (BestCircuitValue < tsp.BestCircuitValue) { //update the best circuit used
      tsp.BestCircuitValue = BestCircuitValue;     // by the heuristic
      for (int i=0;i<tsp.NNodes;i++) {
	BestCircuit[i] = Circuit[i];
	tsp.BestCircuit[i] = Circuit[i];
      }
    }
  }
  return(true);
}

// Convert a double value to string
string DoubleToString2(double x)
{ std::ostringstream oss; oss << x; return(oss.str()); }

void ViewTspCircuit(TSP_Data &tsp)
{
  ListGraph h;
  ListGraph::NodeMap<string> h_vname(h);  // node names
  ListGraph::NodeMap<Node> h_g2h(tsp.g);  // maps a node of g to a node of h
  ListGraph::NodeMap<double> h_posx(h);
  ListGraph::NodeMap<double> h_posy(h);
  ListGraph::NodeMap<int> vcolor(h);   // color of the vertices
  ListGraph::EdgeMap<int> acolor(h);  // color of edges
  ListGraph::EdgeMap<string> aname(h);  // name of edges
  for (ListGraph::NodeIt v(tsp.g); v!=INVALID; ++v) {
    Node hv;
    hv = h.addNode();
    h_g2h[v] = hv;
    h_posx[hv] = tsp.posx[v];
    h_posy[hv] = tsp.posy[v];
    h_vname[hv] = tsp.vname[v];
    vcolor[hv] = BLUE;
  }
  for (int i=0;i<tsp.NNodes;i++) {
    ListGraph::Node u,v;
    ListGraph::Edge a;
    u = tsp.BestCircuit[i]; 
    v = tsp.BestCircuit[(i+1) % tsp.NNodes]; 
    a = h.addEdge(h_g2h[u] , h_g2h[v]);
    aname[a] = "";
    acolor[a] = BLUE;
  }
  ViewListGraph(h,h_vname,aname,h_posx,h_posy,vcolor,acolor,"TSP Circuit with cost "+DoubleToString(tsp.BestCircuitValue));
}


// Only to see if a file exists. It (tries to) open and close the file.
//
// bool FileExists(const char *filename)
// {
//   bool r;
//   ifstream ifile(filename);
//   r = ifile;
//   ifile.close();
//   return r;
// }
//
bool FileExists(const string& filename)
{
    ifstream f(filename.c_str());
    return f.is_open();
}


int main(int argc, char *argv[]) 
{
  int time_limit;
  char name[1000];
  double cutoff=0.0;
  ListGraph g;
  EdgeWeight lpvar(g);
  EdgeWeight weight(g);
  NodeName vname(g);
  ListGraph::NodeMap<double> posx(g),posy(g);
  string filename;

  int seed=1;


  // uncomment one of these lines to change default pdf reader, or insert new one
  //set_pdfreader("open");    // pdf reader for Mac OS X
  //set_pdfreader("xpdf");    // pdf reader for Linux
  //set_pdfreader("evince");  // pdf reader for Linux

  srand48(seed);
  time_limit = 3600; // solution must be obtained within time_limit seconds
  if (argc!=2) {cout<< endl << "Usage: "<< argv[0]<<" <graph_filename>"<<endl << endl <<
      "Example: " << argv[0] << " gr_berlin52" << endl <<
      "         " << argv[0] << " gr_att48" << endl << endl; exit(0);}
  
  else if (!FileExists(argv[1])) {cout<<"File "<<argv[1]<<" does not exist."<<endl; exit(0);}
  filename = argv[1];
  
  // Read the graph
  if (!ReadListGraph(filename,g,vname,weight,posx,posy)) 
    {cout<<"Error reading graph file "<<argv[1]<<"."<<endl;exit(0);}

  TSP_Data tsp(g,vname,posx,posy,weight); 
  ListGraph::EdgeMap<GRBVar> x(g);
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
#if GUROBI_NEWVERSION
  model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
  model.getEnv().set(GRB_IntParam_Seed, seed);
#else
  model.getEnv().set(GRB_IntParam_DualReductions, 0); // Dual reductions must be disabled when using lazy constraints
#endif
  model.set(GRB_StringAttr_ModelName, "Undirected TSP with GUROBI"); // name to the problem
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // is a minimization problem
  
  // Add one binary variable for each edge and also sets its cost in the objective function
  for (ListGraph::EdgeIt e(g); e!=INVALID; ++e) {
    sprintf(name,"x_%s_%s",vname[g.u(e)].c_str(),vname[g.v(e)].c_str());
    x[e] = model.addVar(0.0, 1.0, weight[e],GRB_BINARY,name);
  }
  model.update(); // run update to use model inserted variables

  // Add degree constraint for each node (sum of solution edges incident to a node is 2)
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    GRBLinExpr expr;
    for (ListGraph::IncEdgeIt e(g,v); e!=INVALID; ++e) expr += x[e];
    //aqui model.addConstr(expr == 2 ); what? ignorou!
    model.addConstr(expr == 2 );

  }

  try {
    model.update(); // Process any pending model modifications.
    if (time_limit >= 0) model.getEnv().set(GRB_DoubleParam_TimeLimit,time_limit);

    subtourelim cb = subtourelim(tsp , x);
    model.setCallback(&cb);
    
    tsp.max_perturb2opt_it = 200; //1000; // number of iterations used in heuristic TSP_Perturb2OPT
    TSP_Perturb2OPT(tsp);
    if (tsp.BestCircuitValue < DBL_MAX) cutoff = tsp.BestCircuitValue-BC_EPS; // 
    // optimum value for gr_a280=2579, gr_xqf131=566.422, gr_drilling198=15808.652
    if (cutoff > 0) model.getEnv().set(GRB_DoubleParam_Cutoff, cutoff );
    model.update(); // Process any pending model modifications.
    model.optimize();

    double soma=0.0;
    for (ListGraph::EdgeIt e(g); e!=INVALID; ++e) {
      lpvar[e] = x[e].get(GRB_DoubleAttr_X);
      if (lpvar[e] > 1-BC_EPS ) {
	soma += weight[e];
	if (
	    (vname[g.u(e)] == "243")||(vname[g.v(e)] == "243") ||
	    (vname[g.u(e)] == "242")||(vname[g.v(e)] == "242")
	    ) {
	  cout << "Achei, x("<< vname[g.u(e)] << " , " << vname[g.v(e)] << " = " << lpvar[e] <<"\n";
	}
      }
    }

    cout << "Solution cost = "<< soma << endl;
    Update_Circuit(tsp,x); // Update the circuit in x to tsp circuit variable (if better)
    ViewTspCircuit(tsp);

  }catch (...) {
    if (tsp.BestCircuitValue < DBL_MAX) {
      cout << "Heuristic obtained optimum solution"  << endl;
      ViewTspCircuit(tsp);
      return 0;
    }else {
      cout << "Graph is infeasible"  << endl;
      return 1;
    }
  }
}

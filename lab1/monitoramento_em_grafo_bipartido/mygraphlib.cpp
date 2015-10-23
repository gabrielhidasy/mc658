#include <math.h>
#include "myutils.h"
#include <iostream> 
#include <fstream> 
#include <cstdlib>
#include <cstring>
#include<lemon/math.h>
#include "mygraphlib.h"
#include <tr1/unordered_map>

using namespace std;


int graphlabelfontsize=10;
int graphnodefontsize=50;
int graphedgefontsize=50;

void SetGraphLabelFontSize(int s)
{
  graphlabelfontsize = s;
}


void PulaBrancoComentario(ifstream &ifile)
{ char c;  string line;
  while (!ifile.eof()) {
    c = ifile.get();
    while ((c==' ') && (!ifile.eof())) c = ifile.get();
    if (!ifile.eof()) { if (c=='#')  getline(ifile,line);  else {ifile.unget(); break; } }
  }
}

bool WriteListGraphGraphviz(ListGraph &g,
		   ListGraph::NodeMap<string> &vname, // vertex names
		   ListGraph::EdgeMap<string> &ename,  // edge names
		   ListGraph::NodeMap<int> &vcolor,   // vertex colors
		   ListGraph::EdgeMap<int> &acolor,    // edge colors
		   string filename)
{
  ofstream out;
  string linha;
  
  out.open(filename.c_str());
  if (out.is_open()) return(false);

  out << "graph g {\n";
  out << "\tsize = \"" << graphlabelfontsize << ", " << graphlabelfontsize << "\";\n";
  out << "\tnode [shape = \"circle\", fontsize = \"" << graphnodefontsize << "\"];\n";
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    linha = "\t";    linha += vname[v].c_str();   linha += " [color=";
    linha += ColorName(vcolor[v]);                linha += "];\n";
    out << linha;
  }
  for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
    if (acolor[a]!=WHITE) {
      linha = "\t";
      linha += vname[g.u(a)].c_str() ;
      linha += "  -- ";
      linha += vname[g.v(a)].c_str();
      linha += " [label = \"";  linha += ename[a].c_str();  linha += "\"";
      linha += ", color=\""; linha += ColorName(acolor[a]); linha += "\" ];\n";
      out << linha;
    }
  }
  out << "}\n";
  out.close();
  return(true);
}




/* Dado um texto e um padrao, devolve a posicao no
   texto onde o padrao ocorre. Devolve (-1) se nao achar. */
int gr_busca_padrao(char *texto,char *busca)
{
  int i,j,tam_texto,tam_busca;
  char *p1,*p2;
  bool achou;
  tam_texto = strlen(texto);
  tam_busca = strlen(busca);
  for (i=0;i<tam_texto - tam_busca + 1;i++) {
      p1 = &texto[i];
      p2 = busca;
      for (j=0,achou = true;j<tam_busca && achou == true;j++)
	if (p1[j] != p2[j]) achou = false;
      if (achou) return(i);
  }
  return(-1);
}


bool ReadEuclideanListGraph(string filename,
			    ListGraph &g,
			    ListGraph::NodeMap<string>& vname,
			    ListGraph::EdgeMap<double>& custo,
			    ListGraph::NodeMap<double>& posx,
			    ListGraph::NodeMap<double>& posy)
{
  int i,n,m;
  ListGraph::Node nu,nv;
  ListGraph::Edge a;
  char nomev[100];
  ListGraph::Node v;
  double px,py;
  
  ifstream ifile;  ifile.open(filename.c_str());  if (!ifile) return(false);
  PulaBrancoComentario(ifile);
  // format: <number_of_nodes>   -1
  // The value -1 is to indicate that there is no edge/arc, as edge weights
  // are given by the euclidean distance
  ifile >> n;    ifile >> m;
  if (m!=-1) {
    printf("Wrong format in the euclidean graph of file %s.\n",filename.c_str());
    return(false);
  }
  for (i=0;i<n;i++) {
    ifile >> nomev;  ifile >> px; ifile >> py;
    v = g.addNode();    vname[v] = nomev;    posx[v] = px;    posy[v] = py;
  }

  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    ListGraph::NodeIt u(g);
    u=v;
    for (++u; u!=INVALID; ++u) {
      a = g.addEdge(u,v);
      custo[a] = sqrt((posx[u]-posx[v])*(posx[u]-posx[v]) + 
		      (posy[u]-posy[v])*(posy[u]-posy[v]));
    }
  }
  ifile.close();
  return(true);
}

// This routine use the neato program to generate positions.
bool GenerateVertexPositions(ListGraph &g,
			     ListGraph::EdgeMap<double>& custo,
			     ListGraph::NodeMap<double>& posx,
			     ListGraph::NodeMap<double>& posy)
{
  size_t t=0;
  double x,y;
  char tempname[1000],tempnamedot[1000],tempnameposdot[1000],cmd[1000];
  ofstream out;
  ifstream in;
  string linha,substring;

  (void) custo;// to avoid "non-used" parameter message.
  
  // obtain a temporary file name
  strcpy(tempname,".readgraphtempname");
  strcpy(tempnamedot,tempname);   strcat(tempnamedot,".dot");
  strcpy(tempnameposdot,tempname);strcat(tempnameposdot,"_pos.dot");
  
  out.open(tempnamedot);
  if (!out.is_open()) return(false);

  out << "graph g {\n";
  out << "\tsize = \"" << graphlabelfontsize << ", " << graphlabelfontsize << "\";\n";
  out << "\tnode [shape = \"circle\", fontsize = \"" << graphnodefontsize << "\"];\n";
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    linha = "\t";    linha += IntToString(g.id(v));   linha += ";\n";
    out << linha;
  }
  for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
    linha = "\t";  linha += IntToString(g.id(g.u(a)));
    linha += "  -- ";  linha += IntToString(g.id(g.v(a))); linha += ";\n";
    out << linha;
  }
  out << "}\n";
  out.close();
  sprintf(cmd,"neato -Goverlap=false %s -o %s",tempnamedot,tempnameposdot); 
  //printf("neato -Goverlap=false %s -o %s\n",tempnamedot,tempnameposdot);
  fflush(stdout);
  system(cmd); // gera outro arquivo do neato, mas com posições
  
  in.open(tempnameposdot);
  if (!in.is_open()) return(false);
  while (!in.eof()) {
    getline(in,linha);
    t = linha.find("{");
    if (t!=string::npos) break;
  }
  if (t==string::npos) {
    cout << "Temp Graphviz file is not appropriate for GenerateVertexPositions.\n";
    exit(0);
  }
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    getline(in,linha);
    // avoid info on the graph, node or edges
    while ((!in.eof()) && ((linha.find("graph [")!=string::npos) ||
			   (linha.find("node [")!=string::npos) ||
			   (linha.find("edge [")!=string::npos) || 
			   (linha.find(" -- ")!=string::npos))) {
      while ((!in.eof()) && (linha.find("];")==string::npos)) {
	string linha2;
	getline(in,linha2);
	linha += linha2;
      }
      getline(in,linha);
      while ((!in.eof()) && (linha.find("[")==string::npos)) {
	string linha2;
	getline(in,linha2);
	linha += linha2;
      }
    }
    if (linha.find("[")!=string::npos) {
      while (linha.find("];")==string::npos) {
	string linha2;
	getline(in,linha2);
	linha += linha2;
      }
    }
    t = linha.find("pos=\"");
    if (t!=string::npos) {
      stringstream s;
      int nodeid;
      sscanf(linha.c_str(),"%d",&nodeid);
      substring = linha.substr(t+5);
      sscanf(substring.c_str(),"%lf,%lf",&x,&y);
      //printf("%lf , %lf",x,y);
      for (ListGraph::NodeIt vv(g); vv!=INVALID; ++vv) {
	if (nodeid==g.id(vv)){
	  posx[vv] = x;
	  posy[vv] = y;
	  //printf("interno: %d   ( %lf , %lf )\n",g.id(vv),posx[vv],posy[vv]);
	  break;
	}
      }
    } else {
      printf("GenerateVertexPositions: Error to obtain vertex coordinates.\n");
      return(false);
    }
  }
  //for (ListGraph::NodeIt vv(g); vv!=INVALID; ++vv) {
  // printf(" %d   ( %lf , %lf )\n",g.id(vv),posx[vv],posy[vv]);
  //}
  //printf("-------------------------------------------------------\n");
  return(true);
}

bool ReadListDigraph(string filename,
		     ListDigraph &g,
		     ListDigraph::NodeMap<string>& vname,
		     ListDigraph::ArcMap<double>& weight,
		     ListDigraph::NodeMap<double>& posx,
		     ListDigraph::NodeMap<double>& posy,
		     const bool dupla)
{
  ifstream ifile;
  int i,n,m;
  double peso;
  ListDigraph::Arc a;
  char nomeu[100],nomev[100];
  string STR;
  DiNode u,v;
  std::tr1::unordered_map<string,DiNode> string2node;

  ifile.open(filename.c_str());
  if (!ifile) {cout << "File '" << filename << "' does not exist.\n"; exit(0);}
  PulaBrancoComentario(ifile);
  ifile >> n;    ifile >> m; // first line have number of nodes and number of arcs
  if (m<0 || ifile.eof())
    { cout<<"File "<<filename<<" is not a digraph given by arcs.\n"; exit(0);}

  for (i=0;i<n;i++) {
    getline(ifile,STR);
    if (ifile.eof()) {cout<<"Reached unexpected end of file "<<filename<<".\n";exit(0);}
    while (STR=="") getline(ifile,STR);
    {
      string token;
      istringstream ins; // Declare an input string stream.
      ins.str(STR);        // Specify string to read.
      int nt = 0;
      while( getline(ins, token, ' ') ) {
	// format: <node_name>  <pos_x>  <pos_y>
	if (nt==0) {
	  auto test = string2node.find(token);
	  if (test!=string2node.end()){cout<<"ERROR: Repeated node: "<<nomev<<endl;exit(0);}
	  v = g.addNode(); string2node[token] = v; vname[v] = token;}
	else if (nt==1) { posx[v] = atof(token.c_str());}
	else if (nt==2) { posy[v] = atof(token.c_str());}
	nt++;
      }
    }
  }
  for (i=0;i<m;i++) {
    // format: <node_source>   <node_target>   <arc_weight>
    ifile >> nomeu;  ifile >> nomev; ifile >> peso;
    if (ifile.eof()) 
      {cout << "Reached unexpected end of file " <<filename << ".\n"; exit(0);}
    auto test = string2node.find(nomeu);
    if (test == string2node.end()) {cout<<"ERROR: Unknown node: "<<nomeu<<endl;exit(0);}
    else u = string2node[nomeu];
    
    test = string2node.find(nomev);
    if (test == string2node.end()) {cout<<"ERROR: Unknown node: "<<nomev<<endl;exit(0);}
    else v = string2node[nomev];
    a = g.addArc(u,v); weight[a] = peso;
    if (dupla) {a = g.addArc(v,u);   weight[a] = peso;}
  }
  ifile.close();
  return(true);
}


// Simplificar rotinas de leitura, fazendo rotina que le uma linha, como
// (ifile, vname, double1, double2)    // para ler node, posx, posy
// (ifile, vname, vname, double)       // para ler node_u, node_v, weight
// e outras variacoes para outros problemas.
// Colocar em myutils.cpp
//
// Para ser usado pelas rotinas de leitura de vertices e arestas.
//.....


// To read list of nodes in the format: <node_name>  <double1>  <double2>
void ReadListGraphNodes(ListGraph &g,int nnodes,ifstream & ifile,
			std::tr1::unordered_map<string,Node> & string2node,
			   ListGraph::NodeMap<string>& vname,
			   ListGraph::NodeMap<double>& double1,
			   ListGraph::NodeMap<double>& double2)
{ string STR;  Node u,v;  string token;
  for (int i=0;i<nnodes;i++) {
    getline(ifile,STR);
    if (ifile.eof()) {cout<<"Reached unexpected end of file.\n";exit(0);}
    while (STR=="") getline(ifile,STR);
    { istringstream ins; // Declare an input string stream.
      ins.str(STR);        // Specify string to read.
      for (int p=0; getline(ins, token, ' ') ; p++) {
	if (p==0) {   // For example, to read:   node_name   posx   posy
	  auto test = string2node.find(token);
	  if (test!=string2node.end()){cout<<"ERROR: Repeated node: "<<token<<endl;exit(0);}
	  v = g.addNode(); string2node[token] = v; vname[v] = token;
	  double1[v]=DBL_MAX;  double2[v]=DBL_MAX; }
	else if (p==1) { double1[v] = atof(token.c_str());}
	else if (p==2) { double2[v] = atof(token.c_str());}
      }}
  }
}


void ReadListGraphEdges(ListGraph &g,int nedges,ifstream & ifile,
			std::tr1::unordered_map<string,Node> & string2node,
			EdgeWeight &weight)
{
  Node u,v;
  Edge a;
  string nomeu,nomev;
  double peso;
  for (int i=0;i<nedges;i++) {
    // format: <node_u>   <node_v>   <edge_weight>
    ifile >> nomeu;  ifile >> nomev; ifile >> peso;
    if (ifile.eof()) 
      {cout << "Reached unexpected end of file.\n"; exit(0);}
    auto test = string2node.find(nomeu);
    if (test == string2node.end()) {cout<<"ERROR: Unknown node: "<<nomeu<<endl;exit(0);}
    else u = string2node[nomeu];
    
    test = string2node.find(nomev);
    if (test == string2node.end()) {cout<<"ERROR: Unknown node: "<<nomev<<endl;exit(0);}
    else v = string2node[nomev];
    a = g.addEdge(u,v);
    weight[a] = peso;
  }
}


// this routine is called by ReadListGraph.
bool ReadListGraph2(string filename,
		    ListGraph &g,
		    ListGraph::NodeMap<string>& vname,
		    ListGraph::EdgeMap<double>& weight,
		    ListGraph::NodeMap<double>& posx,
		    ListGraph::NodeMap<double>& posy)
{
  ifstream ifile;
  int n,m;
  ListGraph::Edge a;
  string STR;
  Node u,v;
  std::tr1::unordered_map<string,Node> string2node;

  ifile.open(filename.c_str());
  if (!ifile) {cout << "File '" << filename << "' does not exist.\n"; exit(0);}
  PulaBrancoComentario(ifile);
  ifile >> n;    ifile >> m; // first line have number of nodes and number of arcs
  if (m<0 || ifile.eof())
    { cout<<"File "<<filename<<" is not a digraph given by arcs.\n"; exit(0);}

  // continue to read ifile, and insert information in g for the next n nodes
  ReadListGraphNodes(g,n,ifile,string2node,vname,posx,posy);
  // continue to read ifile, and obtain edge weights
  ReadListGraphEdges(g,m,ifile,string2node,weight);
  // if there exists some node without pre-defined position,
  // generate all node positions
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    if (posx[v]==DBL_MAX || posy[v]==DBL_MAX) {
      GenerateVertexPositions(g,weight,posx,posy);
      break;
    }}
  ifile.close();
  return(true);
}

bool ReadListGraph(string filename,
		   ListGraph &g,
		   ListGraph::NodeMap<string>& nodename,
		   ListGraph::EdgeMap<double>& custo,
		   ListGraph::NodeMap<double>& posx,
		   ListGraph::NodeMap<double>& posy)
{
  int n,m;
  bool r;
  ifstream ifile;  
  char fname[1000];
  strcpy(fname,filename.c_str());

  ifile.open(fname);  if (!ifile) return(false);

  PulaBrancoComentario(ifile);
  ifile >> n;    ifile >> m;
  ifile.close();
  if (n<=0){cout<<"Wrong number of vertices in file "<<filename<<".\n";exit(0);}
  if (m==-1) r = ReadEuclideanListGraph(filename,g,nodename,custo,posx,posy);
  else r = ReadListGraph2(filename,g,nodename,custo,posx,posy);
  return(r);
}

/* --------------------------------------------------------------------------------------------------------------------- */
/* -----------------------------------------------------BAIXO--------------------------------------------------------- */
/* --------------------------------------------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------------------------------------------- */

void ReadListGraphNodes2(ListGraph &g,int nnodes,ifstream & ifile,
			std::tr1::unordered_map<string,Node> & string2node,
			   ListGraph::NodeMap<string>& vname,
			   ListGraph::NodeMap<double>& double1,
			   ListGraph::NodeMap<double>& double2,
			   NodeWeight &weight)
{ string STR;  Node v;  string nomev;
  double px, py, peso;
  for (int i=0;i<nnodes;i++) {
    ifile >> nomev;  
    ifile >> px;
    ifile >> py;
    ifile >> peso;    
    getline(ifile,STR);
    if (ifile.eof()) {cout<<"Reached unexpected end of file.\n";exit(0);}
    
    v = g.addNode(); string2node[nomev] = v; vname[v] = nomev;
	
    double1[v] = px;
    double2[v] = py;
    weight[v] = peso;
  }
}


void ReadListGraphEdges2(ListGraph &g,int nedges,ifstream & ifile,
			std::tr1::unordered_map<string,Node> & string2node)
{
  Node u,v;
  Edge a;
  string nomeu,nomev;
  for (int i=0;i<nedges;i++) {
    // format: <node_u>   <node_v>
    ifile >> nomeu;  ifile >> nomev;

    if (ifile.eof()) 
      {cout << "Reached unexpected end of file.\n"; exit(0);}
    auto test = string2node.find(nomeu);
    if (test == string2node.end()) {cout<<"ERROR: Unknown node X: "<<nomeu<<endl;exit(0);}
    else u = string2node[nomeu];
    
    test = string2node.find(nomev);
    if (test == string2node.end()) {cout<<"ERROR: Unknown nodeY: "<<nomev<<endl;exit(0);}
    else v = string2node[nomev];
    a = g.addEdge(u,v);
  }
}

// this routine is called by ReadListGraph.
bool ReadListGraph4(string filename,
		    ListGraph &g,
		    ListGraph::NodeMap<string>& vname,
		    ListGraph::NodeMap<double>& weight,
		    ListGraph::NodeMap<double>& posx,
		    ListGraph::NodeMap<double>& posy)
{
  ifstream ifile;
  int n,m;
  ListGraph::Edge a;
  string STR;
  Node u,v;
  std::tr1::unordered_map<string,Node> string2node;

  ifile.open(filename.c_str());
  if (!ifile) {cout << "File '" << filename << "' does not exist.\n"; exit(0);}
  PulaBrancoComentario(ifile);
  ifile >> n;    ifile >> m; // first line have number of nodes and number of arcs
  if (m<0 || ifile.eof())
    { cout<<"File "<<filename<<" is not a digraph given by arcs.\n"; exit(0);}

  // continue to read ifile, and insert information in g for the next n nodes
  ReadListGraphNodes2(g,n,ifile,string2node,vname,posx,posy, weight);
  // continue to read ifile, and obtain edge weights
  ReadListGraphEdges2(g,m,ifile,string2node);
  // if there exists some node without pre-defined position,
  // generate all node positions
  /*  
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    if (posx[v]==DBL_MAX || posy[v]==DBL_MAX) {
      GenerateVertexPositions(g,weight,posx,posy);
      break;
    }}
  */  
  ifile.close();
  return(true);
}

bool ReadListGraph3(string filename,
		   ListGraph &g,
		   ListGraph::NodeMap<string>& nodename,
		   ListGraph::NodeMap<double>& custo,
		   ListGraph::NodeMap<double>& posx,
		   ListGraph::NodeMap<double>& posy)
{
  int n,m;
  bool r;
  ifstream ifile;  
  char fname[1000];
  strcpy(fname,filename.c_str());
  ifile.open(fname);  if (!ifile) return(false);

  PulaBrancoComentario(ifile);
  ifile >> n;    ifile >> m;
  ifile.close();
  if (n<=0){cout<<"Wrong number of vertices in file "<<filename<<".\n";exit(0);}
  if (m <=0){ cout<<"Wrong number of edges in file "<<filename<<".\n";exit(0);}
  r = ReadListGraph4(filename,g,nodename,custo,posx,posy);
  return(r);
}


/* -------------------------------------------------------CIMA------------------------------------------------------------- */
/* --------------------------------------------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------------------------------------------- */

//Generate a random complete euclidean ListGraph
bool GenerateRandomEuclideanListDigraph(ListDigraph &g,
			  ListDigraph::NodeMap<string> &vname, // node name
			  ListDigraph::NodeMap<double>& px, // x-position of the node
			  ListDigraph::NodeMap<double>& py, // y-position of the node
			  ListDigraph::ArcMap<double>& weight, // weight of edges
			  int n, // number of nodes
			  double SizeX, // coordinate x is a random number in [0,SizeX)
			  double SizeY) // coordinate y is a random number in [0,SizeY)
{
  DiNode *V;
  V = new DiNode[n];
  if (V==NULL){
    cout << "Memory allocation error, number of nodes " << n << " too large\n";
    exit(0);}
  
  for (int i=0;i<n;i++) {   // insert nodes (random points in [0,100] x [0,100] )
    V[i] = g.addNode();    // new node
    px[V[i]] = SizeX*drand48();
    py[V[i]] = SizeY*drand48();
    vname[V[i]] = IntToString(i+1);// name of node is i+1
  }
  for (int i=0;i<n;i++)
    for (int j=i+1;j<n;j++) {
      Arc e = g.addArc(V[i],V[j]);  // generate new arc from v_i to v_j
      weight[e] = sqrt(pow(px[V[i]]-px[V[j]],2) + pow(py[V[i]]-py[V[j]],2));
      Arc f = g.addArc(V[j],V[i]);  // generate new arc from v_j to v_i
      weight[f] = weight[e];
    }
  delete[] V;
  return(true);
}


// Generate a triangulated ListGraph, building the Delaunay
// triangulation of random points
// Uses the geompack program, available in
// http://people.sc.fsu.edu/~jburkardt/cpp_src/geompack/geompack.html
bool GenerateTriangulatedListGraph(ListGraph &g, // return with generated graph
			  ListGraph::NodeMap<string> &vname, // return with name of the nodes
			  ListGraph::NodeMap<double>& px, // return with x-position of the nodes
			  ListGraph::NodeMap<double>& py, // return with y-position of the nodes
			  ListGraph::EdgeMap<double>& weight, // return with weight of edges
			  int n, // number of nodes 
			  double SizeX, // coordinate x is a random number in [0,SizeX)
			  double SizeY) // coordinate y is a random number in [0,SizeY)
{
  int i; // n=number of nodes
  int ntri; // number of Delaunay triangles
  Node *V = new Node[n];
  double *p = new double[2*n+2];// node coodinates are (x;y) = ( p[2*i] ; p[2*i+1] )
  int *tri = new int[6*n]; // Each 3 elements are the indexes of a triangle
  int *tri_nabe = new int[6*n];
  
  if ((V==NULL)||(p==NULL)||(tri==NULL)||(tri_nabe==NULL)){
    cout << "Memory allocation error, number of nodes " << n << " too large\n";
    exit(0);}

  for (i=0;i<n;i++) {
    V[i] = g.addNode();    // gera um vértice nó do grafo 
    px[V[i]] = SizeX*drand48();  // nodes are random points
    py[V[i]] = SizeY*drand48();
    p[2*i]   = px[V[i]];
    p[2*i+1] = py[V[i]];
    vname[V[i]] = IntToString(i+1);  // name of the node is i+1
  }
  if (r8tris2 ( n, p, &ntri, tri,  tri_nabe )) {  printf("ERROR\n");Pause(); }
  for (i=0;i<ntri;i++) { 
    int a,b,c;
    a = tri[3*i]-1; b = tri[3*i+1]-1; c = tri[3*i+2]-1;
    // each triangle if formed with nodes  V[a] , V[b] , V[c]
    // insert edges without duplications
    if ((findEdge(g,V[a],V[b])==INVALID) && (findEdge(g,V[b],V[a])==INVALID)){
      Edge e = g.addEdge(V[a],V[b]);  
      weight[e] = sqrt(pow(px[V[a]]-px[V[b]],2) + pow(py[V[a]]-py[V[b]],2));
    }
    if ((findEdge(g,V[a],V[c])==INVALID)  && (findEdge(g,V[c],V[a])==INVALID)){
      Edge e = g.addEdge(V[a],V[c]);  
      weight[e] = sqrt(pow(px[V[a]]-px[V[c]],2) + pow(py[V[a]]-py[V[c]],2));
    }
    if ((findEdge(g,V[b],V[c])==INVALID) && (findEdge(g,V[c],V[b])==INVALID)) {
      Edge e = g.addEdge(V[b],V[c]);  
      weight[e] = sqrt(pow(px[V[b]]-px[V[c]],2) + pow(py[V[b]]-py[V[c]],2));
    }
  }
  delete[] V;
  delete[] p;
  delete[] tri;
  delete[] tri_nabe;
  return(true);
}

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
		 double SizeY) // coordinate y is a random number in [0,SizeY)
{
 int i; // n=number of nodes
 ListDigraph::Node *V = new ListDigraph::Node[n];
 double *p = new double[2*n+2];// node coodinates are (x;y) = ( p[2*i] ; p[2*i+1] )
 int *tri = new int[6*n]; // Each 3 elements are the indexes of a triangle
 int ntri;
 int *tri_nabe = new int[6*n];
  
  if ((V==NULL)||(p==NULL)||(tri==NULL)||(tri_nabe==NULL)){
    cout << "Memory allocation error, number of nodes " << n << " too large\n";
    exit(0);}

 for (i=0;i<n;i++) {   // insere os vértices (pontos aleatórios no plano [0,SizeX]x[0,SizeY] )
    V[i] = g.addNode(); // gera um vértice nó do grafo 
    px[V[i]] = SizeX*drand48();  // nodes are random points
    py[V[i]] = SizeY*drand48();
    p[2*i]   = px[V[i]];
    p[2*i+1] = py[V[i]];
    vname[V[i]] = IntToString(i+1);// name of node is i+1
  }
  if (r8tris2 ( n, p, &ntri, tri,  tri_nabe )) {  printf("ERROR\n");Pause(); }
  for (i=0;i<ntri;i++) { 
    int a,b,c;
    a = tri[3*i]-1; b = tri[3*i+1]-1; c = tri[3*i+2]-1;
    // each triangle if formed with nodes  V[a] , V[b] , V[c]
    // insert arcs without duplications in both senses
    if (findArc(g,V[a],V[b])==INVALID){
      ListDigraph::Arc e = g.addArc(V[a],V[b]); 
      weight[e] = sqrt(pow(px[V[a]]-px[V[b]],2) + pow(py[V[a]]-py[V[b]],2));
    }
    if (findArc(g,V[b],V[a])==INVALID){
      ListDigraph::Arc e = g.addArc(V[b],V[a]);  
      weight[e] = sqrt(pow(px[V[b]]-px[V[a]],2) + pow(py[V[b]]-py[V[a]],2));
    }
    if (findArc(g,V[a],V[c])==INVALID){
      ListDigraph::Arc e = g.addArc(V[a],V[c]);  
      weight[e] = sqrt(pow(px[V[a]]-px[V[c]],2) + pow(py[V[a]]-py[V[c]],2));
    }
    if (findArc(g,V[c],V[a])==INVALID){
      ListDigraph::Arc e = g.addArc(V[c],V[a]);  
      weight[e] = sqrt(pow(px[V[c]]-px[V[a]],2) + pow(py[V[c]]-py[V[a]],2));
    }
    if (findArc(g,V[b],V[c])==INVALID) {
      ListDigraph::Arc e = g.addArc(V[b],V[c]);  
      weight[e] = sqrt(pow(px[V[b]]-px[V[c]],2) + pow(py[V[b]]-py[V[c]],2));
    }
    if (findArc(g,V[c],V[b])==INVALID) {
      ListDigraph::Arc e = g.addArc(V[c],V[b]);  
      weight[e] = sqrt(pow(px[V[c]]-px[V[b]],2) + pow(py[V[c]]-py[V[b]],2));
    }
  }
  delete[] V;
  delete[] p;
  delete[] tri;
  delete[] tri_nabe;
  return(true);
}


//Generate a random complete euclidean ListGraph
bool GenerateRandomEuclideanListGraph(ListGraph &g,
		  ListGraph::NodeMap<string> &vname, // node names
		  ListGraph::NodeMap<double>& px, // x-position of the nodes
		  ListGraph::NodeMap<double>& py, // y-position of the nodes
		  ListGraph::EdgeMap<double>& weight, // weight of edges
		  int n, // number of nodes
		  double SizeX, // coordinate x is a random number in [0,SizeX)
		  double SizeY) // coordinate y is a random number in [0,SizeY)
{
  int i,j; // n=number of nodes
  Node *V;
  V = new Node[n];
  if (V==NULL){
    cout << "Memory allocation error, number of nodes " << n << " too large\n";
    exit(0);}
  
  for (i=0;i<n;i++) {   // insert nodes (random points in [0,100] x [0,100] )
    V[i] = g.addNode(); // generate a new node
    px[V[i]] = SizeX*drand48();
    py[V[i]] = SizeY*drand48();
    vname[V[i]] = IntToString(i+1); // name of the node is i+1
  }
  for (i=0;i<n;i++)
    for (j=i+1;j<n;j++) {
      Edge e = g.addEdge(V[i],V[j]);  // generate an edge
      weight[e] = sqrt(pow(px[V[i]]-px[V[j]],2) + pow(py[V[i]]-py[V[j]],2));
    }
  delete[] V;
  return(true);
}


// // View an undirected graph with DOT or NEATO program (both available in the 
// // graphviz package. If you do not have the positions of the vertices of a graph,
// // you can use the routine GenerateVertexPositions, available in the file 
// // readgraph.cpp.
// int ViewListGraph(ListGraph &g,
// 	      int DOT_or_NEATO, // can be VIEW_DOT or VIEW_NEATO
// 	      ListGraph::NodeMap<string> &vname, // name of the vertices
// 	      ListGraph::EdgeMap<string> &ename,  // name of edges
// 	      ListGraph::NodeMap<double>& posx, // position x of the vertices
// 	      ListGraph::NodeMap<double>& posy, // position y of the vertices
// 	      ListGraph::NodeMap<int> &vcolor,   // color of the vertices
// 	      ListGraph::EdgeMap<int> &acolor)  // color of edges
// {
//   char tempname[1000],cmd[1000];
//   FILE *fp;
//   double minposx=10000000,minposy=10000000,maxposx=-10000000,maxposy=-10000000,
//     delta,factor;

//   (void) vname;// to avoid "non-used" parameter.
//   strcpy(tempname,".viewgraphtempname");
//   fp = fopen(tempname,"w+");
//   // fp = stdout;
//   if (fp==NULL) {
//     cout << "Erro ao abrir arquivo para visualizar o grafo.\n";
//     return(0);
//   }
//   for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
//     if (posx[v] < minposx) minposx = posx[v];
//     if (posx[v] > maxposx) maxposx = posx[v];
//     if (posy[v] < minposy) minposy = posy[v];
//     if (posy[v] > maxposy) maxposy = posy[v];
//   }
//   factor = 30;
//   delta = fmax(maxposx - minposx,maxposy - minposy);

//   fprintf(fp,"graph g {\n");
//   // fprintf(fp,"\t\"start\" ;\n");
//   fprintf(fp,"\tsize = \"10, 10\";\n");
//   fprintf(fp,"\tnode [shape = \"circle\"];\n");
//   for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
//     fprintf(fp,"\t%s [color=\"%s\", pos = \"%lf,%lf!\" ];\n",vname[v].c_str(),ColorName(vcolor[v]),factor*(posx[v]-minposx)/delta,factor*(posy[v]-minposy)/delta);
//   }
//   for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
//     if (acolor[a]!=WHITE)
//       fprintf(fp,"\t%s  -- %s [label = \"%s\", color=\"%s\" ];\n",vname[g.u(a)].c_str(),vname[g.v(a)].c_str(),ename[a].c_str(),ColorName(acolor[a]));
//   }
//   fprintf(fp,"}\n");
//   fclose(fp);

//   if (DOT_or_NEATO==VIEW_NEATO) sprintf(cmd,"neato -Tpdf %s -o %s.pdf",tempname,tempname); 
//   else sprintf(cmd,"dot -Tpdf %s -o %s.pdf",tempname,tempname);
//   system(cmd);

//   sprintf(cmd,"%s %s.pdf",pdfreader.c_str(),tempname); system(cmd);
//   return(1);
// }



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
	      string text) // text displayed below the figure
{
  char tempname[1000],cmd[1000],outputname[1000];
  FILE *fp;
  double minpx=DBL_MAX,minpy=DBL_MAX,maxpx=-DBL_MAX,maxpy=-DBL_MAX,delta,factor;

  // obtain a temporary file name
  strcpy(tempname,".viewgraphtempname");
  sprintf(outputname,"%s.pdf",tempname);
  fp = fopen(tempname,"w+");
  if (fp==NULL) {cout << "Error to open temporary file to visualize graph.\n"; return(0);}
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    if (px[v] < minpx) minpx = px[v];
    if (px[v] > maxpx) maxpx = px[v];
    if (py[v] < minpy) minpy = py[v];
    if (py[v] > maxpy) maxpy = py[v];
  }
  factor = 40; // quanto maior, menor o desenho do vértice
  delta = fmax(maxpx - minpx,maxpy - minpy);
  // Generate a text file with the graph format of neato program
  fprintf(fp,"graph g {\n");
  cout << "\tsize = \"" << graphlabelfontsize << ", " << graphlabelfontsize << "\";\n";
  //fprintf(fp,"\tsize = \"10, 10\";\n");
  fprintf(fp,"\tnode [shape = \"circle\", fontsize = \"%d\"];\n",graphnodefontsize);
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    if (vcolor[v]==NOCOLOR) continue;
    fprintf(fp,"\t%s [style = filled, fontsize = \"%d\", color=\"%s\", pos = \"%lf,%lf!\" ];\n",
	    vname[v].c_str(),graphnodefontsize,ColorName(vcolor[v]).c_str(),factor*(px[v]-minpx)/delta,factor*(py[v]-minpy)/delta);
  }
  for (ListGraph::EdgeIt e(g); e!=INVALID; ++e) {
    if (ecolor[e]==NOCOLOR) continue;
    fprintf(fp,"\t%s  -- %s [label = \"%s\", color=\"%s\", fontsize=\"%d\" ];\n",vname[g.u(e)].c_str(),vname[g.v(e)].c_str(),ename[e].c_str(),ColorName(ecolor[e]).c_str(), graphedgefontsize);
  }
  fprintf(fp,"label=\"%s\";\nfontsize=\"%d\";\n",text.c_str(),graphlabelfontsize);
  fprintf(fp,"}\n");
  fclose(fp);
  sprintf(cmd,"neato -Tpdf %s -o %s",tempname,outputname); system(cmd);
  view_pdffile(outputname);
  return(1);
}



// // View an undirected graph with DOT or NEATO program (both available in the 
// // graphviz package. If you do not have the positions of the vertices of a graph,
// // you can use the routine GenerateVertexPositions, available in the file 
// // readgraph.cpp.
// int ViewListGraph(ListGraph &g,
// 	      int DOT_or_NEATO, // can be VIEW_DOT or VIEW_NEATO
// 	      ListGraph::NodeMap<string> &vname, // name of the vertices
// 	      ListGraph::EdgeMap<string> &ename,  // name of edges
// 	      ListGraph::NodeMap<double>& posx, // position x of the vertices
// 	      ListGraph::NodeMap<double>& posy, // position y of the vertices
// 	      ListGraph::NodeMap<int> &vcolor,   // color of the vertices
// 	      ListGraph::EdgeMap<int> &acolor)  // color of edges
// {
//   char tempname[1000],cmd[1000];
//   FILE *fp;
//   double minposx=10000000,minposy=10000000,maxposx=-10000000,maxposy=-10000000,
//     delta,factor;

//   (void) vname;// to avoid "non-used" parameter.
//   strcpy(tempname,".viewgraphtempname");
//   fp = fopen(tempname,"w+");
//   // fp = stdout;
//   if (fp==NULL) {
//     cout << "Erro ao abrir arquivo para visualizar o grafo.\n";
//     return(0);
//   }
//   for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
//     if (posx[v] < minposx) minposx = posx[v];
//     if (posx[v] > maxposx) maxposx = posx[v];
//     if (posy[v] < minposy) minposy = posy[v];
//     if (posy[v] > maxposy) maxposy = posy[v];
//   }
//   factor = 30;
//   delta = fmax(maxposx - minposx,maxposy - minposy);

//   fprintf(fp,"graph g {\n");
//   // fprintf(fp,"\t\"start\" ;\n");
//   fprintf(fp,"\tsize = \"10, 10\";\n");
//   fprintf(fp,"\tnode [shape = \"circle\"];\n");
//   for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
//     fprintf(fp,"\t%s [color=\"%s\", pos = \"%lf,%lf!\" ];\n",vname[v].c_str(),ColorName(vcolor[v]).c_str(),factor*(posx[v]-minposx)/delta,factor*(posy[v]-minposy)/delta);
//   }
//   for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
//     if (acolor[a]!=WHITE)
//       fprintf(fp,"\t%s  -- %s [label = \"%s\", color=\"%s\" ];\n",vname[g.u(a)].c_str(),vname[g.v(a)].c_str(),ename[a].c_str(),ColorName(acolor[a]).c_str());
//   }
//   fprintf(fp,"}\n");
//   fclose(fp);

//   if (DOT_or_NEATO==VIEW_NEATO) sprintf(cmd,"neato -Tpdf %s -o %s.pdf",tempname,tempname); 
//   else sprintf(cmd,"dot -Tpdf %s -o %s.pdf",tempname,tempname);
//   system(cmd);

//   sprintf(cmd,"%s %s.pdf",pdfreader.c_str(),tempname); system(cmd);
//   return(1);
// }



// int ViewListDigraph(ListDigraph &g,
// 	      int DOT_or_NEATO,
// 	      ListDigraph::NodeMap<string> &vname, // nome dos vertices
// 	      ListDigraph::ArcMap<string> &ename,  // nome das arestas (e.g. peso dela)
// 	      ListDigraph::NodeMap<int> &vcolor,   // cor dos vertices
// 	      ListDigraph::ArcMap<int> &acolor)    // cor das arestas
// {
//   char tempname[1000],cmd[1000];
//   FILE *fp;

//   (void) vname;// to avoid "non-used" parameter.
//   strcpy(tempname,".viewdigraphtempname");
//   fp = fopen(tempname,"w+");
//   // fp = stdout;
//   if (fp==NULL) {
//     cout << "Erro ao abrir arquivo para visualizar o grafo.\n";
//     return(0);
//   }
//   fprintf(fp,"digraph g {\n");
//   // fprintf(fp,"\t\"start\" ;\n");
//   // fprintf(fp,"\tgraph;\n");
//   fprintf(fp,"\tedge [fontsize=9];\n");
//   fprintf(fp,"\tnode [style=filled, fontsize=10];\n");
//   fprintf(fp,"\tnodesep = .3;\n");
//   fprintf(fp,"\tranksep = .5;\n");
//   //fprintf(fp,"\tsize = \"expand\";\n");
//   fprintf(fp,"\tnode [shape = circle];\n");
//   for (ListDigraph::NodeIt v(g); v!=INVALID; ++v) {
//     fprintf(fp,"\t%s [color=%s, width=.3, height=.3];\n",vname[v].c_str(),ColorName(vcolor[v]).c_str());
//   }
//   for (ListDigraph::ArcIt a(g); a!=INVALID; ++a) {
//     if (acolor[a]!=WHITE)
//       fprintf(fp,"\t%s  -> %s [label = \"%s\", color=%s ];\n",vname[g.source(a)].c_str(),vname[g.target(a)].c_str(),ename[a].c_str(),ColorName(acolor[a]).c_str());
//   }
//   fprintf(fp,"}\n");
//   fclose(fp);

//   sprintf(cmd,"neato -Tpdf %s -o %s.pdf",tempname,tempname); 
//   system(cmd);
//   sprintf(cmd,"%s %s.pdf",pdfreader.c_str(),tempname); system(cmd);
//   return(1);
// }


// This routine visualize a digraph using a pdf viewer. It uses neato (from
// graphviz.org) to generate a pdf file and a program to view the pdf file. The
// pdf viewer name is given in the viewername parameter.
int ViewListDigraph(ListDigraph &g,
      ListDigraph::NodeMap<string> &vname, // node names
      ListDigraph::NodeMap<double>& px, // x-position of the nodes
      ListDigraph::NodeMap<double>& py, // y-position of the nodes
      ListDigraph::NodeMap<int>&    vcolor, // color of node (see myutils.h)
      ListDigraph::ArcMap<int>&    ecolor, // color of edge 
      string text) // text displayed below the figure
{
  char tempname[1000],cmd[1000];
  FILE *fp;
  double minpx=DBL_MAX,minpy=DBL_MAX,maxpx=-DBL_MAX,maxpy=-DBL_MAX,delta,factor;
  string str;

  // obtain a temporary file name
  strcpy(tempname,".viewdigraphtempname");
  fp = fopen(tempname,"w+");
  if (fp==NULL) {cout << "Error to open temporary file to visualize digraph.\n"; return(0);}
  for (ListDigraph::NodeIt v(g); v!=INVALID; ++v) {
    if (px[v] < minpx) minpx = px[v];
    if (px[v] > maxpx) maxpx = px[v];
    if (py[v] < minpy) minpy = py[v];
    if (py[v] > maxpy) maxpy = py[v];
  }
  factor = 40; // using larger values makes small nodes
  delta = fmax(maxpx - minpx,maxpy - minpy);
  // Generate a text file with the graph format of neato program
  fprintf(fp,"digraph g {\n");
  cout << "\tsize = \"" << graphlabelfontsize << ", " << graphlabelfontsize << "\";\n";
  //fprintf(fp,"\tsize = \"10, 10\";\n");
  fprintf(fp,"\tnode [style = filled, shape = \"circle\", fontsize = \"%d\"];\n",graphnodefontsize);
  for (ListDigraph::NodeIt v(g); v!=INVALID; ++v) {
    if (vcolor[v]==NOCOLOR) continue;
    fprintf(fp,"\t%s [color=\"%s\", pos = \"%lf,%lf!\"];\n",
	    vname[v].c_str(),ColorName(vcolor[v]).c_str(),factor*(px[v]-minpx)/delta,factor*(py[v]-minpy)/delta);
  }
  for (ListDigraph::ArcIt e(g); e!=INVALID; ++e) {
    if (ecolor[e]==NOCOLOR) continue;
    fprintf(fp,"\t%s -> %s [color=\"%s\" ];\n",vname[g.source(e)].c_str(),vname[g.target(e)].c_str(),ColorName(ecolor[e]).c_str());
  }
  fprintf(fp,"label=\"%s\";\nfontsize=%d;\n",text.c_str(),graphlabelfontsize);
  fprintf(fp,"}\n");
  fclose(fp);
  sprintf(cmd,"neato -Tpdf %s -o %s.pdf",tempname,tempname); system(cmd);
  str = tempname;
  str = str + ".pdf";
  view_pdffile(str);
  return(1);
}



// obtain a mincut separating vertices 's' and 't' for undirected graphs
// The input is given by the graph 'g'. Each edge of g has a "weight".
// The returned cut is given by the vector of nodes 'cut' (boolean
// vector, nodes in one side have value false and nodes in the other side
// have value true).
double MinCut(ListGraph &g, EdgeWeight &weight, Node &s,Node &t, CutMap &cut)
{
  Preflow<ListGraph, EdgeWeight> pf(g, weight, s, t); 
  pf.runMinCut();
  pf.minCutMap(cut);
  return (pf.flowValue());
}

// Obtain a mininum cut for directed graphs from s to t.
// The returned cut is given by the vector of nodes 'cut' (boolean
// vector: nodes v in the same side of s have cut[v]=true, otherwise cut[v]=false.
double DiMinCut(ListDigraph &g, ArcWeight &weight, DiNode &s,DiNode &t, DiCutMap &vcut)
{
  Preflow<ListDigraph, ArcWeight> preflow_test(g, weight, s, t); 
  preflow_test.run();
  preflow_test.minCutMap(vcut);
  return (preflow_test.flowValue());
}



/*---------------------------------------------------------------------*/
#define MAXPOINTPOSITION 6000 /* each coordinate belongs to the
				  interval [0..MAXPOINTPOSITION)  */
void getepscolor(char *epscolorname,int cor)
{
  switch(cor) {
  case BLACK: strcpy(epscolorname,"col0"); return;
  case BLUE: strcpy(epscolorname,"col1"); return;
  case GREEN: strcpy(epscolorname,"col2"); return;
  case RED: strcpy(epscolorname,"col4"); return;
  case WHITE: strcpy(epscolorname,"col7"); return;
  }
  strcpy(epscolorname,"col0"); /* sem cor definida fica com preto */
}

void PrintListGraph(ListGraph &g, ListGraph::NodeMap<string> &vname, EdgeWeight &graphweight)
{
  int Nnodes = countNodes(g); // number of nodes in the input graph
  int Nedges = countEdges(g); // number of edges in the input graph

  printf("-------------------------------------------------------\n");
  printf("Number of nodes: %d\n",Nnodes);
  printf("Number of edges: %d\n",Nedges);
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) printf("%s\n",vname[v].c_str()); printf("\n");
  printf("-------------------------------------------------------\n");
  for (ListGraph::EdgeIt a(g); a!=INVALID; ++a)
    printf("%s -- %s  %lf\n",vname[g.u(a)].c_str(),vname[g.v(a)].c_str(), graphweight[a]);
  printf("\n");
}


int ViewEuclideanListGraph(ListGraph &g,
		       ListGraph::NodeMap<string> &vname, // nome dos vertices
		       ListGraph::NodeMap<double> &posx, // coord. x dos vertices
		       ListGraph::NodeMap<double> &posy, // coord. y dos vertices
		       ListGraph::NodeMap<int> &vcolor,  // cor dos vertices
		       ListGraph::EdgeMap<int> &ecolor)  // cor das arestas
{
  char tempname[1000],cmd[1000];
  FILE *fp;
  double gap,maxx, maxy, minx, miny,
    telax,posxu,posxv,posyu,posyv;
  char epscolor[100];
  (void) vname;// to avoid "non-used" parameter.
  string str;
  
  strcpy(tempname,".vieweuclideangraphtempname");
  fp = fopen(tempname,"w+");
  if (fp==NULL) {
    cout << "Erro ao abrir arquivo para visualizar o grafo.\n";
    return(0);
  }

  maxx = -999999;
  maxy = -999999;
  minx =  999999;
  miny =  999999;
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    if (posx[v] > maxx) maxx =  posx[v];
    if (posx[v] < minx) minx =  posx[v];
    if (posy[v] > maxy) maxy =  posy[v];
    if (posy[v] < miny) miny =  posy[v];
  }
  telax = 500;
  fprintf(fp,"%%!PS-Adobe-2.0 EPSF-2.0\n");
  fprintf(fp,"%%%%Title: x.eps\n");
  fprintf(fp,"%%%%Creator: fig2dev Version 3.2 Patchlevel 3c\n");
  fprintf(fp,"%%%%CreationDate: Thu Sep 12 13:02:34 2002\n");
  fprintf(fp,"%%%%For: fkm@hobbes.dcc.unicamp.br ()\n");
  fprintf(fp,"%%%%BoundingBox: 0 0 %d %d\n",(int) telax,(int) telax);
  fprintf(fp,"%%%%Magnification: 1.0000\n");
  fprintf(fp,"%%%%EndComments\n");
  fprintf(fp,"/$F2psDict 200 dict def\n");
  fprintf(fp,"$F2psDict begin\n");
  fprintf(fp,"$F2psDict /mtrx matrix put\n");
  fprintf(fp,"/col-1 {0 setgray} bind def\n");
  fprintf(fp,"/col0 {0.000 0.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col1 {0.000 0.000 1.000 srgb} bind def\n");
  fprintf(fp,"/col2 {0.000 1.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col3 {0.000 1.000 1.000 srgb} bind def\n");
  fprintf(fp,"/col4 {1.000 0.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col5 {1.000 0.000 1.000 srgb} bind def\n");
  fprintf(fp,"/col6 {1.000 1.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col7 {1.000 1.000 1.000 srgb} bind def\n");
  fprintf(fp,"/col8 {0.000 0.000 0.560 srgb} bind def\n");
  fprintf(fp,"/col9 {0.000 0.000 0.690 srgb} bind def\n");
  fprintf(fp,"/col10 {0.000 0.000 0.820 srgb} bind def\n");
  fprintf(fp,"/col11 {0.530 0.810 1.000 srgb} bind def\n");
  fprintf(fp,"/col12 {0.000 0.560 0.000 srgb} bind def\n");
  fprintf(fp,"/col13 {0.000 0.690 0.000 srgb} bind def\n");
  fprintf(fp,"/col14 {0.000 0.820 0.000 srgb} bind def\n");
  fprintf(fp,"/col15 {0.000 0.560 0.560 srgb} bind def\n");
  fprintf(fp,"/col16 {0.000 0.690 0.690 srgb} bind def\n");
  fprintf(fp,"/col17 {0.000 0.820 0.820 srgb} bind def\n");
  fprintf(fp,"/col18 {0.560 0.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col19 {0.690 0.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col20 {0.820 0.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col21 {0.560 0.000 0.560 srgb} bind def\n");
  fprintf(fp,"/col22 {0.690 0.000 0.690 srgb} bind def\n");
  fprintf(fp,"/col23 {0.820 0.000 0.820 srgb} bind def\n");
  fprintf(fp,"/col24 {0.500 0.190 0.000 srgb} bind def\n");
  fprintf(fp,"/col25 {0.630 0.250 0.000 srgb} bind def\n");
  fprintf(fp,"/col26 {0.750 0.380 0.000 srgb} bind def\n");
  fprintf(fp,"/col27 {1.000 0.500 0.500 srgb} bind def\n");
  fprintf(fp,"/col28 {1.000 0.630 0.630 srgb} bind def\n");
  fprintf(fp,"/col29 {1.000 0.750 0.750 srgb} bind def\n");
  fprintf(fp,"/col30 {1.000 0.880 0.880 srgb} bind def\n");
  fprintf(fp,"/col31 {1.000 0.840 0.000 srgb} bind def\n");
  fprintf(fp,"\n");
  fprintf(fp,"end\n");
  fprintf(fp,"save\n");
  fprintf(fp,"newpath 0 %d moveto 0 0 lineto %d 0 lineto %d %d "
	  "lineto closepath clip newpath\n",(int) telax,(int) telax,(int) telax,(int) telax);
  fprintf(fp,"%d %d translate\n",-10,(int) telax+10);
  fprintf(fp,"1 -1 scale\n");
  fprintf(fp,"\n");
  fprintf(fp,"/cp {closepath} bind def\n");
  fprintf(fp,"/ef {eofill} bind def\n");
  fprintf(fp,"/gr {grestore} bind def\n");
  fprintf(fp,"/gs {gsave} bind def\n");
  fprintf(fp,"/sa {save} bind def\n");
  fprintf(fp,"/rs {restore} bind def\n");
  fprintf(fp,"/l {lineto} bind def\n");
  fprintf(fp,"/m {moveto} bind def\n");
  fprintf(fp,"/rm {rmoveto} bind def\n");
  fprintf(fp,"/n {newpath} bind def\n");
  fprintf(fp,"/s {stroke} bind def\n");
  fprintf(fp,"/sh {show} bind def\n");
  fprintf(fp,"/slc {setlinecap} bind def\n");
  fprintf(fp,"/slj {setlinejoin} bind def\n");
  fprintf(fp,"/slw {setlinewidth} bind def\n");
  fprintf(fp,"/srgb {setrgbcolor} bind def\n");
  fprintf(fp,"/rot {rotate} bind def\n");
  fprintf(fp,"/sc {scale} bind def\n");
  fprintf(fp,"/sd {setdash} bind def\n");
  fprintf(fp,"/ff {findfont} bind def\n");
  fprintf(fp,"/sf {setfont} bind def\n");
  fprintf(fp,"/scf {scalefont} bind def\n");
  fprintf(fp,"/sw {stringwidth} bind def\n");
  fprintf(fp,"/tr {translate} bind def\n");
  fprintf(fp,"/tnt {dup dup currentrgbcolor\n");
  fprintf(fp,"  4 -2 roll dup 1 exch sub 3 -1 roll mul add\n");
  fprintf(fp,"  4 -2 roll dup 1 exch sub 3 -1 roll mul add\n");
  fprintf(fp,"  4 -2 roll dup 1 exch sub 3 -1 roll mul add srgb}\n");
  fprintf(fp,"  bind def\n");
  fprintf(fp,"/shd {dup dup currentrgbcolor 4 -2 roll mul 4 -2 roll mul\n");
  fprintf(fp,"  4 -2 roll mul srgb} bind def\n");
  fprintf(fp," /DrawEllipse {\n");
  fprintf(fp,"	/endangle exch def\n");
  fprintf(fp,"	/startangle exch def\n");
  fprintf(fp,"	/yrad exch def\n");
  fprintf(fp,"	/xrad exch def\n");
  fprintf(fp,"	/y exch def\n");
  fprintf(fp,"	/x exch def\n");
  fprintf(fp,"	/savematrix mtrx currentmatrix def\n");
  fprintf(fp,"	x y tr xrad yrad sc 0 0 1 startangle endangle arc\n");
  fprintf(fp,"	closepath\n");
  fprintf(fp,"	savematrix setmatrix\n");
  fprintf(fp,"	} def\n");
  fprintf(fp,"\n");
  fprintf(fp,"/$F2psBegin {$F2psDict begin "
	  "/$F2psEnteredState save def} def\n");
  fprintf(fp,"/$F2psEnd {$F2psEnteredState restore end} def\n");
  fprintf(fp,"\n");
  fprintf(fp,"$F2psBegin\n");
  fprintf(fp,"%%%%Page: 1 1\n");
  fprintf(fp,"10 setmiterlimit\n");
  fprintf(fp," %10.8lf %10.8lf sc\n",(double) MAXPOINTPOSITION/100000.0,(double) MAXPOINTPOSITION/100000.0);
  fprintf(fp,"%%\n");
  fprintf(fp,"%% Fig objects follow\n");
  fprintf(fp,"%%\n");
  fprintf(fp,"25.000 slw\n");
  gap = 300;
  for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
    ListGraph::Node u,v;
    if (ecolor[a]==WHITE) continue;

    u = g.u(a);   v = g.v(a);
    posxu = (int) (MAXPOINTPOSITION*((double)(posx[u] -minx))/
		   ((double) (maxx-minx)))+gap;
    posyu = MAXPOINTPOSITION - (int) (MAXPOINTPOSITION*((double)(posy[u] -miny))/
				      ((double) (maxy-miny)))+gap;
    posxv = (int) (MAXPOINTPOSITION*((double)(posx[v] -minx))/
		   ((double) (maxx-minx)))+gap;
    posyv = MAXPOINTPOSITION - (int) (MAXPOINTPOSITION*((double)(posy[v] -miny))/
				      ((double) (maxy-miny)))+gap;
    getepscolor(epscolor,ecolor[a]);
    fprintf(fp,"n %d %d m\n %d %d l gs %s s gr \n",(int) posxu,(int) posyu,(int) posxv,(int) posyv,epscolor);
  }
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    posxv = (int) (MAXPOINTPOSITION*((double)(posx[v] -minx))/
		  ((double) (maxx-minx)))+gap;
    posyv = MAXPOINTPOSITION - (int) (MAXPOINTPOSITION*((double)(posy[v] -miny))/
				     ((double) (maxy-miny)))+gap;
    getepscolor(epscolor,vcolor[v]);
    /* 45 eh o raio (ambos)*/
    fprintf(fp,"n %d %d 45 45 0 360 DrawEllipse gs %s 1.00 shd ef gr \n", 
	    (int) posxv,(int) posyv,epscolor);
  }
  fprintf(fp,"$F2psEnd\n");
  fprintf(fp,"rs\n");
  fclose(fp);

  sprintf(cmd,"mv %s %s.eps",tempname,tempname);           system(cmd);
  sprintf(cmd,"convert %s.eps %s.pdf",tempname,tempname);  system(cmd);
  //sprintf(cmd,"%s %s.pdf",pdfreader.c_str(),tempname);  system(cmd);
  str = tempname;
  str = str+".pdf";
  view_pdffile(str);
  return(true);
}

// ================================================================
int globalcounter=0;

void ADJMAT_FreeNotNull(void *p){  if (p) free(p);  }

// Define an adjacency matrix, so as we have a fast access for the edges of a graph,
// given a pair of vertices. This is mainly used in the subroutine 2opt, for 
// dense graphs.
// The adjacency matrix is stored in a strictly triangular inferior matrix.

AdjacencyMatrix::AdjacencyMatrix(ListGraph &graph,EdgeWeight &graphweight,double NonEdgValue):
  Node2Index(graph),Edge2Index(graph)
{
  int i;
  g = &graph;
  NonEdgeValue = NonEdgValue;
  weight = &graphweight;
  Nnodes = countNodes(graph); // number of nodes in the input graph
  Nedges = countEdges(graph); // number of edges in the input graph
  Nmatrix = (Nnodes*(Nnodes-1))/2; // no. of edges/elem. in strict. triang. inf. matrix

  AdjMatrix = (double *) malloc(sizeof(double)*Nmatrix);
  Index2Node = (Node *) malloc(sizeof(Node)*Nnodes);
  Index2Edge = (Edge *) malloc(sizeof(Edge)*Nedges);

  if ((AdjMatrix==NULL)||(Index2Node==NULL)||(Index2Edge==NULL)) { 
    cout << "Out of memory in constructor of AdjacencyMatrix\n"; 
    ADJMAT_FreeNotNull(AdjMatrix); ADJMAT_FreeNotNull(Index2Node); ADJMAT_FreeNotNull(Index2Edge);
    exit(0);}

  i = 0;
  for (NodeIt v(*g); v != INVALID; ++v) {
    Index2Node[i] = v;
    AdjacencyMatrix::Node2Index[v]=i;
    i++;
  }

  // Initially all edges have infinity weight
  for (int i=0;i<Nmatrix;i++) AdjMatrix[i] = NonEdgeValue;
  // Then, update the existing edges with the correct weight
  for (EdgeIt e(graph); e != INVALID; ++e) {
    Node u,v;    int i_u,i_v;
    u = graph.u(e);  v = graph.v(e);  // obtain the extremities of e
    i_u = Node2Index[u];
    i_v = Node2Index[v];
    if (i_u > i_v) AdjMatrix[i_u*(i_u-1)/2+i_v] = graphweight[e];
    else if (i_u < i_v) AdjMatrix[i_v*(i_v-1)/2+i_u] = graphweight[e];
    else {
      cout << "Out of memory in constructor of AdjacencyMatrix\n"; 
      exit(0);}
  }
}

double AdjacencyMatrix::Cost(Node u,Node v)
{
  int i_u,i_v;
  i_u = Node2Index[u];
  i_v = Node2Index[v];
  globalcounter++;
  //cout << "globalcounter1 = " << globalcounter << endl;
  try{
    if (i_u > i_v) return(AdjMatrix[i_u*(i_u-1)/2+i_v]);
    else if (i_u < i_v) return(AdjMatrix[i_v*(i_v-1)/2+i_u]);
    else return(NonEdgeValue);
  }catch (...) {
    cout << "Exception: " << globalcounter << endl;
    exit(1);
  }
}

double AdjacencyMatrix::Cost(Edge e)
{
  int i_u,i_v;
  Node u,v;
  globalcounter++;
  //cout << "globalcounter2 = " << globalcounter << endl;
  u = (*g).u(e);  v = (*g).v(e);
  i_u = Node2Index[u];
  i_v = Node2Index[v];
  if (i_u > i_v) return(AdjMatrix[i_u*(i_u-1)/2+i_v]);
  else return(AdjMatrix[i_v*(i_v-1)/2+i_u]);
}


AdjacencyMatrix::~AdjacencyMatrix()
{
  ADJMAT_FreeNotNull(AdjMatrix); ADJMAT_FreeNotNull(Index2Node); ADJMAT_FreeNotNull(Index2Edge);
}




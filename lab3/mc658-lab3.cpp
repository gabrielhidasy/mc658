// Projeto e Análise de Algoritmos
// Flávio Keidi Miyazawa
// Laboratório 1
// Enunciado em www.ic.unicamp.br/~fkm/disciplinas/mc658/2013s1/index.html
#include<lemon/list_graph.h>
#include <lemon/lp.h>
#include <iostream>
//#include <stdlib.h>
#include <string>
#include <cstdio>
using namespace std;
using namespace lemon;
typedef ListGraph Graph;
typedef Graph::Node Node;
typedef Graph::NodeMap<int> Capacidade;
typedef Graph::NodeIt NodeIt;
typedef Graph::EdgeIt EdgeIt;
typedef Graph::IncEdgeIt IncEdgeIt;
typedef Graph::Edge Edge;

#define N 1000


int seleciona_propagandas(int n, double C, double V[N], double P[N],
			  double w[N][N], int S[N], double *Upper_Bound,
			  long maxtime) {
  S[0] = 1;
  *Upper_Bound = (double) S[0];
}

int main()
{
  int semente;
  cout << "Entre com a semente para geracao de numeros aleatorios: ";
  scanf("%d",&semente);
  srand48(semente);
  Graph g;                        // declare a graph g
  Graph::EdgeMap<double> Lucro(g);// declare weights for each edge of the graph
  Graph::NodeMap<int> Capacidade(g);// declare weights for each edge of the graph
  Graph::Node U[1000],P[1000];
  Graph::NodeMap<string> NomeVertice(g);
  Graph::EdgeMap<string> NomeAresta(g);
  int i,j,n,m;
  n = (int) 10*drand48()+1;
  m = (int) 30*drand48()+1;
  // insere usuarios e capacidade de cada usuario
  for (i=0;i<n;i++) {
    char s[100];
    U[i] = g.addNode();
    sprintf(s,"U[%2d]",i+1);
    NomeVertice[U[i]] = s;
    Capacidade[U[i]] = (int) 3*drand48()+1;
  }
  // insere propagandas e capacidade igual a 1
  for (j=0;j<m;j++) {
    char s[100];
    P[j] = g.addNode();
    sprintf(s,"P[%2d]",j+1);
    NomeVertice[P[j]] = s;
    Capacidade[P[j]] = 1;
  }
  // insere as arestas (o grafo é bipartido completo) e os pesos
  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      Edge e;
      e = g.addEdge(U[i],P[j]);
      NomeAresta[e] = NomeVertice[U[i]]+"---"+NomeVertice[P[j]];
      Lucro[e]=(int) 1+10*drand48();
    }
  }

  // Altere o trecho abaixo para imprimir apenas as arestas da solucao
  double LucroTotal=0.0;
  for (EdgeIt a(g); a != INVALID; ++a) {
      cout << NomeVertice[g.u(a)] << "---" << NomeVertice[g.v(a)] << "\n";
      LucroTotal += Lucro[a];
  }
  cout << "Lucro total = " << LucroTotal << "\n";
  return 0;
}

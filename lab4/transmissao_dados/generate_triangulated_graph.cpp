// MC658 - Project and Analysis of Algorithms III
// Proper input generator for lab4 (2015s2 - fkm)
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <lemon/dijkstra.h>
#include <lemon/list_graph.h>
#include <lemon/path.h>
#include "mygraphlib.h"
#include "myutils.h"
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace lemon;
using namespace std;

#define BC_EPS 0.001
#define BC_INF 1000000000000.0
#define EPS 0.00001

#define X 10.0

int main(int argc, char *argv[]) {
    int n;
    double box_width, box_height;
    int num_test = 1;
    ofstream out_f; // Output file

    // C++11 random generator
    random_device rd;
    mt19937 rg(rd()); // Random generator

    srand48(1);

    if (argc == 6) {
        num_test = atoi(argv[5]);
    } else if (argc != 5) {
        cout<< "Usage: " << argv[0] << " <number_of_nodes_in_graph> <number_of_pairs> <box_width> <box_height> <number of tests (optional)>" << endl;
        exit(1);
    }

    n = atoi(argv[1]);
    int kPairs = atoi(argv[2]);
    box_width = atof(argv[3]);
    box_height = atof(argv[4]);

    // Check kPairs
    if (0 >= kPairs || kPairs > n * (n-1)) {
        cerr << argv[0] << ": invalid number of pairs" << endl;
        exit(2);
    }

    // Generate k testcases
    for (int k = 0; k < num_test; k++) {
        ListGraph g;  // graph declaration
        NodeName vname(g);  // name of graph nodes
        ListGraph::NodeMap<double> px(g), py(g);  // xy-coodinates for each node
        EdgeWeight weight(g);   // edge weights
        ListGraph::EdgeMap<double> lat(g), cap(g), cost(g);
        char fname[100];
        vector< pair<Node,Node> > avail_p;                  // Available pairs
        vector< tuple<Node,Node,double,double> > sel_p;     // Selected pairs
        double sol = 0;

        // Make sure vectors are empty (probably unnecessary)
        avail_p.clear();
        sel_p.clear();

        GenerateTriangulatedListGraph(g, vname, px, py, weight, n, box_width, box_height);

        // Copy weight to cost map and generate latencies and capacities
        for (EdgeIt e(g); e!= INVALID; ++e) {
            cost[e] = weight[e];
            lat[e] = 10*drand48();
            cap[e] = 20*drand48();
        }

        // Generate all possible pairs
        int pair_n = 0;
        for (NodeIt u(g); u != INVALID; ++u) {
            for (NodeIt v(g); v != INVALID; ++v) {
                if (u == v) // Vertices must be different
                    continue;

                avail_p.push_back(make_pair(u, v));
                pair_n++;
            }
        }

        // Shuffle available pairs vector
        shuffle(avail_p.begin(), avail_p.end(), rg);

        // Set Dijkstra class and generate pairs
        for (int i = 0; i < kPairs; i++) {
            double Tmax, q;

            // p = (u,v)
            auto p = avail_p.back();
            avail_p.pop_back();

#ifdef KNOWN_SOL
            /* Find path from u to v, and set Tmax and q as:
             *  - Tmax: sum of all latencies + a random value from 0 to X
             *  - q: minimum capacity on path divided by (kPairs + a random value from 0 to X)
             */
            Dijkstra<ListGraph,ListGraph::EdgeMap<double> > dij(g, cost);
            ListGraph::NodeMap<ListGraph::Arc> pred(g);
            double csum; // cost sum

            // Get path from u to v
            dij.init();
            dij.predMap(pred);
            dij.addSource(p.first);
            dij.start(p.second);
            if (dij.reached(p.second)) {
                vector<Edge> path;

                // Compute q and Tmax
                Tmax = 0;
                q = DBL_MAX;
                csum = 0;

                ListGraph::Arc a = pred[p.second];
                for (;;) {
                    // Sum of latencies and costs
                    Tmax += lat[a];
                    csum += cost[a];

                    // Minimum capacity
                    if (q > cap[a])
                        q = cap[a];

                    // Check if found source node
                    if (g.source(a) == p.first)
                        break;
                    else
                        a = pred[g.source(a)];
                }

                // Finish Tmax and q generation
                Tmax += X*drand48();
                q = q / (kPairs + X*drand48());
                sol += csum*q;
            } else {
                cerr << "wtf no path!" << endl;
                i--;
            }
#else
            Tmax = 50*drand48();
            q    = 15*drand48();

            if (!(g.valid(p.first) && g.valid(p.second)))
                cerr << "invalid nodes!\n";

            Edge e = g.addEdge(p.first, p.second);

            cap[e] = q + 1;
            lat[e] = Tmax - 1;
            cost[e] = 1e3;

            if (!g.valid(e))
                cerr << "invalid edge!\n";
#endif /* KNOWN_SOL */

            sel_p.push_back(make_tuple(p.first, p.second, q, Tmax));
        }

        /**** Print testcase (to file) ****/

        int nV = countNodes(g);
        int nA = countEdges(g);

        // Name format: input[_k]_n<#nodes>_p<#pairs>_t<test#>.in (k is present when the solution is known)
#ifdef KNOWN_SOL
        sprintf(fname, "input_k_n%d_p%d_t%d.in", n, kPairs, k+1);
#else
        sprintf(fname, "input_n%d_p%d_t%d.in", n, kPairs, k+1);
#endif /* KNOWN_SOL */
        out_f.open(fname);

        // Print to file
        out_f << nA << " " << nV << " " << kPairs << endl;
        for (NodeIt v(g); v != INVALID; ++v)
            out_f << vname[v] << " " << px[v] << " " << py[v] << endl;
        for (EdgeIt e(g); e != INVALID; ++e)
            //       u                       v                       latency          capacity         cost
            out_f << vname[g.u(e)] << " " << vname[g.v(e)] << " " << lat[e] << " " << cap[e] << " " << cost[e] << endl;
        for (auto p = sel_p.begin(); p != sel_p.end(); ++p) {
            Node s, t;
            double q, Tmax;

            // Unpack tuple
            tie(s, t, q, Tmax) = *p;

            //       s                  t                  q           Tmax
            out_f << vname[s] << " " << vname[t] << " " << q << " " << Tmax << endl;
        }
        out_f.close();

#ifdef KNOWN_SOL
        // Save expected solution
        cout << fname << ": " << sol << endl;
        sprintf(fname, "inputk_n%d_p%d_t%d.sol", n, kPairs, k+1);
        out_f.open(fname);
        out_f << sol << endl;
        out_f.close();
#endif /* KNOWN_SOL */

#ifdef GRAPHVIZ
        //  Save Graphviz file
#ifdef KNOWN_SOL
        sprintf(fname, "input_k_n%d_p%d_t%d.dot", n, kPairs, k+1);
#else
        sprintf(fname, "input_n%d_p%d_t%d.dot", n, kPairs, k+1);
#endif /* KNOWN_SOL */
        out_f.open(fname);
        string line;
        out_f << "graph g {\n";
        out_f << "\tnode [shape = \"circle\"];\n";
        for (NodeIt v(g); v != INVALID; ++v) {
            line = "\t";
            line += vname[v].c_str();
            line += ";\n";
            out_f << line;
        }
        for (EdgeIt e(g); e != INVALID; ++e) {
            char d[100], l[100], c[100];
            sprintf(d, "%.1e", cap[e]); // Capacity
            sprintf(l, "%.1e", lat[e]); // Latency
            sprintf(c, "%.1e", cost[e]); // Cost

            if (!g.valid(e)) {
                cerr << "invalid edge!\n";
                continue;
            }

            line = "\t";
            line += vname[g.u(e)].c_str() ;
            line += "  -- ";
            line += vname[g.v(e)].c_str();

            line += " [label = \"";

            line += "d = ";
            line += d;
            line += "\\n";

            line += "l = ";
            line += l;
            line += "\\n";

            line += "c = ";
            line += c;
            line += "\\n";

            line += "\" ];\n";
            out_f << line;
        }
        out_f << "}\n";
        out_f.close();
#endif /* GRAPHVIZ */
    }

    return 0;
}


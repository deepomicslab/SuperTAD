//
// Created by MSI-PC on 2019/12/14.
//

#include "compare.h"
#include "inputAndOutput.h"


Comparator::Comparator(std::string path1, std::string path2)
{
    init(path1, path2);
}


Comparator::~Comparator()
{
    for (int i=0; i<_n; i++) {
        delete _graph[i];
//        delete _rGraph[i];
    }
    delete _graph;
//    delete _rGraph;
}



void Comparator::init(std::string path1, std::string path2)
{
    SuperTAD::Reader::readBoundariesIntoGraph(path1, path2, _boundaries1, _boundaries2, _graph);
    _n1 = _boundaries1.size();
    _n2 = _boundaries2.size();
    _n = _n1 + _n2 + 2;
//    utils::print2Darray(_graph, _n, _n);
//    _rGraph = new int *[_n];
//    for (int i=0; i<_n; i++)
//        _rGraph[i] = new int [_n]{};
}

void Comparator::execute()
{
    std::vector<std::vector<double>> costMatrix;
    for (int i=0; i<_n1; i++) {
        std::vector<double> tmp;
        for (int j=0; j<_n2; j++) {
            double d = -(double)_graph[i+1][_n1+1+j];
            if (d==0)
                d == std::numeric_limits<double>::infinity();
            tmp.emplace_back(d);
        }
        costMatrix.emplace_back(tmp);
    }
//    printf("costMatrix:\n");
//    for (int i=0; i<costMatrix.size(); i++) {
//        for (int j=0; j<costMatrix[i].size(); j++) {
//            std::cout << costMatrix[i][j] << "\t";
//        }
//        std::cout << "\n";
//    }

    HungarianAlgorithm HungAlgo;
    vector<int> assignment;

    double cost = HungAlgo.Solve(costMatrix, assignment);
//    for (unsigned int x = 0; x < costMatrix.size(); x++)
//        std::cout << x << "," << assignment[x] << "\t";
//    std::cout << "\ncost: " << cost << std::endl;

    double sumX = 0., sumY = 0., sumIntersection = -2 * cost;
    for (int i=0; i<_n1; i++)
        sumX += _boundaries1[i].size;
    for (int i=0; i<_n2; i++)
        sumY += _boundaries2[i].size;
//    for (int i=0; i<_n1; i++) {
//        if (assignment[i] >= 0) {
//            printf("assignment[%d]=%d, cost=%f\n", i, assignment[i], costMatrix[i][assignment[i]]);
//            sumIntersection += costMatrix[i][assignment[i]];
//            printf("sumIntersection=%f\n", sumIntersection);
//        }
//    }
//    sumIntersection = -2*sumIntersection;
    double overlappingRatio = sumIntersection / (sumX + sumY);
//    printf("sumX=%f, sumY=%f, sumInter=%f\n", sumX, sumY, sumIntersection);
    printf("overlapping ratio for given two TAD sets is %f\n", overlappingRatio);
}


//void Comparator::execute()
//{
//    fordFulkerson(_graph, _rGraph, 0, _n-1);
//    int sumX = 0, sumY = 0, sumIntersection = 0;
//    for (int i=0; i<_n1; i++)
//        sumX += _boundaries1[i].size;
//    for (int i=0; i<_n2; i++)
//        sumY += _boundaries2[i].size;
//    for (int i=1; i<_n1+1; i++) {
//        for (int j=_n1+1; j<_n-1; j++) {
//            if (_graph[i][j]>0 && _rGraph[i][j]==0)
//                sumIntersection += _graph[i][j];
//        }
//    }
//    sumIntersection *= 2;
//    double overlappingRatio = (double) sumIntersection / (double) (sumX + sumY);
//    printf("sumX=%d, sumY=%d, sumInter=%d\n", sumX, sumY, sumIntersection);
//    printf("overlapping ratio for given two TAD sets is %f\n", overlappingRatio);
//}


///* Returns true if there is a path from source 's' to sink 't' in
//  residual graph. Also fills parent[] to store the path */
//bool Comparator::bfs(int **rGraph, int s, int t, int *parent)
//{
//    // Create a visited array and mark all vertices as not visited
//    bool *visited = new bool [_n]{};
//    memset(visited, 0, sizeof(visited));
//
//    // Create a queue, enqueue source vertex and mark source vertex
//    // as visited
//    std::queue<int> q;
//    q.push(s);
//    visited[s] = true;
//    parent[s] = -1;
//
//    // Standard BFS Loop
//    while (!q.empty())
//    {
//        int u = q.front();
//        q.pop();
//
//        for (int v=0; v<_n; v++)
//        {
//            if (visited[v]==false && rGraph[u][v] > 0)
//            {
//                q.push(v);
//                parent[v] = u;
//                visited[v] = true;
//            }
//        }
//    }
//
//    // If we reached sink in BFS starting from source, then return true, else false
//    return (visited[t] == true);
//}
//
//// Returns the maximum flow from s to t in the given graph
//int Comparator::fordFulkerson(int **graph, int **rGraph, int s, int t)
//{
//    int u, v;
//
//    // Create a residual graph and fill the residual graph with
//    // given capacities in the original graph as residual capacities
//    // in residual graph
//
////    int rGraph[V][V]; // Residual graph where rGraph[i][j] indicates
//
//    // residual capacity of edge from i to j (if there
//    // is an edge. If rGraph[i][j] is 0, then there is not)
//    for (u = 0; u < _n; u++)
//        for (v = 0; v < _n; v++)
//            rGraph[u][v] = graph[u][v];
//
//    utils::print2Darray(rGraph, _n, _n);
//
//    int *parent = new int [_n]{};  // This array is filled by BFS and to store path
//
//    int max_flow = 0;  // There is no flow initially
//
//    // Augment the flow while tere is path from source to sink
//    while (bfs(rGraph, s, t, parent))
//    {
//        // Find minimum residual capacity of the edges along the
//        // path filled by BFS. Or we can say find the maximum flow
//        // through the path found.
//        int path_flow = std::numeric_limits<int>::max();
//        for (v=t; v!=s; v=parent[v])
//        {
//            u = parent[v];
//            path_flow = std::min(path_flow, rGraph[u][v]);
//        }
//        printf("path_flow=%d\n", path_flow);
//
//        // update residual capacities of the edges and reverse edges
//        // along the path
//        for (v=t; v != s; v=parent[v])
//        {
//            u = parent[v];
//            rGraph[u][v] -= path_flow;
//            rGraph[v][u] += path_flow;
//        }
//
//        // Add path flow to overall flow
//        max_flow += path_flow;
//
//        utils::print2Darray(rGraph, _n, _n);
//    }
//
//    // Return the overall flow
//    return max_flow;
//}

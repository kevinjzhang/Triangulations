#ifndef SIMP_INFO_H
#define SIMP_INFO_H

#include <vector>
#include <algorithm>
#include <iostream>

#include<triangulation/dim3.h>
#include<triangulation/example3.h>
#include<triangulation/detail/triangulation.h>
#include<triangulation/detail/facenumbering.h>

using namespace regina;

// Class assumes triangulation->size() ^ 2 < INT32_MAX
class SimplexInfo {
    private:
        int label;
        std::vector<int> vertexLabel; //Vertex degree = 40 for one-vertex triangulations
        std::vector<int> edgeLabel;
        std::vector<int> edgeCombLabel; //Edge labels with opposite edges combined into a single edge
        std::vector<std::vector<int>> incidentEdgeDegrees;

        int factorial (int n) {
            if (n == 1) {
                return 1;
            }
            return n * factorial(n - 1);
        }

        // <= function for sorted vectors
        static bool compArr(const std::vector<int>& v1, const std::vector<int>& v2) {
            for (int i = 0; i < v1.size(); i++) {
                if (v1[i] > v2[i]) {
                    return false;
                } else if (v1[i] < v2[i]) {
                    return true;
                }
            }            
            return true;
        }

        //A <= rank function for a vertex ordering
       bool compVertex(int i, int j) {
            if (vertexLabel[i] < vertexLabel[j]) {
                return true;
            } else if (vertexLabel[i] > vertexLabel[j]) {
                return false;
            }
            return compArr(incidentEdgeDegrees[i], incidentEdgeDegrees[j]);
        }
        
    public:
        void disp() {
            std::cout << "Vertices:" << std::endl;
            for (int i = 0; i < vertexLabel.size(); i++) {
                std::cout << vertexLabel[i] << " ";
            }
            std::cout << std::endl;
            std::cout << "Incident Edges:" << std::endl;
            for (int i = 0; i < incidentEdgeDegrees.size(); i++) {
                for (int j = 0; j < incidentEdgeDegrees[i].size(); j++) {
                    std::cout << incidentEdgeDegrees[i][j] << " ";
                }
                std::cout << std::endl;
            }
        }

        //Returns all valid permutations 
        std::vector<Perm<4>::Index> getAllPerms() {
            std::vector<Perm<4>::Index> ans;
            for (Perm<4>::Index perm = 0; perm < Perm<4>::nPerms; ++perm) {
                for (int i = 1; i < 4; i++) {
                    //Image or preimage?
                    if(!compVertex(Perm<4>::atIndex(perm)[i - 1], Perm<4>::atIndex(perm)[i])) {
                        break;
                    } 
                }
                ans.push_back(perm);
            }
            return ans;
        }

        SimplexInfo(Simplex<3>* tetrahedra, int size) {
            //Get a labelling for vertices and edges (this is computed non-local to the tetrahedra)
            for (int i = 0; i < 4; i++) {
                vertexLabel.push_back(tetrahedra->face<0>(i)->degree());
            }
            for (int i = 0; i < 6; i++) {
                edgeLabel.push_back(tetrahedra->face<1>(i)->degree());
            }
            //Get a ranking function for vertices based on edges
            incidentEdgeDegrees = std::vector<std::vector<int>>(6, std::vector<int>());
            for (int i = 0; i < 6; i++) {
                incidentEdgeDegrees[Edge<3>::edgeVertex[i][0]].push_back(edgeLabel[i]);
                incidentEdgeDegrees[Edge<3>::edgeVertex[i][1]].push_back(edgeLabel[i]);
            }
            //Sort incident edge degrees
            for (int i = 0; i < 4; i++) {
                std::sort(incidentEdgeDegrees[i].begin(), incidentEdgeDegrees[i].end());
            }
            //Use the information i and 5-i are opposite
            for (int i = 0; i < 3; i++) {
                int hi = std::max(edgeLabel[i], edgeLabel[5 - i]);
                int lo = std::min(edgeLabel[i], edgeLabel[5 - i]);
                edgeCombLabel.push_back(lo * size + hi);
            }
            //Sort combined edge labels
            std::sort(edgeCombLabel.begin(), edgeCombLabel.end());
        }

        int numOrderings() {
            std::vector<int> ordering;
            for (int i = 0; i < 4; i++) {
                ordering.push_back(i);
            }
            std::sort(ordering.begin(), ordering.end(), [this](int i, int j) { 
                return this->compVertex(i, j); });

            //A partition based on equality
            std::vector<int> partition;
            int prev = 0;
            int current = 1;
            //current > prev means the run has ended, current == prev otherwise
            for (int i = 1; i < ordering.size(); i++) {
                if (compVertex(ordering[current], ordering[prev])) {
                    current++;
                } else {
                    partition.push_back(current);
                    prev = i;
                    current = 1;
                }
            }
            partition.push_back(current);
            int combs = 1;
            for (int num : partition) {
                combs *= factorial(num);
            }
            return combs;
        }

        bool operator <(const SimplexInfo & other) {
            return compArr(edgeCombLabel, other.edgeCombLabel); 
        }

        bool operator ==(const SimplexInfo & other) {
            return (edgeLabel == other.edgeLabel) && (vertexLabel == other.vertexLabel);
        }
};

#endif
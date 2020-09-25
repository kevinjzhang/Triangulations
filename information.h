#ifndef SIMP_INFO_H
#define SIMP_INFO_H

#include <vector>
#include <algorithm>
#include <iostream>

#include<triangulation/dim3.h>
#include<triangulation/dim4.h>
#include<triangulation/example3.h>
#include<triangulation/example4.h>
#include<triangulation/detail/triangulation.h>
#include<triangulation/detail/facenumbering.h>

using namespace regina;

// Class assumes triangulation->size() ^ 2 < INT32_MAX
class SimplexInfo {
    private:
        int label;
        std::vector<int> vertexLabel; //Vertex degree 
        std::vector<int> sortedVertex; //Vertex labels sorted
        std::vector<int> edgeCombLabel; //Edge labels with opposite (ordered faces(4d) or unordered edges(3d))
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
        template <int dim>
        bool compVertex(int i, int j) {
            if (vertexLabel[i] < vertexLabel[j]) {
                return true;
            } else if (vertexLabel[i] > vertexLabel[j]) {
                return false;
            }
            return compArr(incidentEdgeDegrees[i], incidentEdgeDegrees[j]);
        }
        
    public:
        int getLabel() {
            return label;
        }

        void disp() {
            std::cout << "Vertices:" << std::endl;
            for (int i = 0; i < vertexLabel.size(); i++) {
                std::cout << vertexLabel[i] << " ";
            }
            std::cout << std::endl;
            std::cout << "Edges:" << std::endl;
            for (int i = 0; i < edgeCombLabel.size(); i++) {
                std::cout << edgeCombLabel[i] << " ";
            }
            std::cout << std::endl;
        }

        //Returns all valid permutations 
        template <int dim>
        std::vector<int> getAllPerms() {
            std::vector<int> ans;
            for (int perm = 0; perm < Perm<dim + 1>::nPerms; ++perm) {
                for (int i = 1; i < dim + 1; i++) {
                    if(!compVertex<dim>(Perm<dim + 1>::atIndex(perm)[i - 1], Perm<dim + 1>::atIndex(perm)[i])) {
                        break;
                    } 
                }
                ans.push_back(perm);
            }
            return ans;
        }

        SimplexInfo(Simplex<3>* tetrahedra, int simpNum, int size) {
            label = simpNum;
            //Get a labelling for vertices and edges (this is computed non-local to the tetrahedra)
            for (int i = 0; i < 4; i++) {
                vertexLabel.push_back(tetrahedra->face<0>(i)->degree());
                sortedVertex.push_back(tetrahedra->face<0>(i)->degree());
            }
            //Get a ranking function for vertices based on edges
            incidentEdgeDegrees = std::vector<std::vector<int>>(4, std::vector<int>());
            for (int i = 0; i < 6; i++) {
                incidentEdgeDegrees[Edge<3>::edgeVertex[i][0]].push_back(tetrahedra->face<1>(i)->degree());
                incidentEdgeDegrees[Edge<3>::edgeVertex[i][1]].push_back(tetrahedra->face<1>(i)->degree());
            }
            //Sort vertex degrees
            std::sort(sortedVertex.begin(), sortedVertex.end());
            //Sort incident edge degrees
            for (int i = 0; i < 4; i++) {
                std::sort(incidentEdgeDegrees[i].begin(), incidentEdgeDegrees[i].end());
            }
            //Use the information i and 5-i are opposite
            for (int i = 0; i < 3; i++) {
                int hi = std::max(tetrahedra->face<1>(i)->degree(), tetrahedra->face<1>(5 - i)->degree());
                int lo = std::min(tetrahedra->face<1>(i)->degree(), tetrahedra->face<1>(5 - i)->degree());
                edgeCombLabel.push_back(lo * size + hi);
            }
            //Sort combined edge labels
            std::sort(edgeCombLabel.begin(), edgeCombLabel.end());
        }

        SimplexInfo(Simplex<4>* pentachora, int simpNum, int size) {
            label = simpNum;
            for (int i = 0; i < 5; i++) {
                vertexLabel.push_back(pentachora->face<0>(i)->degree());
                sortedVertex.push_back(pentachora->face<0>(i)->degree());
            }
            for (int i = 0; i < 10; i++) {
                edgeCombLabel.push_back(pentachora->face<1>(i)->degree() * size + pentachora->face<2>(i)->degree());
            }

            incidentEdgeDegrees = std::vector<std::vector<int>>(5, std::vector<int>());
            for (int i = 0; i < 10; i++) {
                incidentEdgeDegrees[Edge<4>::edgeVertex[i][0]].push_back(edgeCombLabel[i]);
                incidentEdgeDegrees[Edge<4>::edgeVertex[i][1]].push_back(edgeCombLabel[i]);
            }
            for (int i = 0; i < 5; i++) {
                std::sort(incidentEdgeDegrees[i].begin(), incidentEdgeDegrees[i].end());
            }
            //Sort vertex degrees
            std::sort(sortedVertex.begin(), sortedVertex.end());
            //Sort Edge labels
            std::sort(edgeCombLabel.begin(), edgeCombLabel.end());
        }

        template <int dim>
        int numOrderings() {
            std::vector<int> ordering;
            for (int i = 0; i < dim + 1; i++) {
                ordering.push_back(i);
            }
            std::sort(ordering.begin(), ordering.end(), [this](int i, int j) { 
                return this->compVertex<dim>(i, j); 
            });

            //A partition based on equality
            std::vector<int> partition;
            int prev = 0;
            int current = 1;
            //current > prev means the run has ended, current == prev otherwise
            for (int i = 1; i < ordering.size(); i++) {
                if (compVertex<dim>(ordering[current], ordering[prev])) {
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
            if (sortedVertex == other.sortedVertex) {
                return compArr(edgeCombLabel, other.edgeCombLabel); 
            } else {
                return compArr(sortedVertex, other.sortedVertex); 
            }
        }

        bool operator ==(const SimplexInfo & other) {
            return (edgeCombLabel == other.edgeCombLabel) && (sortedVertex == other.sortedVertex);
        }
};

#endif
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>
#include <type_traits>
#include <chrono> 

#include "isosig.h"
#include "search.h"
#include "information.h"

#include<triangulation/dim3.h>
#include<triangulation/example3.h>
#include<triangulation/detail/triangulation.h>
#include<surfaces/normalsurfaces.h>
#include<surfaces/normalcoords.h>
#include<surfaces/nsvectorquad.h>
#include<link/link.h>

#include <maths/ray.h>

// #define STAT
// #define CORRECTNESS

using namespace regina;

int main(int argc, char *argv[]) {
    std::string inFile = std::string(argv[1]);
    std::string outFile = std::string(argv[2]);
    std::ifstream in (inFile, std::ifstream::in);
    std::ofstream out;
    out.open(outFile);
    int number;
    in >> number;
#ifdef STAT
    std::unordered_map<int, int> hist;
    for(int x = 0; x < number; x++) {
        std::string name;
        in >> name;
        Triangulation<3>* triangulation = Triangulation<3>::fromIsoSig(name);
        //Properties is an array of information
        std::vector<SimplexInfo> properties;
        for (int i = 0; i < triangulation->size(); i++) {
            Simplex<3>* tetrahedra = triangulation->simplex(i);
            properties.emplace_back(SimplexInfo(tetrahedra, triangulation->size()));
        }
        std::sort(properties.begin(), properties.end());
        delete triangulation;

        //Iterate through and find shortest run
        int prev = 0;
        int current = 1;
        std::vector<int> partition;
        int combs = INT32_MAX;
        for (int i = 1; i < properties.size(); i++) {
            if (properties[prev] == properties[i]) {
                current++;
            } else {
                partition.push_back(current);
                prev = i;
                current = 1;
            }
        }
        partition.push_back(current);
        //Compute statistics on how many configurations need to be used
        int index = 0;
        for (auto num : partition) {
            int count = 0;
            for (int i = 0; i < num; i++) {
                count += properties[index].numOrderings();
                index++;
            }
            combs = std::min(count, combs);
        }
        hist[combs]++;
    }

    for(auto& it : hist) {
        out << it.first << " " << it.second << std::endl;
    }
#endif
#ifdef CORRECTNESS
    for(int x = 0; x < number; x++) {
        std::string name;
        in >> name;
        Triangulation<3>* triangulation = Triangulation<3>::fromIsoSig(name);
        Triangulation<3>* triCopy = Triangulation<3>::fromIsoSig(IsoSig::computeSignature(triangulation));
        if (name != triCopy->isoSig()) {
            out << "Case Error" << std::endl;
        }
        delete triangulation;
        delete triCopy;
    }
#else
    std::vector<Triangulation<3>*> triangulations;
    for(int x = 0; x < number; x++) {
        std::string name;
        in >> name;
        Triangulation<3>* triangulation = Triangulation<3>::fromIsoSig(name);
        triangulations.push_back(triangulation);
    }
    auto start = std::chrono::high_resolution_clock::now();
    for (auto tri : triangulations) {
        std::string newSig = IsoSig::computeSignature(tri);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    out << "New: " << duration.count() << std::endl;

    for (auto tri : triangulations) {
        std::string oldSig = tri->isoSig();
    }
    auto stop2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - stop);
    out << "Old: " << duration.count() << std::endl;
#endif
    out.close();
    return 0;
}




#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>
#include <type_traits>
#include <chrono> 
#include <mpi.h>

#include "isosig.h"
#include "search.h"
#include "searchParallel.h"
#include "information.h"

#include<triangulation/dim3.h>
#include<triangulation/dim4.h>
#include<triangulation/example3.h>
#include<triangulation/example4.h>
#include<triangulation/detail/triangulation.h>
#include<surfaces/normalsurfaces.h>
#include<surfaces/normalcoords.h>
#include<surfaces/nsvectorquad.h>
#include<link/link.h>

#include <maths/ray.h>

// #define STAT
// #define CORRECTNESS
// #define TIMING
#define SEARCH

template <int dim>
void verify_correctness(int number, std::ifstream& in,  std::ofstream& out) {
    std::vector<std::string> names;
    for(int x = 0; x < number; x++) {
        std::string name;
        in >> name;
        names.push_back(name);
    }
    #pragma omp parallel for
    for (int i = 0; i < names.size(); i++) {
        std::string name = names[i];
        Triangulation<dim>* triangulation = Triangulation<dim>::fromIsoSig(name);
        std::string newName = IsoSig::computeSignature(triangulation);
        //Reorder labels randomly
        std::unordered_set<std::string> s;
        s.insert(newName);
        for (int simp = 0; simp < triangulation->size(); ++simp) {
            for (int perm = 0; perm < Perm<dim + 1>::nPerms; ++perm) {
                std::string curr = IsoSig::isoSigFrom(triangulation, triangulation->simplex(simp)->index(),
                    Perm<dim + 1>::atIndex(perm), (Isomorphism<dim>*) nullptr);
                Triangulation<dim>* relabelled = Triangulation<dim>::fromIsoSig(curr);
                std::string check = IsoSig::computeSignature(relabelled);
                delete relabelled;
                s.insert(check);
            }
        }
        if (s.size() != 1) {
            out << name << std::endl;
        }
        delete triangulation;
    }    
}

template <int dim>
void output_stats(int number, std::ifstream& in,  std::ofstream& out) {
    std::unordered_map<int, int> hist;
    std::vector<std::string> names;
    for(int x = 0; x < number; x++) {
        std::string name;
        in >> name;
        names.push_back(name);
    }

    #pragma omp parallel for
    for(int x = 0; x < number; x++) {
        Triangulation<dim>* triangulation = Triangulation<dim>::fromIsoSig(names[x]);
        //Properties is an array of information
        std::vector<SimplexInfo> properties;
        for (int i = 0; i < triangulation->size(); i++) {
            Simplex<dim>* tetrahedra = triangulation->simplex(i);
            properties.emplace_back(SimplexInfo(tetrahedra, i, triangulation->size()));
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
                count += properties[index].numOrderings<dim>();
                index++;
            }
            combs = std::min(count, combs);
        }
        #pragma omp critical
        hist[combs]++;
    }

    for(auto& it : hist) {
        out << it.first << " " << it.second << std::endl;
    }
}

template <int dim>
void check_perf(int number, std::ifstream& in,  std::ofstream& out) {
    std::vector<Triangulation<dim>*> triangulations;
    for(int x = 0; x < number; x++) {
        std::string name;
        in >> name;
        Triangulation<dim>* triangulation = Triangulation<dim>::fromIsoSig(name);
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
}

int main(int argc, char *argv[]) {
    std::string inFile = std::string(argv[1]);
    std::string outFile = std::string(argv[2]);
    std::ifstream in (inFile, std::ifstream::in);
    std::ofstream out;
    out.open(outFile);
    int number;
    in >> number;
#ifdef STAT
    output_stats<4>(number, in, out);
#endif
#ifdef CORRECTNESS
    //4d correctness check
    verify_correctness<4>(number, in, out);
#endif
#ifdef TIMING
    check_perf<4>(number, in, out);
#endif
#ifdef SEARCH
    std::vector<std::string> names;
    int maxHeight;
    in >> maxHeight;
    for(int x = 0; x < number; x++) {
        std::string name;
        in >> name;
        //Names in new isoSig form
#ifdef REGINASIG
        names.emplace_back(name);
#else
        names.emplace_back(IsoSig::computeSignature(regina::Triangulation<3>::fromIsoSig(name)));
#endif
    }
    SearchParallel::searchExhaustiveParallel<3>(names, maxHeight);
    // Search::searchExhaustive<4>(names, maxHeight, 10000); 
#endif
    out.close();
    return 0;
}




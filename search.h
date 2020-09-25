#ifndef SEARCH_H
#define SEARCH_H
#include <string>
#include <vector>
#include <mutex>
#include <unordered_map>

#include <unistd.h>

#include <mpi.h>
#include <omp.h>

#include<triangulation/dim3.h>
#include<triangulation/example3.h>
#include<triangulation/detail/triangulation.h>
#include<triangulation/detail/isosig-impl.h>


/* Warning: Sigset and processqueue share lock called processqueue
*/

using namespace regina;
class Search {
    //Private variables are temporarily unimplemented features
private:
    Search();

    template <int dim, class T, class U>
    static void processNode(T& sigSet, U& processingQueue, int tLimit) {
        std::string sig;
        #pragma omp critical(processQueue)
        {
            if (processingQueue.size() > 0) {
                sig = processingQueue.front();
                processingQueue.pop();
            }
        }
        if (sig.size() == 0) {
            return;
        }
        Triangulation<dim>* t = sigSet[sig];
        std::vector<Triangulation<dim>*> adj = getPachnerMoves(t, tLimit);
        //Convert all to sigs and add to processingQueue + sigSet
        for (auto tri : adj) {
            std::string s = IsoSig::computeSignature(tri);
            #pragma omp critical(processQueue)
            {
                if (sigSet.count(s) == 0) { //New triangulation
                    sigSet[s] = tri;
                    processingQueue.push(s);
                } else { //Duplicate triangulation found
                    delete tri;
                }
            }
        }
        //Deleting triangulation occurs after it has been processed
        delete t;        
    }
public:   
    static std::vector<Triangulation<3>*> getPachnerMoves (Triangulation<3>* t, int tLimit) {
        std::vector<Triangulation<3>*> adj;
        //Get all copies of tetrahedra made using 3-2 moves
        for (int i = 0; i < t->countEdges(); i++) {
            if (t->pachner(t->edge(i), true, false)) {
                Triangulation<3>* alt = new Triangulation<3>(*t, false);
                alt->pachner(alt->edge(i), false, true);
                adj.push_back(alt);
            }
        }
        //Get all copies of tetrahedra made using 2-3 moves
        if (t->size() < tLimit) {
            for (int i = 0; i < t->countTriangles(); i++) {
                if (t->pachner(t->triangle(i), true, false)) {
                    Triangulation<3>* alt = new Triangulation<3>(*t, false);
                    alt->pachner(alt->triangle(i), false, true);
                    adj.push_back(alt);
                }           
            }     
        }
        return adj;
    }

    static std::vector<Triangulation<4>*> getPachnerMoves (Triangulation<4>* t, int tLimit) {
        std::vector<Triangulation<4>*> adj;
        //5-1 move
        for (int i = 0; i < t->countVertices(); i++) {
            if (t->pachner(t->vertex(i), true, false)) {
                Triangulation<4>* alt = new Triangulation<4>(*t, false);
                alt->pachner(alt->vertex(i), false, true);
                adj.push_back(alt);
            }
        }
        //4-2 move
        for (int i = 0; i < t->countEdges(); i++) {
            if (t->pachner(t->edge(i), true, false)) {
                Triangulation<4>* alt = new Triangulation<4>(*t, false);
                alt->pachner(alt->edge(i), false, true);
                adj.push_back(alt);
            }
        }        
        //3-3 move
        for (int i = 0; i < t->countTriangles(); i++) {
            if (t->pachner(t->triangle(i), true, false)) {
                Triangulation<4>* alt = new Triangulation<4>(*t, false);
                alt->pachner(alt->triangle(i), false, true);
                adj.push_back(alt);
            }
        }       
        //2-4 move
        if (t->size() + 2 <= tLimit) {
            for (int i = 0; i < t->countTetrahedra(); i++) {
                if (t->pachner(t->tetrahedron(i), true, false)) {
                    Triangulation<4>* alt = new Triangulation<4>(*t, false);
                    alt->pachner(alt->tetrahedron(i), false, true);
                    adj.push_back(alt);
                }
            }        
        }
        //1-5 move
        if (t->size() + 4 <= tLimit) {
            for (int i = 0; i < t->size(); i++) {
                Triangulation<4>* alt = new Triangulation<4>(*t, false);
                alt->pachner(alt->pentachoron(i), false, true);
                adj.push_back(alt);
            }        
        }
        return adj;        
    }

    template <int dim>
    static void searchExhaustive(std::vector<std::string> & start, int tLimit) {
        std::unordered_map<std::string, Triangulation<3>*> sigSet;
        //(MPI)Must receive into this queue when message received
        std::queue<std::string> processingQueue;
        for (auto name : start) {
            processingQueue.push(name);
            sigSet[name] = Triangulation<dim>::fromIsoSig(name);
        }
        #pragma omp parallel
        #pragma omp single
        {
            while(!processingQueue.empty()) {
                #pragma omp task
                {
                    processNode<dim>(sigSet, processingQueue, tLimit);
                }
            }
            #pragma omp taskwait
        }
        std::cout << sigSet.size() << std::endl;
    }
};
#endif
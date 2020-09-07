#ifndef SEARCH_H
#define SEARCH_H
#include <string>
#include <vector>
#include <unordered_map>

#include <omp.h>

#include<triangulation/dim3.h>
#include<triangulation/example3.h>
#include<triangulation/detail/triangulation.h>
#include<triangulation/detail/isosig-impl.h>

using namespace regina;
class Search {
    //Private variables are temporarily unimplemented features
private:
    //(MPI) Number of computers, required for distributed hashing
    int nComp;
    //(MPI) A queue for sending batches of data outwards
    std::vector<std::queue<std::string>> sendBatch;
    Search();

    template <class T, class U>
    static void processNode(T& sigSet, U& processingQueue, int tLimit) {
        std::string sig;
        #pragma omp critical 
        {
            if (processingQueue.size() > 0) {
                sig = processingQueue.top();
                processingQueue.pop();
                std::cout << sig << std::endl;
            }
        }
        if (sig.size() == 0) {
            return;
        }
        Triangulation<3>* t = sigSet[sig];
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
        //Convert all to sigs and add to processingQueue + sigSet
        for (auto tri : adj) {
        #ifdef REGINASIG
            std::string s = tri->isoSig();
        #else
            std::string s = IsoSig::computeSignature(tri);
        #endif
            #pragma omp critical 
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
    static void searchExhaustive(std::vector<std::string> & start, int tLimit) {
        //Put items in sigSet before putting in queue
        std::unordered_map<std::string, Triangulation<3>*> sigSet;
        auto comp = []( std::string & s1, std::string & s2 ) { return s2.size() < s1.size(); };
        //(Optimisation) Change processing queue to batch sort pq
        //(MPI)Must receive into this queue when message received
        std::priority_queue<std::string, std::vector<std::string>, decltype(comp)> processingQueue(comp);
        for (auto name : start) {
            processingQueue.push(name);
            sigSet[name] = Triangulation<3>::fromIsoSig(name);
        }
        #pragma omp parallel
        #pragma omp single
        {
            while(!processingQueue.empty()) {
                #pragma omp task
                {
                    processNode(sigSet, processingQueue, tLimit);
                }
            }
            #pragma omp taskwait
        }
        std::cout << sigSet.size() << std::endl;
    }
};
#endif
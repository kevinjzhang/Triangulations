#ifndef SEARCH_H
#define SEARCH_H
#include <string>
#include <vector>
#include <unordered_set>

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
public:
    static void searchExhaustive(std::string start, int tLimit) {
        std::unordered_set<std::string> sigSet;
        sigSet.insert(start);
        auto comp = []( std::string s1, std::string s2 ) { return s1.size() < s2.size(); };
        //(Optimisation) Change processing queue to batch sort pq
        //(MPI)Must receive into this queue when message received
        std::priority_queue<std::string, std::vector<std::string>, decltype(comp)> processingQueue(comp);
        processingQueue.push(start);
        while(!processingQueue.empty()) {
            std::string sig = processingQueue.top();
            processingQueue.pop();
            //cilk_spawn func(sig)
            //(Optimisation) Reconstruct only if necessary
            Triangulation<3>* t = Triangulation<3>::fromIsoSig(sig);
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
                std::string s = IsoSig::computeSignature(tri);
                if (sigSet.count(s) == 0) {
                    sigSet.insert(s);
                    processingQueue.push(s);
                    //(Optimisation) Reconstruct only if necessary
                    delete tri;
                } else {
                    //Use a try to catch errors if they occur
                    delete tri;
                }
            }
            delete t;
        }
    }
};
#endif
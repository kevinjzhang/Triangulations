#ifndef SEARCH_PARALLEL_H
#define SEARCH_PARALLEL_H
#include "search.h"

using namespace regina;
class SearchParallel {
    private:
    //(MPI) Batch size to decide to send
    static const int batchSize = 100;
    static const int maxLength = 100;
    static const int wait = 1000000;
    static const int tag = 0;
    static const int tag1 = 1;
    SearchParallel();

    static void debug_info(std::string s) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::cout << s << " Rank:" << rank << std::endl;
    }

    static void receiveStatus(std::vector<bool>& states, int rank, int nComp) {
        MPI_Status status;
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, tag1, MPI_COMM_WORLD, &flag, &status);
        while (flag) {
            int count;
            int val;
            MPI_Recv(&val, 1, MPI_INT, MPI_ANY_SOURCE, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            states[status.MPI_SOURCE] = (bool)val;
            MPI_Iprobe(MPI_ANY_SOURCE, tag1, MPI_COMM_WORLD, &flag, &status);
        }
    }

    template <int dim, class T, class U>
    static bool receive(T& sigSet, U& processingQueue) {
        MPI_Status status;
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &flag, &status);
        while (flag) {
            int count;
            MPI_Get_count(&status, MPI_CHAR, &count);
            char* buffer = (char*)malloc(sizeof(char) * count);
            MPI_Recv(buffer, count, MPI_CHAR, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::string recvSig = std::string(buffer);
            #pragma omp critical(sig) 
            {
                if (sigSet.count(recvSig) == 0) { 
                    sigSet.insert(recvSig);
                    processingQueue.push(recvSig);
                }
            }
            MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &flag, &status);
        }        
        return flag;
    }

    //Waits for a small amount of time and checks if there is a message from any source
    template <int dim, class T, class U>
    static bool check_status(T& sigSet, U& processingQueue, std::vector<bool>& states, int rank, int nComp) {
        bool res = false;    
        states[rank] = true;
        #pragma omp critical(communication)
        {
            receiveStatus(states, rank, nComp);
            //Make all processors aware of their state
            int size;
            int data = 1;
            MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &size );
            int bufsize = (size + MPI_BSEND_OVERHEAD) * nComp;
            char* b = (char*)malloc(bufsize);
            MPI_Buffer_attach(b, bufsize);
            for (int i = 0; i < nComp; i++) {
                if (i != rank) {
                    MPI_Bsend(&data, 1, MPI_INT, i, tag1, MPI_COMM_WORLD);
                }
            }
            MPI_Buffer_detach(&b, &bufsize);
            free(b);
        }
        for (bool state : states) {
            if (state == false) {
                res = true;
                break;
            }
        }
        return res;
    }

    //Process nodes function in parallel
    template <int dim, class T, class U>
    static void processNodeParallel(T& sigSet, U& processingQueue, int tLimit, std::vector<std::queue<std::string>>& sendBatch,
            std::vector<bool>& states, int rank, int nComp) {
        //MPI receive into queue -> number of times until empty
        #pragma omp critical(communication)
        {
            if (receive<dim>(sigSet, processingQueue)) {
                if (states[rank]) {
                    states[rank] = false;
                    int size;
                    int data = 0;
                    MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &size );
                    int bufsize = (size + MPI_BSEND_OVERHEAD) * nComp;
                    char* b = (char*)malloc(bufsize);
                    MPI_Buffer_attach(b, bufsize);
                    for (int i = 0; i < nComp; i++) {
                        if (i != rank) {
                            MPI_Bsend(&data, 1, MPI_INT, i, tag1, MPI_COMM_WORLD);
                        }
                    }
                    MPI_Buffer_detach(&b, &bufsize);
                    free(b);
                }
            }
        }
        std::string sig;
        #pragma omp critical(sig)
        {
            if (processingQueue.size() > 0) {
                sig = processingQueue.front();
                processingQueue.pop();
                std::cout << sig << std::endl;
            }
        }
        if (sig.size() == 0) {
            return;
        }
        Triangulation<dim>* t = Triangulation<dim>::fromIsoSig(sig);
        std::vector<Triangulation<dim>*> adj = Search::getPachnerMoves(t, tLimit);
        //Convert all to sigs and add to processingQueue + sigSet
        for (auto tri : adj) {
            queueSig<dim>(sigSet, processingQueue, tri, sendBatch, rank);
        }
        //Deleting triangulation occurs after it has been processed
        delete t;      
    }

    //General function for queuing signature (multiple machines)
    template <int dim, class T, class U>
    static void queueSig(T& sigSet, U& processingQueue, Triangulation<dim>* tri, std::vector<std::queue<std::string>>& sendBatch, int rank) {
        std::string s = IsoSig::computeSignature(tri);
        delete tri;
        int hash = std::hash<std::string>{}(s) % sendBatch.size();
        //Compute locally
        if (hash == rank) {
            #pragma omp critical(sig)
            {
                if (sigSet.count(s) == 0) { //New triangulation
                    sigSet.insert(s);
                    processingQueue.push(s);
                } 
            }
        //Send Externally
        } else {
            //Second condition for finishing processing
            #pragma omp critical(communication)
            {
                sendBatch[hash].push(s);
                if (!sendBatch[hash].empty()) {
                    int size;
                    MPI_Pack_size( maxLength, MPI_CHAR, MPI_COMM_WORLD, &size );
                    int bufsize = (size + MPI_BSEND_OVERHEAD) * sendBatch[hash].size();
                    char* b = (char*)malloc(bufsize);
                    MPI_Buffer_attach(b, bufsize);
                    while (sendBatch[hash].size() > 0) {
                        std::string item = sendBatch[hash].front();
                        sendBatch[hash].pop();
                        MPI_Bsend(&item[0], item.size() + 1, MPI_CHAR, hash, tag, MPI_COMM_WORLD);
                    }
                    MPI_Buffer_detach(&b, &bufsize);
                    free(b);
                }
            }
        }
    }
public:
    template <int dim>
    static void searchExhaustiveParallel(std::vector<std::string> & start, int tLimit) {
        MPI_Init(NULL, NULL);
        int nComp;
        int rank;
        MPI_Comm_size(MPI_COMM_WORLD, &nComp);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::vector<std::queue<std::string>> sendBatch(nComp);
        std::unordered_set<std::string> sigSet;
        std::vector<bool> states(nComp, false);
        //Main application here:
        //(MPI)Must receive into this queue when message received
        std::queue<std::string> processingQueue;
        //Slight unneeded overhead for now (recomputes signature)
        if (rank == 0) {
            for (auto name : start) {
                queueSig<dim>(sigSet, processingQueue, Triangulation<dim>::fromIsoSig(name), sendBatch, rank);
            }
        }
        #pragma omp parallel
        #pragma omp single
        {
            while (!processingQueue.empty() || check_status<dim>(sigSet, processingQueue, states, rank, nComp)) {
                #pragma omp task
                {
                    processNodeParallel<dim>(sigSet, processingQueue, tLimit, sendBatch, states, rank, nComp);
                }
            }
            #pragma omp taskwait
        }        
        //Individual sizes
        std::cout << sigSet.size() << " processor:" << rank << std::endl;
        //Gather sizes
        int* res;
        if (rank == 0) {
            res =  (int*) malloc(sizeof(int) * nComp);
        } 
        int count = sigSet.size();
        MPI_Gather(&count, 1, MPI_INT, res, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //Display cumulative size if root processor
        if (rank == 0) {
            int sum = 0;
            for (int i = 0; i < nComp; i++) {
                sum += res[i];
            }
            std::cout << "Cumulative:" << sum << std::endl;
        }
        //Finished
        MPI_Finalize();
        if (rank == 0) {
            free(res);
        }
    }
};

#endif
#ifndef SEARCH_PARALLEL_H
#define SEARCH_PARALLEL_H
#include "search.h"

using namespace regina;
class SearchParallel {
    private:
    //(MPI) Batch size to decide to send
    static const int batchSize = 100;
    static const int maxLength = 100;
    static const int wait = 50000;
    static const int tag = 0;
    SearchParallel();

    static void debug_info(std::string s) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::cout << s << " Rank:" << rank << std::endl;
    }

    template <class T, class U>
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
            std::cout << "Recv: " << recvSig << std::endl;
            #pragma omp critical(processQueue) 
            {
                if (sigSet.count(recvSig) == 0) { 
                    sigSet[recvSig] = Triangulation<3>::fromIsoSig(recvSig);
                    processingQueue.push(recvSig);
                }
            }
            MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &flag, &status);
        }        
        return flag;
    }

    //Waits for a small amount of time and checks if there is a message from any source
    template <class T, class U>
    static bool check_status(T& sigSet, U& processingQueue) {
        int flag;
        MPI_Status status;
        MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &flag, &status);
        usleep(wait);
        bool res;
        #pragma omp critical(receive)
        {
            res = receive(sigSet, processingQueue);
        }
        return res;
    }

    //Process nodes function in parallel
    template <class T, class U>
    static void processNodeParallel(T& sigSet, U& processingQueue, int tLimit, std::vector<std::queue<std::string>>& sendBatch, int rank) {
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
        Triangulation<3>* t = sigSet[sig];
        std::vector<Triangulation<3>*> adj = Search::getPachnerMoves(t, tLimit);
        //Convert all to sigs and add to processingQueue + sigSet
        for (auto tri : adj) {
            queueSig(sigSet, processingQueue, tri, sendBatch, rank);
        }
        //Deleting triangulation occurs after it has been processed
        delete t;      
        //MPI receive into queue -> number of time it is non-empty
        #pragma omp critical(receive)
        {
            receive(sigSet, processingQueue);
        }
    }

    //General function for queuing signature (multiple machines)
    template <class T, class U>
    static void queueSig(T& sigSet, U& processingQueue, Triangulation<3>* tri, std::vector<std::queue<std::string>>& sendBatch, int rank) {
        std::string s = IsoSig::computeSignature(tri);
        int hash = std::hash<std::string>{}(s) % sendBatch.size();
        //Compute locally
        if (hash == rank) {
            #pragma omp critical(processQueue)
            {
                if (sigSet.count(s) == 0) { //New triangulation
                    sigSet[s] = tri;
                    processingQueue.push(s);
                } else { //Duplicate triangulation found
                    delete tri;
                }
            }
        //Send Externally
        } else {
            delete tri;
            //Second condition for finishing processing
            #pragma omp critical(send)
            {
                sendBatch[hash].push(s);
                if (sendBatch[hash].size() > batchSize || processingQueue.size() < batchSize) {
                    int size;
                    MPI_Pack_size( maxLength, MPI_CHAR, MPI_COMM_WORLD, &size );
                    int bufsize = (size + MPI_BSEND_OVERHEAD) * sendBatch[hash].size();
                    char* b = (char*)malloc(bufsize);
                    MPI_Buffer_attach(b, bufsize);
                    while (sendBatch[hash].size() > 0) {
                        std::string item = sendBatch[hash].front();
                        std::cout << "Send: " << item << std::endl;
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
    static void searchExhaustiveParallel(std::vector<std::string> & start, int tLimit) {
        MPI_Init(NULL, NULL);
        int nComp;
        int rank;
        MPI_Comm_size(MPI_COMM_WORLD, &nComp);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::vector<std::queue<std::string>> sendBatch(nComp);
        std::unordered_map<std::string, Triangulation<3>*> sigSet;
        //Main application here:
        //(MPI)Must receive into this queue when message received
        std::queue<std::string> processingQueue;
        //Slight unneeded overhead for now (recomputes signature)
        for (auto name : start) {
            queueSig(sigSet, processingQueue, Triangulation<3>::fromIsoSig(name), sendBatch, rank);
        }
        #pragma omp parallel
        #pragma omp single
        {
            while(!processingQueue.empty() || check_status(sigSet, processingQueue)) {
                #pragma omp task
                {
                    processNodeParallel(sigSet, processingQueue, tLimit, sendBatch, rank);
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
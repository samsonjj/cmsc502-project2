
#include "tsp.cpp"

/**
 * Gives an approximate solution to the TSP. Seperates cities into "blocks" on a
 * grid. Uses MPI processes to solve a complete TSP for each block. Then connects the blocks
 * using a parallel stitching algorithm
 */

using namespace std;


MPI_Comm cartcomm;
int dims[2] = {0, 0};

const int CITIES_PER_BLOCK = 8;


int doThing(int rank, int nthreads) {

    int row = rank / dims[1] + 1;
    int col = rank % dims[1] + 1;

    vector<city> cities = generate_cities(dims[0], row, col, CITIES_PER_BLOCK);

    // Run the sequential tsp for this block
    thread_vars *vars = new thread_vars();
    vars->cities = cities;
    solution sol = startDynamicSolution(vars, vars->cities);


    // TODO
    // if col % 2^(i+1) == 0, this process must receive
    // else, send, then quit UNLESS we are the last block, which is special
    // last block must receive according to the masking algorithm below

    vector<city> cities_ordered_by_path;
    for(int i = 0; i < sol.path.size(); i++) {
        cities_ordered_by_path.push_back(vars->cities[sol.path[i]]);
    }

    int num = 1;
    if(dims[1] == 1) {

    }
    else {
        for(int i = 0; (1<<i) <= dims[1]; i++) {


            if(col % (1<<(i+1)) == 0) {

                int n;
                int src = col-(1<<i)-1+(row-1)*dims[1];
                MPI_Recv(&n, 1, MPI_INT, col-(1<<i) - 1 + (row-1)*dims[1], 4, cartcomm, MPI_STATUS_IGNORE);
                num += n;

                // This processor will receive the cities from 1 other block, in the order of their path, with any previously stitched cities removed
                // Will also receive another list of cities from that block, which is the complete stitched path so far
                // This processor will calculate which cities to stitch from this list, and this processors cities
                // It will then calculate new total_distance, and create a new




                /*
                // get cities from another block in order of path
                vector<city> cities_incoming_by_path;
                int numCities;
                MPI_Recv(&numCities, 1, MPI_INT, 0, 0, cartcomm, MPI_STATUS_IGNORE);
                float x_coordinates[numCities];
                float y_coordinates[numCities];
                int ids[numCities];
                cout << "recieved" << endl;
                MPI_Recv(x_coordinates, numCities, MPI_FLOAT, 0, 1, cartcomm, MPI_STATUS_IGNORE);
                MPI_Recv(y_coordinates, numCities, MPI_FLOAT, 0, 2, cartcomm, MPI_STATUS_IGNORE);
                MPI_Recv(ids, numCities, MPI_INT, 0, 3, cartcomm, MPI_STATUS_IGNORE);

                cout << "recieved" << endl;

                // Init cities vector for this one block
                for (int i = 0; i < numCities; i++) {
                city newCity;
                newCity.x = x_coordinates[i];
                newCity.y = y_coordinates[i];
                newCity.id = ids[i];
                cities_incoming_by_path.push_back(newCity);
                }

                float distance;
                MPI_Recv(&distance, 1, MPI_INT, MPI_ANY_SOURCE, 4, cartcomm, MPI_STATUS_IGNORE);


                int min_swap_cost = numeric_limits<int>::max();
                int min_ul;
                int min_vl;
                int min_ur;
                int min_vr;
                for(int a1 = 0; a1 < cities_incoming_by_path.size(); a1++) {
                for(int a2 = 0; a2 < cities_ordered_by_path; a2++) {

                int b1 = (a1+1) % cities_incoming_by_path.size();
                int b2 = (b1+1) % cities_ordered_by_path.size();

                city ul = cities_incoming_by_path[a1];
                city vl = cities_incoming_by_path[b1];
                city ur = cities_ordered_by_path[a2];
                city vr = cities_incoming_by_path[b2];

                int swap_cost = calcDistance(ul, vr) + calcDistance(vl, ur) - calcDistance(ul, vl) - calcDistance(ur, vr);
                if(swap_cost < min_swap_cost) {
                min_swap_cost = swap_cost;
                min_ul = ul;
                min_vl = vl;
                min_ur = ur;
                min_vr = vr;
                }
                int swap_cost = calcDistance(ul, ur) + calcDistance(vl, vr) - calcDistance(ul, vl) - calcDistance(ur, vr);
                if(swap_cost < min_swap_cost) {
                min_swap_cost = swap_cost;
                min_ul = ul;
                min_vl = vl;
                min_ur = ur;
                min_vr = vr;
                }
                }
                }
                */

            }
            else {

                if(col == dims[1]) {
                    int n;

                    unsigned int mask = 1;
                    while((col & (~mask)) == col) {
                        mask = (mask << 1) + 1;
                    }

                    while((col & (~mask)) != 0) {
                        int src = (col & (~mask)) - 1 + (row-1)*dims[1];
                        MPI_Recv(&n, 1, MPI_INT, (col & (~mask)) - 1 + (row-1)*dims[1], 4, cartcomm, MPI_STATUS_IGNORE);
                        unsigned int temp = (col & (~mask));
                        while((col & (~mask)) == temp) {
                            mask = (mask << 1) + 1;
                        }

                        num+=n;
                    }

                    cout << "row " << row << " NUM " << num << endl;
                    break;
                }
                else {
                    int dest = col + (1 << i) - 1 >= dims[1] - 1 ? dims[1] - 1 + (row - 1) * dims[1] :
                               (col - 1 + (1 << i)) + (row - 1) * dims[1];
                    MPI_Ssend(&num, 1, MPI_INT, dest, 4, cartcomm);

                    break;
                }
            }

        }
    }

    return 0;
}


int main(int argc, char* argv[]) {

    // Timing
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    MPI_Init(&argc, &argv);

    int periods[2] = {0,0};

    int rank, nthreads;
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);


    if(rank == 0) printf("number of processors: %d\n", nthreads);
    if(rank == 0) printf("dims: %d, %d\n", dims[0], dims[1]);


    // Create cartesian topology, store in cartcomm
    MPI_Dims_create(nthreads, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims,
                    periods, 1, &cartcomm);

    MPI_Comm_rank(cartcomm, &rank);
    MPI_Comm_size(cartcomm, &nthreads);

    doThing(rank, nthreads);

    MPI_Finalize();

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    uint64_t delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    printf("took %lu\n", delta_us);

    return 0;
}







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

void send_cities(vector<city> cities, int dest) {
    int numCities = cities.size();
    MPI_Send(&numCities, 1, MPI_INT, dest, 0, cartcomm);
    int x_coordinates[numCities];
    int y_coordinates[numCities];
    for(int i = 0; i < numCities; i++) {
        x_coordinates[i] = cities[i].x;
        y_coordinates[i] = cities[i].y;
    }
    MPI_Send(x_coordinates, numCities, MPI_FLOAT, dest, 1, cartcomm);
    MPI_Send(y_coordinates, numCities, MPI_FLOAT, dest, 2, cartcomm);
}

vector<city> receive_cities(int src) {
    int numCities;
    MPI_Recv(&numCities, 1, MPI_INT, src, 0, cartcomm, MPI_STATUS_IGNORE);

    float x_coordinates[numCities];
    float y_coordinates[numCities];

    MPI_Recv(x_coordinates, numCities, MPI_FLOAT, src, 1, cartcomm, MPI_STATUS_IGNORE);
    MPI_Recv(y_coordinates, numCities, MPI_FLOAT, src, 2, cartcomm, MPI_STATUS_IGNORE);

    vector<city> cities;
    for(int i = 0; i < numCities; i++) {
        city newCity;
        newCity.x = x_coordinates[i];
        newCity.y = y_coordinates[i];
        cities.push_back(newCity);
    }

    return cities;
}

vector<city> stitch_cities_row(vector<city> cities1, vector<city> cities2) {

    int min_swap_cost = numeric_limits<int>::max();
    int min_ul;
    int min_vl;
    int min_ur;
    int min_vr;

    for(int i = 0; i < cities1.size(); i++) {
        for(int j = 0; j < cities2.size(); j++) {

            int uli = i;
            int vli = (i+1) % cities1.size();
            int uri = j;
            int vri = (j+1) % cities2.size();

            city ul = cities1[uli];
            city vl = cities1[vli];
            city ur = cities2[uri];
            city vr = cities2[vri];


            int swap_cost = calcDistance(ul, ur) + calcDistance(vl, vr) - calcDistance(ul, vl) - calcDistance(ur, vr);
            if(swap_cost < min_swap_cost) {
                min_swap_cost = swap_cost;
                min_ul = uli;
                min_vl = vli;
                min_ur = uri;
                min_vr = vri;
            }

            // try other combination
            int temp = uri;
            uri = vri;
            vri = temp;

            swap_cost = calcDistance(ul, ur) + calcDistance(vl, vr) - calcDistance(ul, vl) - calcDistance(ur, vr);
            if(swap_cost < min_swap_cost) {
                min_swap_cost = swap_cost;
                min_ul = uli;
                min_vl = vli;
                min_ur = uri;
                min_vr = vri;
            }
        }
    }

    vector<city> stitched_cities;

    for (int i = 0; i <= min_ul; i++) {
        stitched_cities.push_back(cities1[i]);
    }
    if(min_ur < min_vr) {
        for (int i = min_ur; true; i = (i-1 + cities2.size()) % cities2.size()) {
            stitched_cities.push_back(cities2[i]);
            if(i == min_vr) break;
        }
    }
    else {
        for (int i = min_ur; true; i = (i+1) % cities2.size()) {
            stitched_cities.push_back(cities2[i]);
            if(i == min_vr) break;
        }
    }
    for (int i = min_vl; i < cities1.size(); i++) {
        stitched_cities.push_back(cities1[i]);
    }

    return stitched_cities;
}

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

                int src = col-(1<<i)-1+(row-1)*dims[1];
                vector<city> other_cities = receive_cities(src);

                cities_ordered_by_path = stitch_cities_row(other_cities, cities_ordered_by_path);

                if(row == 1) printf("Received cities, col %d, # of cities is %d\n", col, cities_ordered_by_path.size());


            }
            else {

                if(col == dims[1]) {
                    unsigned int mask = 1;
                    while((col & (~mask)) == col) {
                        mask = (mask << 1) + 1;
                    }

                    while((col & (~mask)) != 0) {
                        int src = (col & (~mask)) - 1 + (row-1)*dims[1];
                        vector<city> other_cities = receive_cities(src);
                        cities_ordered_by_path = stitch_cities_row(other_cities, cities_ordered_by_path);
                        if(row == 1) printf("Received cities, col %d, # of cities is %d\n", col, cities_ordered_by_path.size());

                        unsigned int temp = (col & (~mask));
                        while((col & (~mask)) == temp) {
                            mask = (mask << 1) + 1;
                        }
                    }

                    cout << "row " << row << " NUM " << cities_ordered_by_path.size() << endl;
                    break;
                }
                else {
                    int dest = col + (1 << i) - 1 >= dims[1] - 1 ? dims[1] - 1 + (row - 1) * dims[1] :
                               (col - 1 + (1 << i)) + (row - 1) * dims[1];

                    // SEND CITIES
                    send_cities(cities_ordered_by_path, dest);


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


    // Create cartesian topology, store in cartcomm
    MPI_Dims_create(nthreads, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims,
                    periods, 1, &cartcomm);

    MPI_Comm_rank(cartcomm, &rank);
    MPI_Comm_size(cartcomm, &nthreads);


    if(rank == 0) printf("number of processors: %d\n", nthreads);
    if(rank == 0) printf("dims: %d, %d\n", dims[0], dims[1]);


    doThing(rank, nthreads);

    MPI_Finalize();

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    uint64_t delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    //printf("took %lu\n", delta_us);

    return 0;
}






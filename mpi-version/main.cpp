/*
 * Jonathan Samson
 * CMSC 502-001-2018Fall
 * November 15, 2018
 */

#include "tsp.cpp"
#include "intersect.cpp"

/**
 * Gives an approximate solution to the TSP. Seperates cities into "blocks" on a
 * grid. Uses MPI processes to solve a complete TSP for each block. Then connects the blocks
 * using a parallel stitching algorithm
 */

using namespace std;


MPI_Comm cartcomm;
int dims[2] = {0, 0};

const int CITIES_PER_BLOCK = 8;


void send_cities(vector <city> cities, int dest) {
    int numCities = cities.size();
    MPI_Send(&numCities, 1, MPI_INT, dest, 0, cartcomm);
    float x_coordinates[numCities];
    float y_coordinates[numCities];
    for (int i = 0; i < numCities; i++) {
        x_coordinates[i] = cities[i].x;
        y_coordinates[i] = cities[i].y;
    }
    MPI_Send(x_coordinates, numCities, MPI_FLOAT, dest, 1, cartcomm);
    MPI_Send(y_coordinates, numCities, MPI_FLOAT, dest, 2, cartcomm);
}

vector <city> receive_cities(int src) {
    int numCities;
    MPI_Recv(&numCities, 1, MPI_INT, src, 0, cartcomm, MPI_STATUS_IGNORE);

    float x_coordinates[numCities];
    float y_coordinates[numCities];

    MPI_Recv(x_coordinates, numCities, MPI_FLOAT, src, 1, cartcomm, MPI_STATUS_IGNORE);
    MPI_Recv(y_coordinates, numCities, MPI_FLOAT, src, 2, cartcomm, MPI_STATUS_IGNORE);

    vector <city> cities;
    for (int i = 0; i < numCities; i++) {
        city newCity;
        newCity.x = x_coordinates[i];
        newCity.y = y_coordinates[i];
        cities.push_back(newCity);
    }

    return cities;
}


/*
 * stitch together cities1 and cities2, store the final distance in sol->distance
 */
vector <city> stitch_cities_row(vector <city> cities1, vector <city> cities2, float distance2, solution *sol) {

    // Store min swap_cost and relavent city indexes
    int min_swap_cost = numeric_limits<int>::max();
    int min_ul;
    int min_vl;
    int min_ur;
    int min_vr;

    for (int i = 0; i < cities1.size(); i++) {
        for (int j = 0; j < cities2.size(); j++) {

            int uli = i;
            int vli = (i + 1) % cities1.size();
            int uri = j;
            int vri = (j + 1) % cities2.size();

            city ul = cities1[uli];
            city vl = cities1[vli];
            city ur = cities2[uri];
            city vr = cities2[vri];


            float swap_cost = calcDistance(ul, ur) + calcDistance(vl, vr) - calcDistance(ul, vl) - calcDistance(ur, vr);
            if (swap_cost < min_swap_cost) {
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

            ul = cities1[uli];
            vl = cities1[vli];
            ur = cities2[uri];
            vr = cities2[vri];

            swap_cost = calcDistance(ul, ur) + calcDistance(vl, vr) - calcDistance(ul, vl) - calcDistance(ur, vr);
            if (swap_cost < min_swap_cost) {
                min_swap_cost = swap_cost;
                min_ul = uli;
                min_vl = vli;
                min_ur = uri;
                min_vr = vri;
            }
        }
    }

    if (doIntersect(cities1[min_ul], cities2[min_ur], cities1[min_vl], cities2[min_vr]))
        cout << "the stitch is inverted!";

    for (int i = 0; i < cities1.size(); i++) {
        if (i != min_ul && i != min_vl && (i + 1) % cities1.size() != min_ul && (i + 1) % cities1.size() != min_vl) {
            if (doIntersect(cities1[min_ul], cities2[min_ur], cities1[i], cities1[(i + 1) % cities1.size()]))
                cout << "stitch1 caused an inversion!" << endl;
            if (doIntersect(cities1[min_vl], cities2[min_vr], cities1[i], cities1[(i + 1) % cities1.size()]))
                cout << "stitch2 caused an inversion!" << endl;
        }
    }

    for (int i = 0; i < cities2.size(); i++) {
        if (i != min_ur && i != min_vr && (i + 1) % cities2.size() != min_ur && (i + 1) % cities2.size() != min_vr) {
            if (doIntersect(cities1[min_ul], cities2[min_ur], cities2[i], cities2[(i + 1) % cities2.size()]))
                cout << "stitch3 caused an inversion!" << endl;
            if (doIntersect(cities1[min_vl], cities2[min_vr], cities2[i], cities2[(i + 1) % cities2.size()]))
                cout << "stitch4 caused an inversion!" << endl;
        }
    }


    // Push cities to a new vector, which has them in the order of the new path
    vector <city> stitched_cities;

    int c1 = 0, c2 = 0, c3 = 0;

    for (int i = min_vl; true; i++) {
        stitched_cities.push_back(cities1[i % cities1.size()]);
        c1++;
        if (i % cities1.size() == min_ul) break;
    }
    if ((min_ur - min_vr + cities2.size()) % cities2.size() == cities2.size() - 1) {
        for (int i = min_ur; true; i = (i - 1 + cities2.size()) % cities2.size()) {
            stitched_cities.push_back(cities2[i]);
            c2++;
            if (i == min_vr) break;
        }
    } else {
        for (int i = min_ur; true; i = (i + 1) % cities2.size()) {
            stitched_cities.push_back(cities2[i]);
            c3++;
            if (i == min_vr) break;
        }
    }


    sol->distance += distance2 + min_swap_cost;

    return stitched_cities;
}


/*
 * This is the JUICE of this application
 * Lot of thing happen here
 *
 * FIRST: We generate cities based off of row and col (row and col are indexed starting at 1)
 * SECOND: We stitch by row, every processor besides the last should quit after this step
 * THIRD: We stitch by column, every processor besides the bottom right in the cartesian topology should quit
 * LASTLY: We print the output
 *
 * We check during stitching to see if there are inversions
 */
void compute_tsp(int rank) {

    int row = rank / dims[1] + 1;
    int col = rank % dims[1] + 1;

    vector <city> cities = generate_cities(dims[0], row, col, CITIES_PER_BLOCK);

    // Run the sequential tsp for this block
    thread_vars *vars = new thread_vars();
    vars->cities = cities;
    solution sol = startDynamicSolution(vars, vars->cities);


    vector <city> cities_ordered_by_path;
    for (int i = 0; i < sol.path.size(); i++) {
        cities_ordered_by_path.push_back(vars->cities[sol.path[i]]);
    }


    // STITCH BY ROW
    // if col % 2^(i+1) == 0, this process must receive
    // else, send, then quit UNLESS we are the last block, which is special
    // last block must receive according to the masking algorithm below
    if (dims[0] != 1) {

        for (int i = 0; (1 << i) <= dims[1]; i++) {


            if (col % (1 << (i + 1)) == 0) {
                // WE ARE RECEIVING

                int src = col - (1 << i) - 1 + (row - 1) * dims[1];
                float other_distance;
                vector <city> other_cities = receive_cities(src);
                MPI_Recv(&other_distance, 1, MPI_FLOAT, src, 0, cartcomm, MPI_STATUS_IGNORE);

                cities_ordered_by_path = stitch_cities_row(other_cities, cities_ordered_by_path, other_distance, &sol);
            } else {
                // PROBABLY SENDING, UNLESS WE ARE THE LAST BLOCK IN THE ROW
                if (col == dims[1]) {
                    unsigned int mask = 1;
                    while ((col & (~mask)) == col) {
                        mask = (mask << 1) + 1;
                    }

                    while ((col & (~mask)) != 0) {

                        int src = (col & (~mask)) - 1 + (row - 1) * dims[1];
                        vector <city> other_cities = receive_cities(src);
                        float other_distance;
                        MPI_Recv(&other_distance, 1, MPI_FLOAT, src, 0, cartcomm, MPI_STATUS_IGNORE);

                        cities_ordered_by_path = stitch_cities_row(other_cities, cities_ordered_by_path, other_distance,
                                                                   &sol);

                        unsigned int temp = (col & (~mask));
                        while ((col & (~mask)) == temp) {
                            mask = (mask << 1) + 1;
                        }
                    }

                    break;
                } else {
                    int dest = col + (1 << i) - 1 >= dims[1] - 1 ? dims[1] - 1 + (row - 1) * dims[1] :
                               (col - 1 + (1 << i)) + (row - 1) * dims[1];

                    // SEND CITIES
                    send_cities(cities_ordered_by_path, dest);
                    MPI_Send(&sol.distance, 1, MPI_FLOAT, dest, 0, cartcomm);

                    break;
                }
            }

        }
    }


    MPI_Barrier(cartcomm);

    // STITCH BY COLUMN
    // if col % 2^(i+1) == 0, this process must receive
    // else, send, then quit UNLESS we are the last block, which is special
    // last block must receive according to the masking algorithm below
    if (dims[0] != 1 && col == dims[1]) {

        for (int i = 0; (1 << i) <= dims[0]; i++) {


            if (row % (1 << (i + 1)) == 0) {
                // WE ARE RECEIVING

                //int src = col-1+(row-1-(1<<i))*dims[1];
                int src;
                int coords[] = {(row - 1 - (1 << i)), (col - 1)};
                MPI_Cart_rank(cartcomm, coords, &src);

                float other_distance;
                vector <city> other_cities;

                other_cities = receive_cities(src);
                MPI_Recv(&other_distance, 1, MPI_FLOAT, src, 0, cartcomm, MPI_STATUS_IGNORE);

                cities_ordered_by_path = stitch_cities_row(other_cities, cities_ordered_by_path, other_distance, &sol);
            } else {
                // PROBABLY SENDING, UNLESS WE ARE THE LAST BLOCK IN THE COLUMN
                if (row == dims[0]) {
                    unsigned int mask = 1;
                    while ((row & (~mask)) == row) {
                        mask = (mask << 1) + 1;
                    }

                    while ((row & (~mask)) != 0) {

                        int src = col - 1 + ((row & (~mask)) - 1) * dims[0];
                        vector <city> other_cities = receive_cities(src);
                        float other_distance;
                        MPI_Recv(&other_distance, 1, MPI_FLOAT, src, 0, cartcomm, MPI_STATUS_IGNORE);

                        cities_ordered_by_path = stitch_cities_row(other_cities, cities_ordered_by_path, other_distance,
                                                                   &sol);

                        unsigned int temp = (row & (~mask));
                        while ((row & (~mask)) == temp) {
                            mask = (mask << 1) + 1;
                        }
                    }

                    break;
                } else {
                    int dest = row + (1 << i) - 1 >= dims[0] - 1 ? dims[0] * dims[1] - 1 :
                               (col - 1) + (row + (1 << i) - 1) * dims[0];

                    // SEND CITIES
                    send_cities(cities_ordered_by_path, dest);
                    MPI_Send(&sol.distance, 1, MPI_FLOAT, dest, 0, cartcomm);

                    break;
                }
            }

        }
    }

    MPI_Barrier(cartcomm);

    if(row == dims[0] && col == dims[1]) {
        cout << endl << "Distance: " << sol.distance << endl << endl;

        printf("%12s %12s\n\n", "x", "y");

        for (int i = 0; i < cities_ordered_by_path.size(); i++) {
            printf("%12.2f %12.2f\n", cities_ordered_by_path[i].x, cities_ordered_by_path[i].y);
        }
    }

    return;
}


int main(int argc, char *argv[]) {

    srand();

    // Timing
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    MPI_Init(&argc, &argv);

    int periods[2] = {0, 0};

    int rank, nthreads;
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);


    // Create cartesian topology, store in cartcomm
    MPI_Dims_create(nthreads, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims,
                    periods, 1, &cartcomm);

    MPI_Comm_rank(cartcomm, &rank);
    MPI_Comm_size(cartcomm, &nthreads);


    if (rank == 0) printf("number of processors: %d\n", nthreads);
    if (rank == 0) printf("dims: %d, %d\n", dims[0], dims[1]);


    compute_tsp(rank);

    MPI_Finalize();

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    uint64_t delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    //printf("took %lu\n", delta_us);

    return 0;
}






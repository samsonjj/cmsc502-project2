/*
 * Jonathan Samson
 * CMSC 502-001-2018Fall
 * November 15, 2018
 */

#include "../mpi-version/tsp.cpp"
#include "../mpi-version/intersect.cpp"

/**
 * Gives an approximate solution to the TSP. Seperates cities into "blocks" on a
 * grid. Uses pthreads to solve a complete TSP for each block. Then connects the blocks
 * using a parallel stitching algorithm
 */

using namespace std;

const int CITIES_PER_BLOCK = 10;
int nthreads;
int b;

// Shared memory between threads
pthread_t **threadArray;
pthread_mutex_t **mutexArray;
vector <city> **all_cities;
solution **all_solutions;


/*
 * stitch together cities1 and cities2, store the final distance in sol->distance
 */
vector <city>
stitch_cities_row(vector <city> cities1, vector <city> cities2, float distance2, solution *sol, int fromRow,
                  int fromCol) {

    // Store min swap_cost and relavent city indexes
    int min_swap_cost = numeric_limits<int>::max();
    int min_ul;
    int min_vl;
    int min_ur;
    int min_vr;

    // Find cities to swap
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

//    UNCOMMENT TO TEST FOR INVERSIONS
//    if (doIntersect(cities1[min_ul], cities2[min_ur], cities1[min_vl], cities2[min_vr])) {
//        cout << "the stitch is inverted!" << from << "," << to << endl;
//    }
//
//    for (int i = 0; i < cities1.size(); i++) {
//        if (i != min_ul && i != min_vl && (i + 1) % cities1.size() != min_ul && (i + 1) % cities1.size() != min_vl) {
//            if (doIntersect(cities1[min_ul], cities2[min_ur], cities1[i], cities1[(i + 1) % cities1.size()]))
//                cout << "stitch1 caused an inversion! " << from << "," << to << endl;
//            if (doIntersect(cities1[min_vl], cities2[min_vr], cities1[i], cities1[(i + 1) % cities1.size()]))
//                cout << "stitch2 caused an inversion!" << from << "," << to << endl;
//        }
//    }
//
//    for (int i = 0; i < cities2.size(); i++) {
//        if (i != min_ur && i != min_vr && (i + 1) % cities2.size() != min_ur && (i + 1) % cities2.size() != min_vr) {
//            if (doIntersect(cities1[min_ul], cities2[min_ur], cities2[i], cities2[(i + 1) % cities2.size()])) {
//                cout << "stitch3 caused an inversion!" << from << "," << to << endl;
//                cout << cities1[min_ul].x << "," << cities1[min_ul].y << endl;
//                cout << cities2[min_ur].x << "," << cities2[min_ur].y << endl;
//                cout << cities2[i].x << "," << cities2[i].y << endl;
//                cout << cities2[(i + 1) % cities2.size()].x << "," << cities2[(i + 1) % cities2.size()].y << endl;
//            }
//            if (doIntersect(cities1[min_vl], cities2[min_vr], cities2[i], cities2[(i + 1) % cities2.size()]))
//                cout << "stitch4 caused an inversion!" << from << "," << to << endl;
//        }
//    }


    // Push cities to a new vector, which has them in the order of the new path
    vector <city> stitched_cities;
    for (int i = min_vl; true; i++) {
        stitched_cities.push_back(cities1[i % cities1.size()]);
        if (i % cities1.size() == min_ul) break;
    }
    if ((min_ur - min_vr + cities2.size()) % cities2.size() == cities2.size() - 1) {
        for (int i = min_ur; true; i = (i - 1 + cities2.size()) % cities2.size()) {
            stitched_cities.push_back(cities2[i]);
            if (i == min_vr) break;
        }
    } else {
        for (int i = min_ur; true; i = (i + 1) % cities2.size()) {
            stitched_cities.push_back(cities2[i]);
            if (i == min_vr) break;
        }
    }

    // Fix inversions (Yes, there are inversions)
    bool fixed = false;
    int count = 0;
    while (!fixed && count < 15) {
        fixed = true;
        for (int i = 0; i < stitched_cities.size(); i++) {
            for (int j = i + 2; j < stitched_cities.size(); j++) {
                if ((j + 1) % stitched_cities.size() == i) continue;
                if (doIntersect(stitched_cities[i], stitched_cities[i + 1], stitched_cities[j],
                                stitched_cities[(j + 1) % stitched_cities.size()])) {

                    vector <city> inverted_fixed;
                    for (int k = 0; k <= i; k++) {
                        inverted_fixed.push_back(stitched_cities[k]);
                    }

                    for (int k = j; k >= i + 1; k--) {
                        inverted_fixed.push_back(stitched_cities[k]);
                    }

                    for (int k = j + 1; k != 0 && k < stitched_cities.size(); k++) {
                        inverted_fixed.push_back(stitched_cities[k]);
                    }

                    stitched_cities = inverted_fixed;
                    fixed = false;
                }

            }
        }
        count++;
    }

    // Calculate new total distance, since inversion swapping can make it weird
    float sum = 0;
    for (int i = 0; i < stitched_cities.size(); i++) {
        sum += calcDistance(stitched_cities[i], stitched_cities[(i + 1) % stitched_cities.size()]);
    }

    sol->distance = sum;

    return stitched_cities;
}


void *methodForThreads(void *pointer) {

    int threadNum = *((int *) pointer);

    int row = threadNum / b + 1;
    int col = threadNum % b + 1;

    all_cities[row - 1][col - 1] = generate_cities(b, row, col, CITIES_PER_BLOCK);


    // Run the sequential tsp for this block
    thread_vars *vars = new thread_vars();
    vars->cities = all_cities[row - 1][col - 1];
    all_solutions[row - 1][col - 1] = startDynamicSolution(vars, vars->cities);


    // Get the full path of cities
    vector <city> cities_ordered_by_path;
    for (int i = 0; i < all_solutions[row - 1][col - 1].path.size(); i++) {
        cities_ordered_by_path.push_back(vars->cities[all_solutions[row - 1][col - 1].path[i]]);
    }

    // Save them in shared memory
    all_cities[row - 1][col - 1] = cities_ordered_by_path;




    // STITCH BY ROW
    // if col % 2^(i+1) == 0, this process must receive
    // else, send, then quit UNLESS we are the last block, which is special
    // last block must receive according to the masking algorithm below
    if (nthreads != 1) {

        for (int i = 0; (1 << i) <= b; i++) {


            if (col % (1 << (i + 1)) == 0) {
                // WE ARE RECEIVING

                int srcCol = col - (1 << i) - 1;
                int srcRow = row - 1;

                pthread_mutex_lock(&mutexArray[srcRow][srcCol]);
                float other_distance = all_solutions[srcRow][srcCol].distance;
                vector <city> other_cities = all_cities[srcRow][srcCol];

                vector <city> stitched = stitch_cities_row(other_cities, all_cities[row - 1][col - 1],
                                                           other_distance, &all_solutions[row - 1][col - 1], srcRow,
                                                           srcCol);
                all_cities[row - 1][col - 1] = stitched;
            } else {

                // PROBABLY SENDING, UNLESS WE ARE THE LAST BLOCK IN THE ROW

                if (col == b) {

                    unsigned int mask = 1;
                    while ((col & (~mask)) == col) {
                        mask = (mask << 1) + 1;
                    }

                    while ((col & (~mask)) != 0) {

                        //int src = (col & (~mask)) - 1 + (row - 1) * dims[1];
                        int srcRow = row - 1;
                        int srcCol = (col & (~mask)) - 1;


                        pthread_mutex_lock(&mutexArray[srcRow][srcCol]);
                        float other_distance = all_solutions[srcRow][srcCol].distance;
                        vector <city> other_cities = all_cities[srcRow][srcCol];


                        vector <city> stitched = stitch_cities_row(other_cities, all_cities[row - 1][col - 1],
                                                                   other_distance, &all_solutions[row - 1][col - 1],
                                                                   srcRow, srcCol);
                        all_cities[row - 1][col - 1] = stitched;

                        unsigned int temp = (col & (~mask));
                        while ((col & (~mask)) == temp) {
                            mask = (mask << 1) + 1;
                        }
                    }

                    break;
                } else {

                    // SEND CITIES, we send by unlocking the sender's mutex lock, this allows the receiver to continue by locking it again
                    pthread_mutex_unlock(&mutexArray[row - 1][col - 1]);

                    break;
                }
            }

        }




        // STITCH BY COLUMN
        // if col % 2^(i+1) == 0, this process must receive
        // else, send, then quit UNLESS we are the last block, which is special
        // last block must receive according to the masking algorithm below
        if (col == b && nthreads != 1) {

            for (int i = 0; (1 << i) <= b; i++) {


                if (row % (1 << (i + 1)) == 0) {
                    // WE ARE RECEIVING

                    int srcCol = col - 1;
                    int srcRow = row - (1 << i) - 1;

                    pthread_mutex_lock(&mutexArray[srcRow][srcCol]);
                    float other_distance = all_solutions[srcRow][srcCol].distance;
                    vector <city> other_cities = all_cities[srcRow][srcCol];


                    vector <city> stitched = stitch_cities_row(other_cities, all_cities[row - 1][col - 1],
                                                               other_distance, &all_solutions[row - 1][col - 1],
                                                               srcCol + 1, col);
                    all_cities[row - 1][col - 1] = stitched;

                } else {

                    // PROBABLY SENDING, UNLESS WE ARE THE LAST BLOCK IN THE COL

                    if (row == b) {

                        unsigned int mask = 1;
                        while ((row & (~mask)) == row) {
                            mask = (mask << 1) + 1;
                        }

                        while ((row & (~mask)) != 0) {

                            //int src = (col & (~mask)) - 1 + (row - 1) * dims[1];
                            int srcRow = (row & (~mask)) - 1;
                            int srcCol = col - 1;

                            pthread_mutex_lock(&mutexArray[srcRow][srcCol]);
                            float other_distance = all_solutions[srcRow][srcCol].distance;
                            vector <city> other_cities = all_cities[srcRow][srcCol];


                            vector <city> stitched = stitch_cities_row(other_cities, all_cities[row - 1][col - 1],
                                                                       other_distance, &all_solutions[row - 1][col - 1],
                                                                       srcCol + 1, col);
                            all_cities[row - 1][col - 1] = stitched;

                            unsigned int temp = (row & (~mask));
                            while ((row & (~mask)) == temp) {
                                mask = (mask << 1) + 1;
                            }
                        }

                        break;
                    } else {

                        // SEND CITIES
                        pthread_mutex_unlock(&mutexArray[row - 1][col - 1]);

                        break;
                    }
                }

            }
        }
    }

    // Print
    if (row == b && col == b) {
        cout << endl << "Distance: " << all_solutions[row - 1][col - 1].distance << endl << endl;
        cout << endl << "Num: " << all_cities[row - 1][col - 1].size() << endl;

        printf("%12s %12s\n\n", "x", "y");

        for (int i = 0; i < all_cities[row - 1][col - 1].size(); i++) {
            printf("%12.2f %12.2f\n", all_cities[row - 1][col - 1][i].x, all_cities[row - 1][col - 1][i].y);
        }
    }


    return 0;
}


// Enter a square number of threads as the first command line argument
int main(int argc, char *argv[]) {

    srand(time(NULL));

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);


    if (argc > 1 && atoi(argv[1]) > 0) nthreads = atoi(argv[1]);
    else {
        cout << "Enter number of threads: ";
        scanf("%d", &nthreads);
    }

    if ((int) (sqrt(nthreads)) * (int) (sqrt(nthreads)) != nthreads) {
        cout << "Number of threads has to be square" << endl;
        return 1;
    }

    b = (int) sqrt(nthreads);

    all_cities = new vector <city> *[b];
    all_solutions = new solution *[b];
    threadArray = new pthread_t *[b];
    mutexArray = new pthread_mutex_t *[b];

    for (int i = 0; i < b; i++) {


        threadArray[i] = new pthread_t[b];
        mutexArray[i] = new pthread_mutex_t[b];
        all_cities[i] = new vector<city>[b];
        all_solutions[i] = new solution[b];

        for (int j = 0; j < sqrt(nthreads); j++) {

            int *thread_num = (int *) malloc(sizeof(*thread_num));
            *thread_num = i * b + j;

            pthread_create(&threadArray[i][j], NULL, methodForThreads, (void *) thread_num);
            pthread_mutex_init(&mutexArray[i][j], NULL);
            pthread_mutex_lock(&mutexArray[i][j]);
        }
    }

    for (int i = 0; i < b; i++) {
        for (int j = 0; j < b; j++) {
            void *status;
            pthread_join(threadArray[i][j], &status);
        }
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    uint64_t delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    printf("took %f seconds\n", delta_us / 1000000.0);

    return 0;
}

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

const int CITIES_PER_BLOCK = 4;
int nthreads;
int b;

pthread_t **threadArray;
pthread_mutex_t **mutexArray;

vector <city> **all_cities;
solution **all_solutions;


/*
 * stitch together cities1 and cities2, store the final distance in sol->distance
 */
vector <city>
stitch_cities_row(vector <city> cities1, vector <city> cities2, float distance2, solution *sol, int fromRow, int fromCol) {

//    bool fixed = false;
//    while(!fixed) {
//        fixed = true;
//        for (int i = 0; i < cities1.size(); i++) {
//            for (int j = i + 2; j < cities1.size(); j++) {
//                if((j+1) % cities1.size() == i) continue;
//                if (doIntersect(cities1[i], cities1[i + 1], cities1[j],
//                                cities1[(j + 1) % cities1.size()])) {
//
//                    vector <city> inverted_fixed;
//                    for (int k = 0; k <= i; k++) {
//                        inverted_fixed.push_back(cities1[k]);
//                    }
//
//                    for (int k = j; k >= i + 1; k--) {
//                        inverted_fixed.push_back(cities1[k]);
//                    }
//
//                    for (int k = j + 1; k != 0 && k < cities1.size(); k++) {
//                        inverted_fixed.push_back(cities1[k]);
//                    }
//
//
//                    cout << cities1[i].x << "," << cities1[i].y << endl;
//                    cout << cities1[i+1].x << "," << cities1[i+1].y << endl;
//                    cout << cities1[j].x << "," << cities1[j].y << endl;
//                    cout << cities1[(j + 1) % cities1.size()].x << "," << cities1[(j + 1) % cities1.size()].y << endl;
//                    cout << cities1[0].x << "," << cities1[0].y << endl;
//                    cout << cities1[1].x << "," << cities1[1].y << endl;
//                    cout << cities1[2].x << "," << cities1[2].y << endl;
//                    cout << cities1[3].x << "," << cities1[3].y << endl;
//
//                    cities1 = inverted_fixed;
//                    fixed = false;
//                    cout << "Oh hell no" << cities1.size() << " " << i << " " << i+1 << " "<< j << " " << (j+1) % cities1.size() << endl;
//                    for(int p = 0; p < all_solutions[fromRow][fromCol].path.size(); p++) {
//                        cout << all_solutions[fromRow][fromCol].path[p] << ", ";
//                    }
//                }
//            }
//        }
//    }

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


    bool fixed = false;
    while(!fixed) {
        fixed = true;
        for (int i = 0; i < stitched_cities.size(); i++) {
            for (int j = i + 2; j < stitched_cities.size(); j++) {
                if((j+1) % stitched_cities.size() == i) continue;
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
    }


    sol->distance += distance2 + min_swap_cost;

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

    vector <city> cities_ordered_by_path;
    for (int i = 0; i < all_solutions[row-1][col-1].path.size(); i++) {
        cities_ordered_by_path.push_back(vars->cities[all_solutions[row - 1][col - 1].path[i]]);
    }

    all_cities[row - 1][col - 1] = cities_ordered_by_path;

//    if(row == 1 && col == 1) {
//        for(int i = 0; i < all_cities[row-1][col-1].size(); i++) {
//            cout << all_cities[row-1][col-1][i].x << ", " << all_cities[row-1][col-1][i].y << endl;
//        }
//    }



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

                if (row == 1) cout << col << " waiting to receive from " << srcCol + 1 << endl;

                pthread_mutex_lock(&mutexArray[srcRow][srcCol]);
                float other_distance = all_solutions[srcRow][srcCol].distance;
                vector <city> other_cities = all_cities[srcRow][srcCol];

                if (row == 1) cout << col << " received from " << srcCol + 1 << endl;

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

//                        cout << "about to lock at " << srcRow << "," << srcCol << " by " << threadNum << endl;
                        if (row == 1) cout << col << " waiting to receive from " << (srcCol + 1) << endl;

                        pthread_mutex_lock(&mutexArray[srcRow][srcCol]);
//                        cout << "locked at " << srcRow << "," << srcCol << "by " << threadNum << endl;
                        float other_distance = all_solutions[srcRow][srcCol].distance;
                        vector <city> other_cities = all_cities[srcRow][srcCol];
                        if (row == 1) cout << col << " received from " << srcCol + 1 << endl;


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
//                    int dest = col + (1 << i) - 1 >= dims[1] - 1 ? dims[1] - 1 + (row - 1) * dims[1] :
//                               (col - 1 + (1 << i)) + (row - 1) * dims[1];
//                    int destRow = row-1;
//                    int destCol = col + (1 << i) - 1 >= b ? b-1 : col + (1 << 1) - 1;
//
                    // SEND CITIES
                    pthread_mutex_unlock(&mutexArray[row - 1][col - 1]);
//                    cout << "sent! from " << threadNum << endl;

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

//                    if(row==1) cout << col << " waiting to receive from " << srcCol+1 << endl;
                    cout << row << " Waiting to receive from " << (srcRow + 1) << endl;
                    pthread_mutex_lock(&mutexArray[srcRow][srcCol]);
                    float other_distance = all_solutions[srcRow][srcCol].distance;
                    vector <city> other_cities = all_cities[srcRow][srcCol];
                    cout << row << " Received from " << srcRow + 1 << endl;

//                    if(row==1) cout << col << " received from " << srcCol+1 << endl;

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

//                        cout << "about to lock at " << srcRow << "," << srcCol << " by " << threadNum << endl;
                            cout << row << " Waiting to receive from " << (srcRow + 1) << endl;

                            pthread_mutex_lock(&mutexArray[srcRow][srcCol]);
//                        cout << "locked at " << srcRow << "," << srcCol << "by " << threadNum << endl;
                            float other_distance = all_solutions[srcRow][srcCol].distance;
                            vector <city> other_cities = all_cities[srcRow][srcCol];
                            cout << row << " Received from " << srcRow + 1 << endl;


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
//                    int dest = col + (1 << i) - 1 >= dims[1] - 1 ? dims[1] - 1 + (row - 1) * dims[1] :
//                               (col - 1 + (1 << i)) + (row - 1) * dims[1];
//                    int destRow = row-1;
//                    int destCol = col + (1 << i) - 1 >= b ? b-1 : col + (1 << 1) - 1;
//
                        // SEND CITIES
                        cout << "Sending from " << col << endl;
                        pthread_mutex_unlock(&mutexArray[row - 1][col - 1]);
//                    cout << "sent! from " << threadNum << endl;

                        break;
                    }
                }

            }
        }







//
//        if(col == b) {
//            cout << all_cities[row-1][col-1].size() << ", " << all_solutions[row-1][col-1].distance << endl;
//        }

    }

    if (row == b && col == b) {
        cout << endl << "Distance: " << all_solutions[row - 1][col - 1].distance << endl << endl;
        cout << endl << "Num: " << all_cities[row-1][col-1].size() << endl;

        printf("%12s %12s\n\n", "x", "y");

        for (int i = 0; i < all_cities[row - 1][col - 1].size(); i++) {
            printf("%12.2f %12.2f\n", all_cities[row - 1][col - 1][i].x, all_cities[row - 1][col - 1][i].y);
        }
    }


    return 0;
}


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

    cout << "DONE" << endl;

/*
    // create threads
    thread_vars *vars;
    vector<city> cities_by_id;
    // store the "solutions" here, each containing the distance, first_city, and last_city in the path
    solution **solutionArray = new solution*[blocks.size()];
    pthread_t **threadArray = new pthread_t*[GRID_LENGTH];
    // for city ids
    int count = 0;
    for(int i = 0; i < blocks.size(); i++) { solutionArray[i] = new solution[blocks[i].size()];
        threadArray[i] = new pthread_t[GRID_LENGTH];
        for(int j = 0; j < blocks[i].size(); j++) {
            for(int k = 0; k < blocks[i][j].size(); k++) {
                blocks[i][j][k].id = count;
                count++;
                cities_by_id.push_back(blocks[i][j][k]);
            }
            vector<city> cities = blocks[i][j];
            //if(cities.size() == 0) continue;
            thread_vars *vars = new thread_vars();
            vars->cities = cities;
            vars->solutionArray = solutionArray;
            vars->i = i;
            vars->j = j;
            pthread_create(&threadArray[i][j], NULL, findAndStoreSolution, (void *) vars);
        }
    }

    for(int i = 0; i < GRID_LENGTH; i++) {
        for(int j = 0; j < GRID_LENGTH; j++) {
            void *status;
            pthread_join(threadArray[i][j], &status);
        }
    }

    if(logLevel == 1) {
        cout << endl << "Cities: " << endl;
        for(int i = 0; i < cities_by_id.size(); i++) {
            cout << setw(10) << left << cities_by_id[i].id << setw(8) << cities_by_id[i].x << setw(8) << cities_by_id[i].y << endl;
        }
    }



    // use cities_by_id to go through cities
    // store what blocks are visited
    // cycle through the 2 blocks which are visited, but not fully connected, and find their closest non-visited city
    // then link with them, calculate distance and add to total
    solution blocksAtEnd[2];
    blocksAtEnd[0] = solutionArray[0][0];

    // store which blocks have not been visited
    // should start without the first block
    vector<solution> unvisited;
    for(int i = 0; i < GRID_LENGTH; i++) {
        for(int j = 0; j < GRID_LENGTH; j++) {
            if(i==0 && j==0) continue;
            unvisited.push_back(solutionArray[i][j]);
        }
    }


    float totalDistance = solutionArray[0][0].distance;
    // go for all unvisited blocks
    for(int i = 0; i < GRID_LENGTH*GRID_LENGTH-1; i++) {
        float min_dist = numeric_limits<float>::max();
        int min_block_index;
        int min_city;
        // decide which block we will start from
        solution fromBlock = blocksAtEnd[i%2];
        // decide which city we are starting from
        city fromCity = fromBlock.visited_cities==FIRST ? cities_by_id[fromBlock.last_city] : cities_by_id[fromBlock.first_city];
        // now go through all unvisited blocks and find the shortest distance to any first or last city on any of them
        for(int j = 0; j < unvisited.size(); j++) {
            city city1 = cities_by_id[unvisited[j].first_city];
            city city2 = cities_by_id[unvisited[j].last_city];
            float dist1 = calcDistance(fromCity.x, fromCity.y, city1.x, city1.y);
            float dist2 = calcDistance(fromCity.x, fromCity.y, city2.x, city2.y);
            if(dist1 >= min_dist && dist2 >= min_dist) {
                // neither of these is an estimated optimal route
                continue;
            }
            // at this point we know we must change previous minimum
            min_block_index = j;
            min_city = dist1<=dist2 ? FIRST : LAST;
            min_dist = min(dist1, dist2);
        }

        // now travel to minimum city
        // remove minimum block from unvisited
        // set block_at_end_* to min_block
        // set block_*_visited_city to min_city
        // add the distance onto the total
        totalDistance += min_dist;
        totalDistance += unvisited[min_block_index].distance;
        if(i == 0) {
            // first loop
            blocksAtEnd[1] = unvisited[min_block_index];
            blocksAtEnd[1].visited_cities = min_city;
            blocksAtEnd[0].visited_cities = FIRST;
        } else {
            blocksAtEnd[i%2] = unvisited[min_block_index];
            blocksAtEnd[i%2].visited_cities = min_city;
        }
        unvisited.erase(unvisited.begin() + min_block_index);


    }
    // dont forget to add the final closing distance at the end
    city fromCity = blocksAtEnd[0].visited_cities==FIRST ? cities_by_id[blocksAtEnd[0].last_city] :
        cities_by_id[blocksAtEnd[0].first_city];
    city toCity = blocksAtEnd[1].visited_cities==FIRST ? cities_by_id[blocksAtEnd[0].last_city] :
        cities_by_id[blocksAtEnd[0].first_city];
    float dist = calcDistance(fromCity.x, fromCity.y, toCity.x, toCity.y);
    totalDistance += dist;



    if(logLevel==1) {
        cout << endl << "block solutions: " << endl;
        for(int i = 0; i < blocks.size(); i++) {
            for(int j = 0; j < blocks[i].size(); j++) {

                cout << setw(4) << i << setw(4) << j << setw(10) << solutionArray[i][j].distance << endl;

            }
        }
    }

    if(logLevel==1) cout << endl << "Final estimated distance: ";
    cout << totalDistance << endl;



    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    uint64_t delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    printf("took %lu\n", delta_us);
    */
    return 0;
}

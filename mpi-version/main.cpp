
#include "tsp.cpp"

/**
 * Gives an approximate solution to the TSP. Seperates cities into "blocks" on a
 * grid. Uses pthread solve a partial TSP for each block. Then connects the blocks
 * using a stitching algorithm.
 */

using namespace std;

int GRID_LENGTH;
int logLevel;

MPI_Comm cartcomm;
int dims[2] = {0, 0};

/**
 * @return returns a vector of 'city' objects
 */
vector<city> readInCities() {

    vector<city> all_cities;

    ifstream file;
    file.open("cities.data");
    if (!file.is_open()) {
        return all_cities;
    }

    string x_pos, y_pos;
    while(file >> x_pos >> y_pos) {
        city newCity;
        newCity.x = atof(x_pos.c_str());
        newCity.y = atof(y_pos.c_str());
        all_cities.push_back(newCity);
    }

    return all_cities;
}

int execSlave(int rank, int nthreads, vector<city> cities) {

    /*
    if(rank > GRID_LENGTH*GRID_LENGTH) {
        cout << "unused rank: " << rank << endl;
        // we don't need this process
        // we already have enough to cover whole grid
        return -1;
    }
    */

    //cout << "about to receive " << rank << endl;


    if(cities.size() == 0) {
        // get cities for this block
        int numCities;
        MPI_Recv(&numCities, 1, MPI_INT, 0, 0, cartcomm, MPI_STATUS_IGNORE);
        float x_coordinates[numCities];
        float y_coordinates[numCities];
        int ids[numCities];
        MPI_Recv(x_coordinates, numCities, MPI_FLOAT, 0, 1, cartcomm, MPI_STATUS_IGNORE);
        MPI_Recv(y_coordinates, numCities, MPI_FLOAT, 0, 2, cartcomm, MPI_STATUS_IGNORE);
        MPI_Recv(ids, numCities, MPI_INT, 0, 3, cartcomm, MPI_STATUS_IGNORE);
        /* Init cities vector for this one block */
        for (int i = 0; i < numCities; i++) {
            city newCity;
            newCity.x = x_coordinates[i];
            newCity.y = y_coordinates[i];
            newCity.id = ids[i];
            cities.push_back(newCity);
        }
    }

    //cout << "received cities" << endl;



    // Run the sequential tsp for this block
    thread_vars *vars = new thread_vars();
    vars->cities = cities;
    solution sol = startDynamicSolution(vars, vars->cities);

    int row = rank / dims[1] + 1;
    int col = rank % dims[1] + 1;
    //cout << "(" << row << "," << col << ") " << sol.distance << endl;

    /*
    // send back the results
    MPI_Send(&sol.distance, 1, MPI_FLOAT, 0, 4, MPI_COMM_WORLD);
    MPI_Send(&sol.first_city, 1, MPI_INT, 0, 5, MPI_COMM_WORLD);
    MPI_Send(&sol.last_city, 1, MPI_INT, 0, 6, MPI_COMM_WORLD);
    */

    // TODO
    // increment i to figure out iteration
    // ONLY SEND if all the ones to the left have sent
    // every process should send only once
    // except for last process, which just receives (some or all of iterations)
    // if rank % 2^(i+1) = 0, this process must receive
    // if haven't sent, and not receiving, try to send
    // receive from any

    cout << rank << ", (" << row << ", " << col << ")" << endl;
    MPI_Barrier(cartcomm);

    vector<city> cities_ordered_by_path;
    for(int i = 0; i < sol.path.size(); i++) {
        cities_ordered_by_path.push_back(vars->cities[sol.path[i]]);
    }

    // TODO this can be random, since its possible a processor will receive from a different order of blocks ex 6 can receive from 5 first, or 4 then 5 if 5 is slow
    int num = 1;
    if(dims[1] == 1) {

    }
    else if(col == dims[1]) {
        // cout << "col: " << col << endl;

        int final_sender = 1;
        while(final_sender < dims[1]) {
            final_sender = final_sender<<1;
        }
        final_sender = final_sender>>1;

        if(row == 5) cout << "final sender: " << final_sender << endl;

        MPI_Status status;
        int n;

        unsigned int mask = 1;
        cout << "start: " << (col & (~mask)) << endl;
        while((col & (~mask)) != 0) {
            if(row == 5) cout << "src: " << (col & (~mask)) << endl;
            MPI_Recv(&n, 1, MPI_INT, (col & (~mask)) - 1 + (row-1)*dims[1], 4, cartcomm, &status);
            cout << "mask before: " << mask << endl;
            unsigned int temp = (col & (~mask));
            while((col & (~mask)) == temp) {
                mask = (mask << 1) + 1;
            }
            cout << "temp: " << temp << (col & (~mask)) << endl;

            cout << "mask after: " << mask << endl;

            if(row == 5) cout << "received from " << status.MPI_SOURCE % dims[1] + 1 << endl;
            num+=n;
        }
        /*while(status.MPI_SOURCE % dims[1] + 1 != final_sender) {
            MPI_Recv(&n, 1, MPI_INT, MPI_ANY_SOURCE, 4, cartcomm, &status);
            if(row == 5) cout << "received from " << status.MPI_SOURCE % dims[1] + 1 << endl;
            num+=n;
        }*/
        cout << "row " << row << " NUM " << num << endl;
    }
    else {
        for(int i = 0; (1<<i) <= dims[1]; i++) {

            /*
            if(col == dims[1]) {
            if(row == 5) cout << "waiting to receive " << col << endl;

            int n;
            MPI_Recv(&n, 1, MPI_INT, MPI_ANY_SOURCE, 4, cartcomm, MPI_STATUS_IGNORE);
            num += n;
            if(row == 5) cout << "received! " << col << endl;

            if(row == 5) cout << "NUM " << num << endl;

            break;
            }*/
            if(col % (1<<(i+1)) == 0) {

                if(row == 5) cout << "waiting to receive " << col << endl;
                int n;
                MPI_Recv(&n, 1, MPI_INT, MPI_ANY_SOURCE, 4, cartcomm, MPI_STATUS_IGNORE);
                num += n;
                if(row == 5) cout << "received! " << col << endl;

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
                int dest = col + (1<<i) - 1 >= dims[1]-1 ? dims[1]-1 + (row-1)*dims[1] : (col - 1 + (1<<i)) + (row-1)*dims[1];
                if(row == 5) cout << "sending to " << dest+1 << " " << col << endl;
                MPI_Ssend(&num, 1, MPI_INT, dest, 4, cartcomm);
                break;
            }

        }
    }

    /*
    if(rank == 0) {
        for(int i = 0; i < cities.size(); i++) {
            cout << cities[i].x << " " << cities[i].y << endl;
        }
        cout <<  "ready?" << endl << std::flush;
        cout << "SIZE " << sol.path.size() << endl;
        for(int i = 0; i < sol.path.size(); i++) {
            cout << sol.path[i] << endl;
        }
    }*/
    return 0;
}



int execMain(int rank, int nthreads, int argc, char* argv[]) {

    // Log level
    if(argc > 1 && strcmp(argv[1], "v")==0) logLevel = 1; // v for verbose
    else logLevel = 0;


    cout << "number of processors: " << nthreads << endl;
    cout << "dims: " << dims[0] << ", " << dims[1] << endl;


    // Take input
    vector< vector <vector<city> > > blocks = generate_cities(dims[0], dims[1], 8);


    // Timing
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);


    /**
     * We now have cities generated for each block, already with ids
     * We have a grid shape of dim[0] by dim[1], which should be square if the number of processors used is a square
     */

    /*
    // Seperate cities into rows
    sort(all_cities.begin(), all_cities.end(), compareCitiesByY);
    vector< vector<city> > rows;
    // ensure that we get each city exactly once
    int chunk_size = floor(1.0 * all_cities.size() / GRID_LENGTH);
    int extra = all_cities.size() % GRID_LENGTH;
    for(int i = 0; i < GRID_LENGTH; i++) {
        vector<city>::iterator start = all_cities.begin() + chunk_size * i + min(i, extra);
        vector<city>::iterator end = all_cities.begin() + chunk_size * (i+1) + min(i+1, extra);
        vector<city> newVector(start, end);
        rows.push_back(newVector);
    }

    // Seperate rows into columns
    for(int i = 0; i < rows.size(); i++) {
        sort(rows[i].begin(), rows[i].end(), compareCitiesByX);
    }
    vector< vector< vector<city> > > blocks(GRID_LENGTH);
    for(int i = 0; i < GRID_LENGTH; i++) {
        // ensure that we get each city exactly once
        chunk_size = floor(1.0 * rows[i].size() / GRID_LENGTH);
        extra = rows[i].size() % GRID_LENGTH;
        for(int j = 0; j < GRID_LENGTH; j++) {
            vector<city>::iterator start = rows[i].begin() + chunk_size * j + min(j, extra);
            vector<city>::iterator end = rows[i].begin() + chunk_size * (j+1) + min(j+1, extra);
            vector<city> newVector(start, end);
            blocks[i].push_back(newVector);
        }
    }
     */


    // print how many cities are in each block
    if(logLevel == 1) {
        cout << endl << "city distribution:" << endl;
        for(int i = 0; i < blocks.size(); i++) {
            for(int j = 0; j < blocks[i].size(); j++) {
                if(logLevel==1) cout << blocks[i][j].size() << " ";
            }
            if(logLevel==1) cout << endl;
        }
    }


    // vector<city> cities_by_id;
    // add cities to vector, and set ids based on order in the grid
    // vector will come out sorted
    int count = 0;
    for(int i = 0; i < blocks.size(); i++) {

        // loop through each column, unless we are in the first row. Postpone the first block for last
        for(int j = (i==0 ? 1 : 0); j < blocks[i].size(); j++) {

            /*
            // Set city ids and add to vector
            for(int k = 0; k < blocks[i][j].size(); k++) {
                blocks[i][j][k].id = count;
                count++;
                cities_by_id.push_back(blocks[i][j][k]);
            }
            */

            // Prepare to send data to the processor that will handle this block

            int numCities = blocks[i][j].size();
            float x_coordinates[numCities];
            float y_coordinates[numCities];
            int ids[numCities];
            for(int k = 0; k < numCities; k++) {
                x_coordinates[k] = blocks[i][j][k].x; 
                y_coordinates[k] = blocks[i][j][k].y; 
                ids[k] = blocks[i][j][k].id;
            }
            int dest = i*dims[0]+j;

            // Send data to processor that will handle this block
            MPI_Send(&numCities, 1, MPI_INT, dest, 0, cartcomm);
            MPI_Send(x_coordinates, numCities, MPI_FLOAT, dest, 1, cartcomm);
            MPI_Send(y_coordinates, numCities, MPI_FLOAT, dest, 2, cartcomm);
            MPI_Send(ids, numCities, MPI_INT, dest, 3, cartcomm);


        }
    }

    cout << "done sending cities" << endl;

    /*

    // TODO
    // 1) Make the tsp solution return the calculated path
    // 2) Remove first_city, last_city
    // 3) Implement some nice parallel stitching using new algorithm
    // 4) Compile the paths together into a single vector or structure of some sort
    
    // store the "solutions" here, each containing the distance, first_city, and last_city in the path
    solution **solutionArray = new solution*[blocks.size()];
    for(int i = 0; i < GRID_LENGTH; i++) {
        solutionArray[i] = new solution[blocks[i].size()];
        for(int j = 0; j < GRID_LENGTH; j++) {
            int src = i*GRID_LENGTH+j+1;
            int numCities = blocks[i][j].size();

            float distance;
            int first_city;
            int last_city;

            MPI_Recv(&distance, 1, MPI_FLOAT, src, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&first_city, 1, MPI_INT, src, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&last_city, 1, MPI_INT, src, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 

            solution sol;
            sol.distance = distance;
            sol.first_city = first_city;
            sol.last_city = last_city;

            solutionArray[i][j] = sol;
        }
    }

     */
    /*
    if(logLevel == 1) {
        cout << endl << "Cities: " << endl;
        for(int i = 0; i < cities_by_id.size(); i++) {
            cout << setw(10) << left << cities_by_id[i].id << setw(8) << cities_by_id[i].x << setw(8) << cities_by_id[i].y << endl;
        }
    }
     */


    /*
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
    city fromCity = blocksAtEnd[0].visited_cities==FIRST ? cities_by_id[blocksAtEnd[0].last_city] : cities_by_id[blocksAtEnd[0].first_city];
    city toCity = blocksAtEnd[1].visited_cities==FIRST ? cities_by_id[blocksAtEnd[0].last_city] : cities_by_id[blocksAtEnd[0].first_city];
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
     */

    execSlave(rank, nthreads, blocks[0][0]);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    uint64_t delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    printf("took %lu\n", delta_us);


    return 0;

}

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int periods[2] = {0,0};

    int rank, nthreads;
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);

    MPI_Dims_create(nthreads, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims,
                    periods, 1, &cartcomm);

    MPI_Comm_rank(cartcomm, &rank);
    MPI_Comm_size(cartcomm, &nthreads);

    // Calulate Grid Length
    GRID_LENGTH = floor(sqrt(nthreads));


    vector<city> cities;
    if(rank==0) {
        execMain(rank, nthreads, argc, argv);
    } else {
        execSlave(rank, nthreads, cities);
    }

    MPI_Finalize();

    return 0;
}






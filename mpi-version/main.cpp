
#include "tsp.cpp"

/**
 * Gives an approximate solution to the TSP. Seperates cities into "blocks" on a
 * grid. Uses pthread solve a partial TSP for each block. Then connects the blocks
 * using a stitching algorithm.
 */

using namespace std;

int GRID_LENGTH;
int logLevel;

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


int execMain(int rank, int nthreads, int argc, char* argv[]) {

    // Log level
    if(argc > 1 && strcmp(argv[1], "v")==0) logLevel = 1; // v for verbose
    else logLevel = 0;

    cout << "number of processors: " << nthreads << endl;
    if(logLevel==1) cout << "Grid Length: " << GRID_LENGTH << endl;

    // Timing
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    // Take input
    vector<city> all_cities = readInCities();



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

    vector<city> cities_by_id;
    // add cities to vector, and set ids based on order in the grid
    // vector will come out sorted
    int count = 0;
    for(int i = 0; i < blocks.size(); i++) {
        for(int j = 0; j < blocks[i].size(); j++) {

            // Set city ids and add to vector
            for(int k = 0; k < blocks[i][j].size(); k++) {
                blocks[i][j][k].id = count;
                count++;
                cities_by_id.push_back(blocks[i][j][k]);
            }

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
            int dest = i*GRID_LENGTH+j+1;

            // Send data to processor that will handle this block
            MPI_Send(&numCities, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
            MPI_Send(x_coordinates, numCities, MPI_FLOAT, dest, 1, MPI_COMM_WORLD);
            MPI_Send(y_coordinates, numCities, MPI_FLOAT, dest, 2, MPI_COMM_WORLD);
            MPI_Send(ids, numCities, MPI_INT, dest, 3, MPI_COMM_WORLD);
            
        }
    }


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



    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    uint64_t delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    printf("took %lu\n", delta_us);


    return 0;

}

int execSlave(int rank) {

    if(rank > GRID_LENGTH*GRID_LENGTH) {
        cout << "unused rank: " << rank << endl;
        // we don't need this process
        // we already have enough to cover whole grid
        return -1;
    }

    // need to get this from the main process
    int numCities;

    MPI_Recv(&numCities, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    float x_coordinates[numCities];
    float y_coordinates[numCities];
    int ids[numCities];

    MPI_Recv(x_coordinates, numCities, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(y_coordinates, numCities, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(ids, numCities, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


    vector<city> cities;
    for(int i = 0; i < numCities; i++) {
        city newCity;
        newCity.x = x_coordinates[i];
        newCity.y = y_coordinates[i];
        newCity.id = ids[i];
        cities.push_back(newCity);
    }
    thread_vars *vars = new thread_vars();
    vars->cities = cities;
    solution sol = startDynamicSolution(vars, vars->cities);

    for(int i = 0; i < sol.path.size(); i++) {
        cout << sol.path[i] << ",";
    }
    cout << endl;

    // send back the results
    MPI_Send(&sol.distance, 1, MPI_FLOAT, 0, 4, MPI_COMM_WORLD);
    MPI_Send(&sol.first_city, 1, MPI_INT, 0, 5, MPI_COMM_WORLD);
    MPI_Send(&sol.last_city, 1, MPI_INT, 0, 6, MPI_COMM_WORLD);

    return 0;
}


int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int rank, nthreads;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);

    // Calulate Grid Length
    // Necessary to decide what each thread should do, so do it before
    GRID_LENGTH = floor(sqrt(nthreads-1));

    if(rank==0) {
        execMain(rank, nthreads, argc, argv);
    } else {
        execSlave(rank);
    }

    MPI_Finalize();

    return 0;
}






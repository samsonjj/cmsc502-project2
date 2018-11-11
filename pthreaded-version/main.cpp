
#include "tsp.cpp"

/**
 * Gives an approximate solution to the TSP. Seperates cities into "blocks" on a
 * grid. Uses pthread solve a partial TSP for each block. Then connects the blocks
 * using a stitching algorithm.
 */

using namespace std;

const int MAX_CITIES_PER_BLOCK = 10;
int GRID_LENGTH;
int logLevel;


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


bool compareCitiesByX(city a, city b) {
    return a.x < b.x;
}

bool compareCitiesByY(city a, city b) {
    return a.y < b.y;
}


void *findAndStoreSolution(void* v) {
    thread_vars* vars = (thread_vars*) v;
    solution s = startDynamicSolution(vars, vars->cities);
    vars->solutionArray[vars->i][vars->j] = s;
}


int main(int argc, char* argv[]) {

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    if(argc > 1 && strcmp(argv[1], "v")==0) logLevel = 1; // v for verbose
    else logLevel = 0;


    vector<city> all_cities = readInCities();

    if(argc > 1 && atoi(argv[1]) > 0) GRID_LENGTH = atoi(argv[1]);
    else GRID_LENGTH = ceil(sqrt(1.0 * all_cities.size() / MAX_CITIES_PER_BLOCK));

    if(logLevel==1) cout << "Grid Length: " << GRID_LENGTH << endl;


    // seperate cities into rows
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

    // split each row up into "columns"
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

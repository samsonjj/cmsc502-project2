#include <stdio.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <limits>
#include <cstdint>


/**
 * This program will use a sequential dynamic TSP solution. Function is given S=unvisited cities, and s=source city, c=current city,
 * and must find the shortest path from c back to source through unvisited cities
 */


using namespace std;


struct city {
    float x;
    float y;
};

// algorithm variables
int source;
int current;
vector <city> cities;
vector<int> unvisited;

// distances array
float **distances_array;

float getDist(int city1, int city2) {
    return distances_array[city1][city2];
}

float dynamicSolution() {
    if (unvisited.size() == 0) {
        return getDist(current, source);
    } else {
        // explore all unvisited, ignore source
        float store_current = current;
        int nextCity;
        float min_dist = numeric_limits<float>::max();
        for (int i = 0; i < unvisited.size(); i++) {

            // set the new current
            nextCity = unvisited[i];
            current = unvisited[i];

            // set the city as visited, then check, then set as unvisited again afterwards
            unvisited.erase(unvisited.begin() + i);
            float dist = getDist(store_current, nextCity) + dynamicSolution();
            unvisited.insert(unvisited.begin() + i, nextCity);
            if (dist < min_dist) min_dist = dist;
        }
        current = store_current;
        return min_dist;
    }
}

float startDynamicSolution() {

    unvisited.clear();

    // start at 1 since 0 will be source, and will count as visited
    for (unsigned i = 1; i < cities.size(); i++) {
        unvisited.push_back(i);
    }

    // start solving
    source = 0;
    current = 0;
    return dynamicSolution();
}

int readInCities() {
    ifstream file;
    file.open("cities.data");
    if (!file.is_open()) {
        return 1;
    }

    string x_pos;
    string y_pos;
    while (file >> x_pos >> y_pos) {
        city newCity;
        newCity.x = atof(x_pos.c_str());
        newCity.y = atof(y_pos.c_str());
        cities.push_back(newCity);
    }
}

int computeDistancesArray() {
    distances_array = new float *[cities.size()];
    for (int i = 0; i < cities.size(); i++) {
        distances_array[i] = new float[cities.size()];
    }

    for (int i = 0; i < cities.size(); i++) {
        for (int j = 0; j < cities.size(); j++) {
            float distance = sqrt((cities[i].x - cities[j].x) * (cities[i].x - cities[j].x) +
                                  (cities[i].y - cities[j].y) * (cities[i].y - cities[j].y));
            distances_array[i][j] = distance;
        }
    }
}


int main() {

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);


    readInCities();

    cout << "number of cities: " << cities.size() << endl;

    computeDistancesArray();

    float dist = startDynamicSolution();

    cout << "minimum distance is " << dist << endl;


    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    uint64_t delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    printf("took %f seconds\n", delta_us / 1000000.0);


    return 0;
}

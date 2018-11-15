
#include "header.h"

using namespace std;

int computeDistanceArray(thread_vars *vars, std::vector <city> cities) {
    vars->distance_array = new float *[cities.size()];
    for (int i = 0; i < cities.size(); i++) {
        vars->distance_array[i] = new float[cities.size()];
        for (int j = 0; j < cities.size(); j++) {
            vars->distance_array[i][j] = calcDistance(cities[i].x, cities[i].y, cities[j].x, cities[j].y);
        }
    }
}

// I need to return the minimum path, as well as minimum distance
// just need to return the ids of the cities in order of the path

solution dynamicSolution(thread_vars *vars) {
    if (vars->unvisited.size() == 0) {
        vector<int> min_path;
        min_path.push_back(vars->currentPath.back());
        solution sol(vars->distance_array[vars->currentPath.back()][vars->currentPath[0]],
                     vars->cities[vars->currentPath.back()].id, min_path);
        return sol;
    } else {
        // explore all unvisited, ignore source
        float min_dist = std::numeric_limits<float>::max();
        int min_last_city;
        vector<int> min_path;
        for (int i = 0; i < vars->unvisited.size(); i++) {

            // set the new current
            vars->currentPath.push_back(vars->unvisited[i]);

            // set the city as visited, then check, then set as unvisited again afterwards
            vars->unvisited.erase(vars->unvisited.begin() + i);
            solution dSol = dynamicSolution(vars);
            vars->unvisited.insert(vars->unvisited.begin() + i, vars->currentPath.back());
            vars->currentPath.pop_back();

            float dist = vars->distance_array[vars->currentPath.back()][vars->unvisited[i]] + dSol.distance;
            if (dist < min_dist) {
                min_dist = dist;
                min_last_city = dSol.last_city;
                min_path = dSol.path;
            }
        }
        min_path.insert(min_path.begin(), vars->cities[vars->currentPath.back()].id);
        solution sol(min_dist, min_last_city, min_path);
        return sol;
    }
}

solution startDynamicSolution(thread_vars *vars, std::vector <city> cities) {

    vars->unvisited.clear();

    // start at 1 since 0 will be source, and will count as visited
    for (unsigned i = 1; i < cities.size(); i++) {
        vars->unvisited.push_back(i);
    }
    vars->currentPath.push_back(0);

    computeDistanceArray(vars, cities);

    solution theSolution = dynamicSolution(vars);
    theSolution.first_city = vars->cities[0].id;
    return theSolution;
}



vector<city> generate_cities(int b, int row, int col, int num_cities) {

    vector<city> cities;
    for (int k = 0; k < num_cities; k++) {
        // Generate cities which lay inside the ith by jth block (range of whole grid is 500x500)
        // Generate up to 2 decimal places
        float x = ((rand() % (50000 / b)) + (50000 / b) * (col-1)) / 100.0;
        float y = ((rand() % (50000 / b)) + (50000 / b) * (row-1)) / 100.0;

        x = floor(x * 100) / 100;
        y = floor(y * 100) / 100;

        city c(x, y, k);

        cities.push_back(c);
    }
    return cities;
}

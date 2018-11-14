
#include <stdio.h>
#include <mpi.h>

#include <vector>
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <pthread.h>

#include <bits/stdc++.h>
#include <sys/time.h>
#include <cstdint>

struct city {
    float x;
    float y;
    int id;

    city(void) {}
    city(int xi, int yi) : x(xi), y(yi) { }
    city(int xi, int yi, int idi) : x(xi), y(yi), id(idi) { }
};

struct solution {
    int first_city;
    int last_city;
    std::vector<int> path;
    float distance;
    bool visited;

    // used in stitching
    int visited_cities;
    
    solution(void) { }
    solution(float d, int last) : distance(d), last_city(last) { }
    solution(float d, int last, std::vector<int> p) : distance(d), last_city(last), path(p) { }
};

struct thread_vars {
    std::vector<int> currentPath;
    std::vector<int> unvisited;
    float** distance_array;
    std::vector<city> cities;

    // for writing
    // solution **solutionArray;
    // int i;
    // int j;
};

int NONE = 0;
int FIRST= 1;
int LAST = 2;
int BOTH = 3;

float calcDistance(float x1, float y1, float x2, float y2) {
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

float calcDistance(city city1, city city2) {
    return calcDistance(city1.x, city1.y, city2.x, city2.y);
}

void printMatrix(float** matrix, int r, int c) {

    for(int i = 0; i < r; i++) {
        for(int j = 0; j < c; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

bool compareCitiesByX(city a, city b) {
    return a.x < b.x;
}

bool compareCitiesByY(city a, city b) {
    return a.y < b.y;
}
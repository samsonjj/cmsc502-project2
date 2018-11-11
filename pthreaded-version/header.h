
#include <stdio.h>
//#include <mpi.h>

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
};

struct solution {
    int first_city;
    int last_city;
    float distance;
    bool visited;

    // used in stitching
    int visited_cities;
    
    solution(void) { }
    solution(float d, int last) : distance(d), last_city(last) { }
};

struct thread_vars {
    std::vector<int> currentPath;
    std::vector<int> unvisited;
    float** distance_array;
    std::vector<city> cities;

    // for writing
    solution **solutionArray;
    int i;
    int j;
};

int NONE = 0;
int FIRST= 1;
int LAST = 2;
int BOTH = 3;

float calcDistance(float x1, float y1, float x2, float y2) {
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

void printMatrix(float** matrix, int r, int c) {

    for(int i = 0; i < r; i++) {
        for(int j = 0; j < c; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}



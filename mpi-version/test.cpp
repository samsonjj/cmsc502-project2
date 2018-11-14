
#include "tsp.cpp"

using namespace std;

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

int main() {
//    vector<city> v = generate_cities(10, 7, 8, 10);
//
//    for(int i = 0; i < v.size(); i++) {
//        cout << "(" << v[i].x << "," << v[i].y << ") ";
//    }
//
//    thread_vars *vars = new thread_vars();
//    vars->cities = v;
//
//    solution sol = startDynamicSolution(vars, vars->cities);
//
////    cout << sol.distance << endl;
//
//    //for(int i = 0; i < sol.path.size(); i++) {
//    //    cout << sol.path[i] << ", ";
//    //}
//
//    cout << (-1 % 15) << endl;

    vector<city> cities = readInCities();

    float sum = 0;
    for(int i = 0; i < cities.size(); i++) {
        sum += calcDistance(cities[i], cities[(i+1)%cities.size()]);
        cout << calcDistance(cities[i], cities[(i+1)%cities.size()]) << endl;
        cout << cities[i].x << ", " << cities[i].y << " " << cities[(i+1)%cities.size()].x << ", " << cities[(i+1)%cities.size()].y << endl;
    }

    cout << "sum : " << sum << endl;
}
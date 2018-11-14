
#include "tsp.cpp"

using namespace std;

int main() {
    vector<city> v = generate_cities(10, 7, 8, 10);

    for(int i = 0; i < v.size(); i++) {
        cout << "(" << v[i].x << "," << v[i].y << ") ";
    }

    thread_vars *vars = new thread_vars();
    vars->cities = v;

    solution sol = startDynamicSolution(vars, vars->cities);

    cout << sol.distance << endl;

    for(int i = 0; i < sol.path.size(); i++) {
        cout << sol.path[i] << ", ";
    }

    cout << (-1 % 15) << endl;
}
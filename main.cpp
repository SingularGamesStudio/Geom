#include <iostream>

#include "geometry.h"

using namespace std;

int main() {
    int n;
    cin >> n;
    vector<Point> points(n);
    for (int i = 0; i < n; i++) {
        cin >> points[i].x >> points[i].y;
    }
    cout << (Polygon(points)).isConvex();
}
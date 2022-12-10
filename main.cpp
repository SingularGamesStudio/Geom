#include <iostream>

#include "geometry.h"

using namespace std;

int main() {
    Point a, b, c, d;
    cin >> a.x >> a.y >> b.x >> b.y >> c.x >> c.y >> d.x >> d.y;
    Point res = Line(a, b).intersect(Line(c, d));
    cout << res.x << " " << res.y;
}
#include <iostream>

const double eps = 1e-7;
bool equals(double x, double y) {
    return abs(x - y) < eps;
}

class Point {
   public:
    double x;
    double y;
    Point() : x(0), y(0) {}
    Point(double x, double y) : x(x), y(y) {}
};
bool operator==(const Point& x, const Point& y) {
    return equals(x.x, y.x) && equals(x.y, y.y);
}
bool operator!=(const Point& x, const Point& y) {
    return !(x == y);
}
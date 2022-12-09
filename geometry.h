#include <cmath>
#include <iostream>
#include <vector>

using std::vector, std::abs, std::sqrt;
const double eps = 1e-7;
const double pi = std::acos(-1);

bool equal(double x, double y) {
    return abs(x - y) < eps;
}

class Point {
   public:
    double x;
    double y;
    Point() : x(0), y(0) {}
    Point(double x, double y) : x(x), y(y) {}
    Point(const Point& another) : x(another.x), y(another.y) {}

    bool operator==(const Point& another) const {
        return equal(x, another.x) && equal(y, another.y);
    }

    bool operator!=(const Point& another) const {
        return !(*this == another);
    }

    Point& operator*=(double a) {
        x *= a;
        y *= a;
        return *this;
    }

    Point& operator+=(const Point& another) {
        x += another.x;
        y += another.y;
        return *this;
    }

    Point& operator-=(const Point& another) {
        x -= another.x;
        y -= another.y;
        return *this;
    }

    double operator*(const Point& another) const {
        return (x * another.x + y * another.y);
    }

    double operator%(const Point& another) const {
        return (x * another.y - y * another.x);
    }

    Point operator*(double a) const {
        Point res(*this);
        res *= a;
        return res;
    }
    Point operator+(const Point& another) const {
        Point res(*this);
        res += another;
        return res;
    }
    Point operator-(const Point& another) const {
        Point res(*this);
        res -= another;
        return res;
    }

    double abs() const {
        return sqrt(x * x + y * y);
    }

    void rotate(const Point& center, double angle) {
        *this -= center;
        double temp = x;
        x = x * cos(angle) - y * sin(angle);
        y = temp * sin(angle) + y * cos(angle);
        *this += center;
    }

    void reflect(const Line& axis);

    void scale(const Point& center, double coefficient) {
        Point dir = *this - center;
        *this = center + dir * coefficient;
    }

    void reflect(const Point& center) {
        scale(center, -1);
    }
};

class Line {};

class Shape {
   public:
    virtual double perimeter() const = 0;

    virtual double area() const = 0;

    virtual bool equals(const Shape&) const {
        return false;
    }

    bool operator==(const Shape& another) {
        return equals(another);
    }

    virtual bool isCongruentTo(const Shape&) const {
        return false;
    }

    virtual bool isSimilarTo(const Shape&) const {
        return false;
    }

    virtual bool containsPoint(const Point&) const = 0;

    virtual void rotate(const Point&, double) = 0;

    virtual void reflect(const Line&) = 0;

    virtual void scale(const Point&, double) = 0;

    void reflect(const Point& center) {
        scale(center, -1);
    }
};

class Polygon : Shape {
   private:
    Point getEdge(int i) const {
        return points[i] - points[(i - 1 + points.size()) % points.size()];
    }

    double getSin(int i) const {
        Point e1 = getEdge(i);
        Point e2 = getEdge((i + 1) % points.size());
        return (e1 % e2) / e1.abs() / e2.abs();
    }

    double getCos(int i) const {
        Point e1 = getEdge(i);
        Point e2 = getEdge((i + 1) % points.size());
        return (e1 * e2) / e1.abs() / e2.abs();
    }

   public:
    vector<Point> points;

    Polygon() {}

    Polygon(vector<Point> points) : points(points) {}

    size_t verticesCount() const {
        return points.size();
    }
    const vector<Point>& getVertices() const {
        return points;
    }

    bool isConvex() const {
        double angle0 = getSin(0);
        for (size_t i = 1; i < points.size(); i++) {
            if ((angle0 > 0) != (getSin(i) > 0)) {
                return 0;
            }
        }
        return 1;
    }

    double perimeter() const override {
        double ans = 0;
        for (size_t i = 0; i < points.size(); i++) {
            ans += getEdge(i).abs();
        }
        return ans;
    }

    double area() const override {
        Point center = Point(0, 0);
        double ans = (points[points.size() - 1] - center) % (points[0] - center);
        for (size_t i = 1; i < points.size(); i++) {
            ans += (points[i - 1] - center) % (points[i] - center);
        }
        return ans / 2.0;
    }

    bool equals(const Shape& another) const override {
        if (dynamic_cast<const Polygon*>(&another) == nullptr)
            return false;
        const Polygon& other = *(dynamic_cast<const Polygon*>(&another));
        if (other.points.size() != points.size())
            return false;
        for (size_t start = 0; start < points.size(); start++) {
            for (int i = 0; i < points.size(); i++) {
                if (points[i] != other.points[(start + i) % points.size()])
                    break;
                if (i == points.size() - 1)
                    return true;
            }
        }
        return false;
    }

    bool isCongruentTo(const Shape& another) const override {
        if (dynamic_cast<const Polygon*>(&another) == nullptr)
            return false;
        const Polygon& other = *(dynamic_cast<const Polygon*>(&another));
        if (other.points.size() != points.size())
            return false;
        for (size_t start = 0; start < points.size(); start++) {
            for (int i = 0; i < points.size(); i++) {
                if (getEdge(i).abs() != other.getEdge((start + i) % points.size()).abs() &&
                    getCos(i) != other.getCos((start + i) % points.size()) &&
                    getSin(i) != other.getSin((start + i) % points.size()))
                    break;
                if (i == points.size() - 1)
                    return true;
            }
            for (int i = points.size() - 1; i >= 0; i--) {
                if (getEdge(i).abs() != other.getEdge((start + i) % points.size()).abs() &&
                    getCos(i) != other.getCos((start + i) % points.size()) &&
                    getSin(i) != other.getSin((start + i) % points.size()))
                    break;
                if (i == 0)
                    return true;
            }
        }
        return false;
    }

    bool isSimilarTo(const Shape& another) const override {
        if (dynamic_cast<const Polygon*>(&another) == nullptr)
            return false;
        const Polygon& other = *(dynamic_cast<const Polygon*>(&another));
        if (other.points.size() != points.size())
            return false;
        for (size_t start = 0; start < points.size(); start++) {
            double coef = getEdge(0).abs() / other.getEdge((start + 0) % points.size()).abs();
            for (int i = 0; i < points.size(); i++) {
                if (getEdge(i).abs() != coef * other.getEdge((start + i) % points.size()).abs() &&
                    getCos(i) != other.getCos((start + i) % points.size()) &&
                    getSin(i) != other.getSin((start + i) % points.size()))
                    break;
                if (i == points.size() - 1)
                    return true;
            }
            coef = getEdge(points.size() - 1).abs() / other.getEdge((start + points.size() - 1) % points.size()).abs();
            for (int i = points.size() - 1; i >= 0; i--) {
                if (getEdge(i).abs() != coef * other.getEdge((start + i) % points.size()).abs() &&
                    getCos(i) != other.getCos((start + i) % points.size()) &&
                    getSin(i) != other.getSin((start + i) % points.size()))
                    break;
                if (i == 0)
                    return true;
            }
        }
        return false;
    }

    void rotate(const Point& center, double angle) override {
        for (auto& p : points) {
            p.rotate(center, angle);
        }
    }

    void reflect(const Line& axis) override {
        for (auto& p : points) {
            p.reflect(axis);
        }
    }

    void scale(const Point& center, double coefficient) override {
        for (auto& p : points) {
            p.scale(center, coefficient);
        }
    }
};

class Ellipse : Shape {
   private:
    Point f1, f2;
    double d;

    double a() const {
        return d / 2.0;
    }

    double c() const {
        return (f1 - f2).abs() / 2.0;
    }

    double b() const {
        double a1 = a();
        double c1 = c();
        return sqrt(a1 * a1 - c1 * c1);
    }

   public:
    Ellipse() {}

    Ellipse(Point f1, Point f2, double d) : f1(f1), f2(f2), d(d) {}

    std::pair<Point, Point> focuses() const {
        return {f1, f2};
    }

    Point center() const {
        return (f1 + f2) * 0.5;
    }

    double eccentricity() const {
        return c() / a();
    }

    double perimeter() const override {
        return pi * (a() + b());
    }

    double area() const override {
        return pi * a() * b();
    }

    bool equals(const Shape& another) const override {
        if (dynamic_cast<const Ellipse*>(&another) == nullptr)
            return false;
        const Ellipse& other = *(dynamic_cast<const Ellipse*>(&another));
        return equal(d, other.d) && ((f1 == other.f1 && f2 == other.f2) || (f2 == other.f1 && f1 == other.f2));
    }

    virtual bool isCongruentTo(const Shape& another) const override {
        if (dynamic_cast<const Ellipse*>(&another) == nullptr)
            return false;
        const Ellipse& other = *(dynamic_cast<const Ellipse*>(&another));
        return equal(a(), other.a()) && equal(b(), other.b());
    }

    virtual bool isSimilarTo(const Shape&) const override {
        if (dynamic_cast<const Ellipse*>(&another) == nullptr)
            return false;
        const Ellipse& other = *(dynamic_cast<const Ellipse*>(&another));
    }

    virtual bool containsPoint(const Point&) const override {
        return false;
    }

    virtual void rotate(const Point&, double) override = 0;

    virtual void reflect(const Line&) override = 0;

    virtual void scale(const Point&, double) override = 0;
};

/*
Класс Line - прямая. Прямую можно задать двумя точками, можно двумя числами (угловой коэффициент и сдвиг), можно точкой и числом (угловой коэффициент). Линии можно сравнивать операторами == и !=.
Можно сконструировать многоугольник из точек, передаваемых в качестве параметров через запятую (т.е. неуказанное число аргументов).
Ellipse::std::pair<Line, Line> directrices() - пару его директрис;
Point::reflect(line)

Класс Circle - круг. Круг - частный случай эллипса. У круга можно спросить double radius() - радиус. Круг можно задать точкой и числом (центр и радиус).
Класс Rectangle - прямоугольник. Прямоугольник - частный случай многоугольника. У прямоугольника можно спросить Point center() - его центр; std::pair<Line, Line> diagonals() - пару его диагоналей. Прямоугольник можно сконструировать по двум точкам (его противоположным вершинам) и числу (отношению смежных сторон), причем из двух таких прямоугольников выбирается тот, у которого более короткая сторона расположена по левую сторону от диагонали, если смотреть от первой заданной точки в направлении второй.
Класс Square - квадрат. Квадрат - частный случай прямоугольника. У квадрата можно спросить Circle circumscribedCircle(), Circle inscribedCircle(). Квадрат можно задать двумя точками - противоположными вершинами.
Класс Triangle - треугольник. Треугольник - частный случай многоугольника. У треугольника можно спросить Circle circumscribedCircle(), Circle inscribedCircle(), Point centroid() - его центр масс, Point orthocenter() - его ортоцентр, Line EulerLine() - его прямую Эйлера, Circle ninePointsCircle() - его окружность Эйлера.
У любой фигуры можно спросить:

double perimeter() - периметр;
double area() - площадь;
bool operator==(const Shape& another) - совпадает ли эта фигура с другой как множество точек. (В частности, треугольник ABC равен треугольнику BCA.)
bool isCongruentTo(const Shape& another) - равна ли эта фигура другой в геометрическом смысле, то есть можно ли совместить эти фигуры движением плоскости. Движение – это отображение плоскости на себя, сохраняющее расстояния.
bool isSimilarTo(const Shape& another) - подобна ли эта фигура другой, то есть можно ли перевести одну фигуру в другую преобразованием подобия. (Определение преобразования подобия, кто не знает, можно посмотреть в Википедии.)
bool containsPoint(const Point& point) - находится ли точка внутри фигуры.

rotate(const Point& center, double angle) - поворот на угол (в градусах, против часовой стрелки) относительно точки;
reflect(const Line& axis) - симметрию относительно прямой;
scale(const Point& center, double coefficient) - гомотетию с коэффициентом coefficient и центром center.
*/
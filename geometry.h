#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using std::vector, std::abs, std::sqrt;
const double eps = 1e-5;
const double pi = std::acos(-1);
const double deg2rad = 2.0 * pi / 360.0;

bool equal(double x, double y) {
    return abs(x - y) < eps;
}

bool lessorequal(double x, double y) {
    return x < (y + eps);
}

int sign(double d) {
    if (equal(d, 0))
        return 0;
    return (d < 0 ? -1 : 1);
}
class Point {
   public:
    double x;
    double y;
    Point() : x(0), y(0) {}
    Point(double x, double y) : x(x), y(y) {}
    Point(const Point& another) : x(another.x), y(another.y) {}

    Point& operator=(const Point& anoter) {
        x = anoter.x;
        y = anoter.y;
        return *this;
    }

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

    Point& operator/=(double a) {
        x /= a;
        y /= a;
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
    Point operator/(double a) const {
        Point res(*this);
        res /= a;
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

    Point perp() const {
        return Point(-y, x);
    }

    Point normalized() const {
        return Point(x / abs(), y / abs());
    }

    bool within(Point a, Point b) const {
        if (lessorequal(std::min(a.x, b.x), x) && lessorequal(x, std::max(a.x, b.x)))
            if (lessorequal(std::min(a.y, b.y), y) && lessorequal(y, std::max(a.y, b.y)))
                return true;
        return false;
    }

    double abs() const {
        return sqrt(x * x + y * y);
    }

    void rotate(const Point& center, double angle) {
        angle *= deg2rad;
        *this -= center;
        double temp = x;
        x = x * cos(angle) - y * sin(angle);
        y = temp * sin(angle) + y * cos(angle);
        *this += center;
    }

    void scale(const Point& center, double coefficient) {
        Point dir = *this - center;
        *this = center + dir * coefficient;
    }

    void reflect(const Point& center) {
        scale(center, -1);
    }
};

class Line {
   protected:
    Point p1;
    Point p2;

   public:
    Line() {}
    Line(Point p1, Point p2) : p1(p1), p2(p2) {}
    Line(Point p1, double k) : p1(p1), p2(Point(p1.x + 1, p1.y + k)) {}
    Line(double k, double b) : p1(Point(0, b)), p2(Point(1, k + b)) {}

    bool parralel(const Line& another) const {
        return equal((p1 - p2).x * (another.p1 - another.p2).y, (p1 - p2).y * (another.p1 - another.p2).x);
    }

    Point project(const Point& point) const {
        Point diff = (p2 - p1);
        return p1 + diff * ((point - p1) * diff) / diff.abs() / diff.abs();
    }

    Point intersect(const Line& another) const {
        Point pr1 = project(another.p1);
        Point pr2 = project(another.p2);
        if (pr2 == pr1)
            return pr1;
        double a = -(another.p1 - pr1).abs() * sign((another.p1 - pr1) % (pr2 - pr1));
        double b = (another.p2 - pr2).abs() * sign((another.p2 - pr2) % (pr2 - pr1));
        return pr1 + (pr2 - pr1) * (a / (a + b));
    }

    bool operator==(const Line& another) const {
        Line intersection;
        if (p1 != another.p1)
            intersection = Line(p1, another.p1);
        else
            intersection = Line(p1, another.p2);
        return parralel(intersection) && parralel(another);
    }

    bool operator!=(const Line& another) const {
        return !(*this == another);
    }
};

void refl(Point& point, const Line& axis) {
    Point projection = axis.project(point);
    point = point + (projection - point) * 2;
}

class Shape {
   public:
    virtual double perimeter() const = 0;

    virtual double area() const = 0;

    virtual bool equals(const Shape&) const {
        return false;
    }

    bool operator==(const Shape& another) const {
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
    virtual ~Shape() {}
};

class Polygon : public Shape {
   protected:
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

    template <typename T>
    Polygon(T i) {
        static_assert(std::is_same<T, Point>::value);
        points.push_back(i);
    }

    template <typename T, typename... R>
    Polygon(T i, R... r) : Polygon(r...) {
        static_assert(std::is_same<T, Point>::value);
        points.push_back(i);
    }

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
            for (size_t i = 0; i < points.size(); i++) {
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
            for (size_t i = 0; i < points.size(); i++) {
                if (!equal(getEdge(i).abs(), other.getEdge((start + i) % points.size()).abs()) ||
                    !equal(getCos(i), other.getCos((start + i) % points.size())) ||
                    !equal(getSin(i), other.getSin((start + i) % points.size())))
                    break;
                if (i == points.size() - 1)
                    return true;
            }
            for (int i = 0; i < static_cast<int>(points.size()); i++) {
                if (!equal(getEdge(i).abs(), other.getEdge((start - i + 1 + points.size()) % points.size()).abs()) ||
                    !equal(getCos(i), other.getCos((start - i + points.size()) % points.size())) ||
                    !equal(getSin(i), other.getSin((start - i + points.size()) % points.size())))
                    break;
                if (i == static_cast<int>(points.size()) - 1)
                    return true;
            }
        }
        return false;
    }

    bool isSimilarTo(const Shape& another) const override {
        if (dynamic_cast<const Polygon*>(&another) == nullptr)
            return false;
        const Polygon& other = *(dynamic_cast<const Polygon*>(&another));
        if (other.points.size() != points.size()) return false;
        for (size_t start = 0; start < points.size(); start++) {
            double coef = getEdge(0).abs() / other.getEdge(start).abs();
            for (size_t i = 0; i < points.size(); i++) {
                if (!equal(getEdge(i).abs(), coef * other.getEdge((start + i) % points.size()).abs()) ||
                    !equal(getCos(i), other.getCos((start + i) % points.size())) ||
                    !equal(getSin(i), other.getSin((start + i) % points.size()))) {
                    break;
                }
                if (i == points.size() - 1)
                    return true;
            }
            coef = getEdge(0).abs() / other.getEdge((start + 1 + points.size()) % points.size()).abs();
            for (int i = 0; i < static_cast<int>(points.size()); i++) {
                if (!equal(getEdge(i).abs(), coef * other.getEdge((start - i + 1 + points.size()) % points.size()).abs()) ||
                    !equal(getCos(i), other.getCos((start - i + points.size()) % points.size())) ||
                    !equal(getSin(i), -other.getSin((start - i + points.size()) % points.size()))) {
                    break;
                }
                if (i == static_cast<int>(points.size()) - 1)
                    return true;
            }
        }
        return false;
    }

    bool containsPoint(const Point& point) const override {
        for (size_t i = 0; i < points.size(); i++) {
            if (points[i] == point)
                return true;
        }
        for (size_t i = 0; i < points.size(); i++) {
            Line edge = Line(points[i], points[i] - getEdge(i));
            if (equal(0, (edge.project(point) - point).abs()) && point.within(points[i], points[i] - getEdge(i))) {
                return true;
            }
        }
        std::mt19937 rnd(42);
        Line raycast = Line(point, Point(point.x + static_cast<double>(rnd()) / 10000000.0 + 1, point.y + static_cast<double>(rnd()) / 10000000.0 + 1));
        bool inside = false;
        for (size_t i = 0; i < points.size(); i++) {
            Line edge = Line(points[i], points[i] - getEdge(i));
            if (!edge.parralel(raycast)) {
                Point intersection = edge.intersect(raycast);
                if (intersection.x > point.x && intersection.within(points[i], points[i] - getEdge(i)))
                    inside = !inside;
            }
        }
        return inside;
    }

    void rotate(const Point& center, double angle) override {
        for (auto& p : points) {
            p.rotate(center, angle);
        }
    }

    void reflect(const Line& axis) override {
        for (auto& p : points) {
            refl(p, axis);
        }
    }

    void scale(const Point& center, double coefficient) override {
        for (auto& p : points) {
            p.scale(center, coefficient);
        }
    }
};

class Ellipse : public Shape {
   protected:
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

    std::pair<Line, Line> directrices() {
        Point cent = center();
        Point delta = (f1 - f2).normalized() * a() * a() / c();
        return {Line(cent + delta, cent + delta + delta.perp()), Line(cent - delta, cent - delta + delta.perp())};
    }

    double eccentricity() const {
        return c() / a();
    }

    double perimeter() const override {
        // std::cerr << f1.x << " " << f1.y << " " << f2.x << " " << f2.y << " " << d << " " << a() << " " << b() << " " << 2.0 * pi * sqrt((a() * a() + b() * b()) / 2.0) << "\n";
        return pi * (3 * (a() + b()) - sqrt((3 * a() + b()) * (a() + 3 * b())));
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

    virtual bool isSimilarTo(const Shape& another) const override {
        if (dynamic_cast<const Ellipse*>(&another) == nullptr)
            return false;
        const Ellipse& other = *(dynamic_cast<const Ellipse*>(&another));
        return equal(eccentricity(), other.eccentricity());
    }

    virtual bool containsPoint(const Point& point) const override {
        return lessorequal((point - f1).abs() + (point - f2).abs(), d);
    }

    virtual void rotate(const Point& center, double angle) override {
        f1.rotate(center, angle);
        f2.rotate(center, angle);
    }

    virtual void reflect(const Line& axis) override {
        refl(f1, axis);
        refl(f2, axis);
    }

    virtual void scale(const Point& center, double coefficient) override {
        f1.scale(center, coefficient);
        f2.scale(center, coefficient);
        d *= abs(coefficient);
    }
};

class Circle : public Ellipse {
   public:
    Circle() {}
    Circle(Point center, double r) : Ellipse(center, center, r * 2) {}

    std::pair<Line, Line> directrices() {
        assert(0);
        return {Line(), Line()};
    }

    double radius() {
        return d / 2;
    }
};

class Rectangle : public Polygon {
   public:
    Rectangle() {}
    Rectangle(Point p1, Point p2, double d) {
        if (d < 1)
            d = 1 / d;
        double angle = ((pi - std::atan(d)) * 2) / deg2rad;
        Point mid = (p1 + p2) / 2;
        points.push_back(p1);
        p1.rotate(mid, -angle);
        points.push_back(p1);
        points.push_back(p2);
        p1.reflect(mid);
        points.push_back(p1);
    }
    Point center() const {
        return (points[0] + points[2]) / 2;
    }
    std::pair<Line, Line> diagonals() const {
        return {Line(points[0], points[2]), Line(points[1], points[3])};
    }
};

class Square : public Rectangle {
   public:
    Square() {}
    Square(const Point& p1, const Point& p2) : Rectangle(p1, p2, 1) {}

    Circle circumscribedCircle() const {
        return Circle(center(), (points[0] - center()).abs());
    }

    Circle inscribedCircle() const {
        return Circle(center(), ((points[0] + points[1]) / 2 - center()).abs());
    }
};

class Triangle : public Polygon {
   private:
    Point serpercenter() const {
        Line serper1 = Line(points[0] + (points[1] - points[0]) / 2, points[0] + (points[1] - points[0]) / 2 + (points[1] - points[0]).perp());
        Line serper2 = Line(points[0] + (points[2] - points[0]) / 2, points[0] + (points[2] - points[0]) / 2 + (points[2] - points[0]).perp());
        return serper1.intersect(serper2);
    }

    Point bissectralcenter() const {
        Line bissectrice1 = Line(points[0], points[0] + (points[1] - points[0]).normalized() + (points[2] - points[0]).normalized());
        Line bissectrice2 = Line(points[1], points[1] + (points[2] - points[1]).normalized() + (points[0] - points[1]).normalized());
        return bissectrice1.intersect(bissectrice2);
    }

   public:
    using Polygon::Polygon;
    Triangle() {}

    Circle circumscribedCircle() const {
        Point center = serpercenter();
        return Circle(center, (center - points[0]).abs());
    }

    Circle inscribedCircle() const {
        Point center = bissectralcenter();
        return Circle(center, (Line(points[1], points[0]).project(center) - center).abs());
    }

    Point centroid() const {
        return (points[0] + points[1] + points[2]) / 3;
    }

    Point orthocenter() const {
        Line height1 = Line(points[0], Line(points[1], points[2]).project(points[0]));
        Line height2 = Line(points[1], Line(points[0], points[2]).project(points[1]));
        return height1.intersect(height2);
    }

    Line EulerLine() const {
        return Line(centroid(), orthocenter());
    }

    Circle ninePointsCircle() const {
        Point center = (orthocenter() + serpercenter()) / 2;
        return Circle(center, (center - (points[1] + points[0]) / 2).abs());
    }
};

/*
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
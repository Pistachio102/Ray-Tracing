//
// Created by nisa on 25/06/21.
//

#ifndef RAY_TRACING_1605102_CLASSES_H
#define RAY_TRACING_1605102_CLASSES_H


#define AMBIENT 0
#define DIFFUSE 1
#define SPECULAR 2
#define REFLECTION 3
extern int recursionLevel;

class Point {
public:
    double x, y, z;

    Point() = default;

    Point(double x1, double y1, double z1) {
        x = x1;
        y = y1;
        z = z1;
    }

    Point operator+(Point point) {
        return {this->x + point.x, this->y + point.y, this->z + point.z};
    }

    Point operator-(Point point) {
        return {this->x - point.x, this->y - point.y, this->z - point.z};
    }

    Point operator*(double mult) {
        return {this->x * mult, this->y * mult, this->z * mult};
    }

    Point operator/(double div) {
        return {this->x / div, this->y / div, this->z / div};
    }

    Point normalize() {
        return *this / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    }

    double length() const {
        return sqrt(x * x + y * y + z * z);
    }

    friend double dotProduct(Point point1, Point point2) {
        return point1.x * point2.x + point1.y * point2.y + point1.z * point2.z;
    }

    friend Point crossProduct(Point point1, Point point2) {
        double a = point1.y * point2.z - point1.z * point2.y;
        double b = point1.z * point2.x - point1.x * point2.z;
        double c = point1.x * point2.y - point1.y * point2.x;

        return {a, b, c};
    }

    friend Point afterRotation(Point v, Point r, double directionSign) {
        Point temp, normalized ;
        double angle = 0.05 * directionSign;

        normalized.x = r.y * v.z - r.z * v.y;
        normalized.y = r.z * v.x - r.x * v.z;
        normalized.z = r.x * v.y - r.y * v.x;

        temp.x = v.x * cos(angle) + normalized.x * sin(angle);
        temp.y = v.y * cos(angle) + normalized.y * sin(angle);
        temp.z = v.z * cos(angle) + normalized.z * sin(angle);

        return temp;
    }

    void print() const {
        cout << x << ", " << y << ", " << z << endl;
    }
};


class Ray {

public:
    Point start;
    Point dir;

    Ray(Point start, Point dir) {
        this->start = start;
        this->dir = dir.normalize();
    }
};


class Object {

public:
    char type = 'o';
    Point referencePoint;
    double height, width, length;
    double color[3];
    double coEfficients[4];
    double shine;
    double sourceFactor = 1.0, eta = 0.0;

    Object() {}

    virtual void draw() = 0;
    
    void setColor(double red, double green, double blue) {
        color[0] = red;
        color[1] = green;
        color[2] = blue;
    }
    void setShine(int shine1) {
        shine = shine1;
    }
    void setCoefficients(double a, double d, double s, double r) {
        coEfficients[AMBIENT] = a;
        coEfficients[DIFFUSE] = d;
        coEfficients[SPECULAR] = s;
        coEfficients[REFLECTION] = r;
    }
    
    virtual double getIntersection(Ray &ray) = 0;

    virtual Point getNormal(Point intersection) = 0;

    void intersect(Ray &ray, double t, double currentColor[3], int level);



    static Point reflectedRayDirection(Ray &ray, Point normal) {
        Point temp = ray.dir - normal * 2.0 * dotProduct(ray.dir, normal);
        return temp.normalize();
    }

    Point refractedRayDirection(Ray &ray, Point normal) const {
        double nDotI = dotProduct(normal, ray.dir);
        double k = 1.0 - eta * eta * (1.0 - nDotI * nDotI);

        if (k < 0.0)
            return {0.0, 0.0, 0.0};

        Point temp = ray.dir * eta - normal * (eta * nDotI + sqrt(k));
        return temp.normalize();
    }
};


class Sphere : public Object {
public:


    Sphere(Point center, double radius) {
        this->referencePoint = center;
        this->length = radius;
    }

    void draw() {
        glColor3f(color[0], color[1], color[2]);
        glPushMatrix();
        glTranslatef(referencePoint.x, referencePoint.y, referencePoint.z);
        glutSolidSphere(length, 100, 100);
        glPopMatrix();
    }

    double getIntersection(Ray &ray) {

        Point cen = referencePoint;
        Point ori = ray.start;
        Point l = ray.dir;

        Point d = ori - cen;
        double discriminant = dotProduct(l, d) * dotProduct(l, d) - dotProduct(d, d) + length * length;

        if (discriminant < 0)
            return -1;

        double res = sqrt(discriminant);
        double t1 = -dotProduct(l, d) + res;
        double t2 = -dotProduct(l, d) - res;

        return min(t1, t2);
    }

    Point getNormal(Point intersectingPoint) {
        Point normal = intersectingPoint - referencePoint;
        return normal.normalize();
    }
};


class Triangle : public Object {
public:
    double a, b, c, d;
    Point p1, p2, p3;


    Triangle(Point point1, Point point2, Point point3) {
        p1 = point1;
        p2 = point2;
        p3 = point3;

        Point AB = p2 - p1;
        Point AC = p3 - p1;
        Point temp = crossProduct(AB, AC);

        a = temp.x;
        b = temp.y;
        c = temp.z;
        d = a * p1.x + b * p1.y + c * p1.z;
    }

    void draw() {
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(p1.x, p1.y, p1.z);
            glVertex3f(p2.x, p2.y, p2.z);
            glVertex3f(p3.x, p3.y, p3.z);
        }
        glEnd();
    }

    Point getNormal(Point intersectingPoint) {
        Point normal = Point(a, b, c);
        return normal.normalize();
    }

    double getIntersection(Ray &ray) {
        Point normal = getNormal(Point(0, 0, 0));
        double t;
        bool boolean = false;
        double denominator = dotProduct(normal, ray.dir);

        if (denominator < 0.0) {
            normal = normal * -1.0;
            denominator = dotProduct(normal, ray.dir);
        }
        
        if (abs(denominator) < 0.0000001)
            return -1;

        t = (dotProduct(normal, p1 - ray.start)) / denominator;
        if (t >= 0)
            boolean = true;

        if (!boolean)
            return -1;

        Point intersectingPoint = ray.start + ray.dir * t;

        Point N = crossProduct(p2 - p1, p3 - p1);
        double area2 = N.length();

        Point area, side1, side2;

        side1 = p2 - p1;
        side2 = intersectingPoint - p1;
        area = crossProduct(side1, side2);
        if (dotProduct(N, area) < 0)
            return -1;

        side1 = p3 - p2;
        side2 = intersectingPoint - p2;
        area = crossProduct(side1, side2);
        if (dotProduct(N, area) < 0)
            return -1;

        side1 = p1 - p3;
        side2 = intersectingPoint - p3;
        area = crossProduct(side1, side2);
        if (dotProduct(N, area) < 0)
            return -1;

        return t; 
    }
};

class GeneralQuadratic : public Object {

public:
    double A1, A2, A3, A4, A5, A6, A7, A8, A9, A10;

    GeneralQuadratic(double coefficients[10], Point referencePoint, double width, double height, double length) {
        A1 = coefficients[0];
        A2 = coefficients[1];
        A3 = coefficients[2];
        A4 = coefficients[3];
        A5 = coefficients[4];
        A6 = coefficients[5];
        A7 = coefficients[6];
        A8 = coefficients[7];
        A9 = coefficients[8];
        A10 = coefficients[9];
        this->referencePoint = referencePoint;
        this->height = height;
        this->width = width;
        this->length = length;
    }

    void draw() {}

    Point getNormal(Point intersection) {

        double a = 2 * A1 * intersection.x + A4 * intersection.y + A6 * intersection.z + A7;
        double b = 2 * A2 * intersection.y + A4 * intersection.x + A5 * intersection.z + A8;
        double c = 2 * A3 * intersection.z + A5 * intersection.y + A6 * intersection.x + A9;

        Point normal(a, b, c);
        return normal.normalize();
    }

    double getIntersection(Ray &ray) {
        double a = A1 * ray.dir.x * ray.dir.x + A2 * ray.dir.y * ray.dir.y + A3 * ray.dir.z * ray.dir.z;
        double b = 2 * (A1 * ray.start.x * ray.dir.x + A2 * ray.start.y * ray.dir.y + A3 * ray.start.z * ray.dir.z);
        double c = A1 * ray.start.x * ray.start.x + A2 * ray.start.y * ray.start.y + A3 * ray.start.z * ray.start.z;

        a = a + (A4 * ray.dir.x * ray.dir.y + A5 * ray.dir.y * ray.dir.z + A6 * ray.dir.z * ray.dir.x);
        b = b + (A4 * (ray.start.x * ray.dir.y + ray.dir.x * ray.start.y) +
              A5 * (ray.start.y * ray.dir.z + ray.dir.y * ray.start.z) +
              A6 * (ray.start.z * ray.dir.x + ray.dir.z * ray.start.x));
        c = c + (A4 * ray.start.x * ray.start.y + A5 * ray.start.y * ray.start.z + A6 * ray.start.z * ray.start.x);
        b = b + (A7 * ray.dir.x + A8 * ray.dir.y + A9 * ray.dir.z);
        c = c + (A7 * ray.start.x + A8 * ray.start.y + A9 * ray.start.z + A10);

        double discriminant = b * b - 4 * a * c;

        if (discriminant < 0) {
            return -1;
        }

        double t1 = (-b + sqrt(discriminant)) / (2.0 * a);
        double t2 = (-b - sqrt(discriminant)) / (2.0 * a);
       // cout << t1 << " " << t2 << endl;

        Point intersectingPoint1 = ray.start + ray.dir * t1;
        Point intersectingPoint2 = ray.start + ray.dir * t2;

        double x_min = referencePoint.x;
        double x_max = x_min + length;

        double y_min = referencePoint.y;
        double y_max = y_min + width;

        double z_min = referencePoint.z;
        double z_max = z_min + height;

        double x1 = intersectingPoint1.x;
        double y1 = intersectingPoint1.y;
        double z1 = intersectingPoint1.z;

        double x2 = intersectingPoint2.x;
        double y2 = intersectingPoint2.y;
        double z2 = intersectingPoint2.z;

        bool bool1 = length > 0 && (x1 < x_min || x1 > x_max) ||
                     width > 0 && (y1 < y_min || y1 > y_max) ||
                     height > 0 && (z1 < z_min || z1 > z_max);

        bool bool2 = length > 0 && (x2 < x_min || x2 > x_max) ||
                     width > 0 && (y2 < y_min || y2 > y_max) ||
                     height > 0 && (z2 < z_min || z2 > z_max);

        if (bool1 && bool2)
            return -1;
        else if (bool1)
            return t2;
        else if (bool2)
            return t1;
        else
            return min(t1, t2);
    }
};

class Floor : public Object {

public:
    Point origin;
    int numberOfTiles;
    double floorWidth, tileWidth;


    Floor(double floorWidth1, double tileWidth1) {
        floorWidth = floorWidth1;
        tileWidth = tileWidth1;
        origin = Point(-floorWidth / 2.0, -floorWidth / 2.0, 0.0);
        numberOfTiles = floorWidth / tileWidth;

    }

    void draw() {
        glBegin(GL_QUADS);
        {
            for (int i = 0; i < numberOfTiles; i++)
                for (int j = 0; j < numberOfTiles; j++) {
                    bool c = (i + j) % 2;
                    glColor3f(c, c, c);

                    glVertex3f(origin.x + tileWidth * i, origin.y + tileWidth * j, origin.z);
                    glVertex3f(origin.x + tileWidth * (i + 1), origin.y + tileWidth * j, origin.z);
                    glVertex3f(origin.x + tileWidth * (i + 1), origin.y + tileWidth * (j + 1), origin.z);
                    glVertex3f(origin.x + tileWidth * i, origin.y + tileWidth * (j + 1), origin.z);
                }
        }
        glEnd();
    }

    Point getNormal(Point intersection) { return Point(0, 0, 1); }

    double getIntersection(Ray &ray) {
        if (ray.dir.z == 0)
            return -1;

        double t = -(ray.start.z / ray.dir.z);
        Point intersectingPoint = ray.start + ray.dir * t;

        double x = intersectingPoint.x - origin.x;
        double y = intersectingPoint.y - origin.y;
        if ((x < 0.0 || x > floorWidth) || (y < 0.0 || y > floorWidth))
            return -1;

        int pixelX = (intersectingPoint.x - origin.x) / tileWidth;
        int pixelY = (intersectingPoint.y - origin.y) / tileWidth;
        if ((pixelX < 0.0 || pixelX > numberOfTiles) || (pixelY < 0.0 || pixelY > numberOfTiles))
            return -1;
        int col = (pixelX + pixelY) % 2;
        color[0] = double(col);
        color[1] = double(col);
        color[2] = double(col);
        
        return t;
    }
};
typedef struct{double r, g, b;} Color;

class Light {
public:
    Point lightPos;
    double color[3];

    Light(Point pos, double red, double green, double blue) {
        lightPos = pos;
        color[0] = red;
        color[1] = green;
        color[2] = blue;
    }

    void draw() {
        glColor3f(color[0], color[1], color[2]);
        glPushMatrix();
        glTranslatef(lightPos.x, lightPos.y, lightPos.z);
        glutSolidSphere(1.0, 30, 30);
        glPopMatrix();

        glBegin(GL_LINES);
        {
            for (int i = -1; i <= 1; i++)
                for (int j = -1; j <= 1; j++)
                    for (int k = -1; k <= 1; k++) {
                        glVertex3f(lightPos.x, lightPos.y, lightPos.z);
                        glVertex3f(lightPos.x + i * 2.0, lightPos.y + j * 2.0, lightPos.z + k * 2.0);
                    }
        }
        glEnd();
    }
};


vector<Object *> objects;
vector<Light *> lights;


pair<double, double> getNearest(Ray &ray) {
    int nearest = -1;
    double minimum = 99999999;

    for (int i = 0; i < objects.size(); i++) {
        double j = objects[i]->getIntersection(ray);

        if (j <= 0) {
            continue;
        } else if (j < minimum) {
            minimum = j;
            nearest = i;
        }
    }

    return make_pair(nearest, minimum);
}

void Object::intersect(Ray &ray, double t, double currentColor[3], int level) {
    Point intersectingPoint = ray.start + ray.dir * t;
    Point normal = getNormal(intersectingPoint);
    Point reflection = reflectedRayDirection(ray, normal);
    Point refraction = refractedRayDirection(ray, normal);

    for (int i = 0; i < lights.size(); i++) {
        double ambient = coEfficients[AMBIENT], lambert = 0.0, phong = 0.0;

        Point dir = (lights[i]->lightPos - intersectingPoint).normalize();

        Point start = intersectingPoint + dir * 0.0000001;

        Ray L(start, dir);
        Point R = (normal * (dotProduct(L.dir, normal) * 2.0) - L.dir).normalize();
        Point V = (intersectingPoint * -1.0).normalize();

        bool boolean = false;

        int totalObjects = objects.size();
        for (int j = 0; j < totalObjects; j++) {
            double t = objects[j]->getIntersection(L);
            if (t > 0) {
                boolean = true;
                break;
            }
        }

        if (!boolean) {
            lambert = sourceFactor * coEfficients[DIFFUSE] * dotProduct(L.dir, normal);
            phong = coEfficients[SPECULAR] * pow(dotProduct(R, V), shine);

            if (lambert < 0.0) lambert = 0.0;
            if (lambert > 1.0) lambert = 1.0;


            if (phong < 0.0) phong = 0.0;
            if (phong > 1.0) phong = 1.0;
        }
        for (int k = 0; k < 3; k++) {
            currentColor[k] += ((ambient + lambert * lights[i]->color[k] + phong * lights[i]->color[k]) * color[k]);
        }

        if (level < recursionLevel) {
            Point start = intersectingPoint + reflection * 0.0000001;

            Ray reflectionRay(start, reflection);
            double reflectedColor[3] = {0.0, 0.0, 0.0};

            pair<double, double> pair = getNearest(reflectionRay);
            int nearest = pair.first;
            double t_min = pair.second;

            if (nearest != -1) {
                objects[nearest]->intersect(reflectionRay, t_min, reflectedColor, level + 1);

                for (int k = 0; k < 3; k++) {
                    currentColor[k] += reflectedColor[k] * coEfficients[REFLECTION];
                }
            }

            start = intersectingPoint + refraction * 0.0000001;

            Ray refractionRay(start, refraction);
            double refractedColor[3] = {0.0, 0.0, 0.0};

            pair = getNearest(refractionRay);
            nearest = pair.first;
            t_min = pair.second;

            if (nearest != -1) {
                objects[nearest]->intersect(refractionRay, t_min, refractedColor, level + 1);

                for (int k = 0; k < 3; k++) {
                    currentColor[k] += refractedColor[k] * eta;
                }
            }
        }

        for (int k = 0; k < 3; k++) {
            if (currentColor[k] < 0.0) currentColor[k] = 0.0;
            if (currentColor[k] > 1.0) currentColor[k] = 1.0;
        }
    }
}

#endif //RAY_TRACING_1605102_CLASSES_H

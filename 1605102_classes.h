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

    Point(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Point operator+(Point p) {
        return {this->x + p.x, this->y + p.y, this->z + p.z};
    }

    Point operator-(Point p) {
        return {this->x - p.x, this->y - p.y, this->z - p.z};
    }

    Point operator*(double scale) {
        return {this->x * scale, this->y * scale, this->z * scale};
    }

    Point operator/(double scale) {
        return {this->x / scale, this->y / scale, this->z / scale};
    }

    Point normalize() {
        return *this / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    }

    double length() const {
        return sqrt(x * x + y * y + z * z);
    }

    friend double dotProduct(Point p1, Point p2) {
        return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
    }

    friend Point crossProduct(Point p1, Point p2) {
        double a = p1.y * p2.z - p1.z * p2.y;
        double b = p1.z * p2.x - p1.x * p2.z;
        double c = p1.x * p2.y - p1.y * p2.x;

        return {a, b, c};
    }

    friend Point rotate(Point vec, Point ref, double directionSign) {
        Point temp, normalized ;
        double angle = 0.05 * directionSign;

        normalized.x = ref.y * vec.z - ref.z * vec.y;
        normalized.y = ref.z * vec.x - ref.x * vec.z;
        normalized.z = ref.x * vec.y - ref.y * vec.x;

        temp.x = vec.x * cos(angle) + normalized.x * sin(angle);
        temp.y = vec.y * cos(angle) + normalized.y * sin(angle);
        temp.z = vec.z * cos(angle) + normalized.z * sin(angle);

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
    char type = 'x';
    double sourceFactor = 1.0, eta = 0.0;
    Point referencePoint;
    double height, width, length;
    double color[3];
    double coEfficients[4];
    double shine;

    Object() {}

    virtual void draw() = 0;
    
    void setColor(double r, double g, double b) {
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }
    void setShine(int shine) {
        this->shine = shine;
    }
    void setCoefficients(double a, double d, double s, double r) {
        coEfficients[AMBIENT] = a;
        coEfficients[DIFFUSE] = d;
        coEfficients[SPECULAR] = s;
        coEfficients[REFLECTION] = r;
    }
    
    virtual double getIntersectionT(Ray &ray) = 0;

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

    double getIntersectionT(Ray &ray) {

        Point c = referencePoint;
        Point o = ray.start;
        Point l = ray.dir;

        Point d = o - c;
        double discriminant = dotProduct(l, d) * dotProduct(l, d) - dotProduct(d, d) + length * length;

        if (discriminant < 0)
            return -1;

        double sqrt_disc = sqrt(discriminant);
        double t1 = -dotProduct(l, d) + sqrt_disc;
        double t2 = -dotProduct(l, d) - sqrt_disc;

        return min(t1, t2);
    }

    Point getNormal(Point intersectingPoint) {
        Point normal = intersectingPoint - referencePoint;
        return normal.normalize();
    }
};


class Triangle : public Object {
public:
    Point p1, p2, p3;
    double a, b, c, d;

    Triangle(Point p1, Point p2, Point p3) {
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;

        Point AB = p2 - p1;
        Point AC = p3 - p1;
        Point temp = crossProduct(AB, AC);

        this->a = temp.x;
        this->b = temp.y;
        this->c = temp.z;
        this->d = this->a * p1.x + this->b * p1.y + this->c * p1.z;
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

    double getIntersectionT(Ray &ray) {
        Point normal = getNormal(Point(0, 0, 0));

        double t;
        bool flag = false;
        double denominator = dotProduct(normal, ray.dir);

        if (denominator < 0.0) {
            normal = normal * -1.0;
            denominator = dotProduct(normal, ray.dir);
        }
        
        if (abs(denominator) < 0.0000001)
            return -1;

        t = (dotProduct(normal, p1 - ray.start)) / denominator;
        if (t >= 0)
            flag = true;

        if (!flag)
            return -1;

        bool b1, b2, b3;
        Point intersectingPoint = ray.start + ray.dir * t;

        Point N = crossProduct(p2 - p1, p3 - p1);
        double area2 = N.length();

        Point A3, edge1, edge2;

        edge1 = p2 - p1;
        edge2 = intersectingPoint - p1;
        A3 = crossProduct(edge1, edge2);
        if (dotProduct(N, A3) < 0)
            return -1;

        edge1 = p3 - p2;
        edge2 = intersectingPoint - p2;
        A3 = crossProduct(edge1, edge2);
        if (dotProduct(N, A3) < 0)
            return -1;

        edge1 = p1 - p3;
        edge2 = intersectingPoint - p3;
        A3 = crossProduct(edge1, edge2);
        if (dotProduct(N, A3) < 0)
            return -1;

        return t; 
    }
};

class GeneralQuadratic : public Object {

public:
    double A1, A2, A3, A4, A5, A6, A7, A8, A9, A10;

    GeneralQuadratic(double coeff[10], Point referencePoint, double length, double width, double height) {
        this->A1 = coeff[0];
        this->A2 = coeff[1];
        this->A3 = coeff[2];
        this->A4 = coeff[3];
        this->A5 = coeff[4];
        this->A6 = coeff[5];
        this->A7 = coeff[6];
        this->A8 = coeff[7];
        this->A9 = coeff[8];
        this->A10 = coeff[9];
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

    double getIntersectionT(Ray &ray) {
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

        bool flag1 = length > 0 && (x1 < x_min || x1 > x_max) ||
                     width > 0 && (y1 < y_min || y1 > y_max) ||
                     height > 0 && (z1 < z_min || z1 > z_max);

        bool flag2 = length > 0 && (x2 < x_min || x2 > x_max) ||
                     width > 0 && (y2 < y_min || y2 > y_max) ||
                     height > 0 && (z2 < z_min || z2 > z_max);

        if (flag1 && flag2)
            return -1;
        else if (flag1)
            return t2;
        else if (flag2)
            return t1;
        else
            return min(t1, t2);
    }
};

class Floor : public Object {

public:
    Point origin;
    double floorWidth, tileWidth;
    int numberOfTiles;
    bitmap_image texture;

    Floor(double floorWidth, double tileWidth) {
        this->floorWidth = floorWidth;
        this->tileWidth = tileWidth;
        this->origin = Point(-floorWidth / 2.0, -floorWidth / 2.0, 0.0);
        this->numberOfTiles = floorWidth / tileWidth;

        texture = bitmap_image("1605102_bitmap.bmp");
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

    double getIntersectionT(Ray &ray) {
        if (ray.dir.z == 0)
            return -1;

        double t = -(ray.start.z / ray.dir.z);
        Point intersectingPoint = ray.start + ray.dir * t;

        double x = intersectingPoint.x - origin.x;
        double y = intersectingPoint.y - origin.y;
        if ((x < 0.0 || x > floorWidth) || (y < 0.0 || y > floorWidth))
            return -1;

        int pixel_x = (intersectingPoint.x - origin.x) / tileWidth;
        int pixel_y = (intersectingPoint.y - origin.y) / tileWidth;
        if ((pixel_x < 0.0 || pixel_x > numberOfTiles) || (pixel_y < 0.0 || pixel_y > numberOfTiles))
            return -1;
        int c = (pixel_x + pixel_y) % 2;

        unsigned char r, g, b;
        int i, j;
        i = (texture.width() - 1.0) / 1000.0 * x;
        j = (texture.height() - 1.0) / 1000.0 * y;
        texture.get_pixel(i, texture.height() - j - 1, r, g, b);

        double tile_portion = 0.8;
        double texture_portion = 1.0 - tile_portion;

        color[0] = double(c) * tile_portion + (double(r) / 255.0) * texture_portion;
        color[1] = double(c) * tile_portion + (double(g) / 255.0) * texture_portion;
        color[2] = double(c) * tile_portion + (double(b) / 255.0) * texture_portion;

        return t;
    }
};
typedef struct{double r, g, b;} Color;

class Light {
public:
    Point lightPos;
    double color[3];

    Light(Point pos, double r, double g, double b) {
        this->lightPos = pos;
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void draw() {
        glColor3f(color[0], color[1], color[2]);
        glPushMatrix();
        glTranslatef(lightPos.x, lightPos.y, lightPos.z);
        glutSolidSphere(1.0, 50, 50);
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
        double j = objects[i]->getIntersectionT(ray);

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

        bool flag = false;

        int totalObjects = objects.size();
        for (int j = 0; j < totalObjects; j++) {
            double t = objects[j]->getIntersectionT(L);
            if (t > 0) {
                flag = true;
                break;
            }
        }

        if (!flag) {
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

//
// Created by nisa on 25/06/21.
//

#ifndef RAY_TRACING_1605102_CLASSES_H
#define RAY_TRACING_1605102_CLASSES_H


#define amb 0
#define dif 1
#define spec 2
#define refl 3
extern int recursion_level;

class Point
{
public:
    double x, y, z;

    Point() {}

    Point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Point operator+(Point p)
    {
        return Point(this->x + p.x, this->y + p.y, this->z + p.z);
    }

    Point operator-(Point p)
    {
        return Point(this->x - p.x, this->y - p.y, this->z - p.z);
    }

    Point operator*(double scale)
    {
        return Point(this->x * scale, this->y * scale, this->z * scale);
    }

    Point operator/(double scale)
    {
        return Point(this->x / scale, this->y / scale, this->z / scale);
    }

    Point normalize()
    {
        return * this / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    }

    double length()
    {
        return sqrt(x * x + y * y + z * z);
    }

    friend double dot(Point p1, Point p2)
    {
        return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
    }

    friend Point cross(Point p1, Point p2)
    {
        double a = p1.y * p2.z - p1.z * p2.y;
        double b = p1.z * p2.x - p1.x * p2.z;
        double c = p1.x * p2.y - p1.y * p2.x;

        return Point(a, b, c);
    }

    friend Point rotate(Point vec, Point ref, double sign)
    {
        Point n, ret;
        double angle = 0.05 * sign;

        n.x = ref.y * vec.z - ref.z * vec.y;
        n.y = ref.z * vec.x - ref.x * vec.z;
        n.z = ref.x * vec.y - ref.y * vec.x;

        ret.x = vec.x * cos(angle) + n.x * sin(angle);
        ret.y = vec.y * cos(angle) + n.y * sin(angle);
        ret.z = vec.z * cos(angle) + n.z * sin(angle);

        return ret;
    }

    void print()
    {
        cout << x << ", " << y << ", " << z << endl;
    }
};


class Ray
{

public:
    Point start;
    Point dir;

    Ray(Point start, Point dir)
    {
        this->start = start;
        this->dir = dir.normalize();
    }
};


class Object
{

public:
    char type = 'x';
    double source_factor = 1.0, eta = 0.0;
    double shine;
    double color[3];
    double co_efficients[4];

    Object() {}

    virtual void draw() = 0;
    virtual double getIntersectionT(Ray &ray) = 0;
    virtual Point getNormal(Point intersection) = 0;

    void fill_color(Ray &ray, double t, double current_color[3], int level);

    void setColor(double r, double g, double b)
    {
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void setShine(int shine)
    {
        this->shine = shine;
    }

    void set_lighting_coefficients(double a, double d, double s, double r)
    {
        co_efficients[amb] = a;
        co_efficients[dif] = d;
        co_efficients[spec] = s;
        co_efficients[refl] = r;
    }

    Point get_reflected_ray_direction(Ray &ray, Point normal)
    {
        Point out_dir = ray.dir - normal * 2.0 * dot(ray.dir, normal);
        return out_dir.normalize();
    }

    Point get_refracted_ray_direction(Ray &ray, Point normal)
    {
        double n_dot_i = dot(normal, ray.dir);
        double k = 1.0 - eta * eta * (1.0 - n_dot_i * n_dot_i);

        if (k < 0.0)
            return Point(0.0, 0.0, 0.0);

        Point out_dir = ray.dir * eta - normal * (eta * n_dot_i + sqrt(k));
        return out_dir.normalize();
    }
};



class Sphere : public Object
{
public:
    Point center;
    double radius;

    Sphere(Point center, double radius)
    {
        this->center = center;
        this->radius = radius;
    }

    void draw()
    {
        glColor3f(color[0], color[1], color[2]);
        glPushMatrix();
        glTranslatef(center.x, center.y, center.z);
        glutSolidSphere(radius, 100, 100);
        glPopMatrix();
    }

    double getIntersectionT(Ray &ray)
    {

        Point c = center;
        Point o = ray.start;
        Point l = ray.dir;

        Point d = o - c;
        double discriminant = dot(l, d) * dot(l, d) - dot(d, d) + radius * radius;

        if (discriminant < 0)
            return -1;

        double sqrt_disc = sqrt(discriminant);
        double t1 = -dot(l, d) + sqrt_disc;
        double t2 = -dot(l, d) - sqrt_disc;

        return min(t1, t2);
    }

    Point getNormal(Point intersection)
    {
        Point normal = intersection - center;
        return normal.normalize();
    }
};


class Triangle : public Object
{
public:
    Point p1, p2, p3;
    // let the equation of the plane of the triangle be ax + by + cz = d
    double a, b, c, d;

    Triangle(Point p1, Point p2, Point p3)
    {
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;

        // calculating a, b, c, d
        // using parametric form
        Point AB = p2 - p1;
        Point AC = p3 - p1;

        Point V = cross(AB, AC);

        this->a = V.x;
        this->b = V.y;
        this->c = V.z;
        this->d = this->a * p1.x + this->b * p1.y + this->c * p1.z;
    }

    void draw()
    {
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(p1.x, p1.y, p1.z);
            glVertex3f(p2.x, p2.y, p2.z);
            glVertex3f(p3.x, p3.y, p3.z);
        }
        glEnd();
    }

    Point getNormal(Point intersection)
    {
        Point normal = Point(a, b, c);
        return normal.normalize();
    }

    double getIntersectionT(Ray &ray)
    {
        Point normal = getNormal(Point(0, 0, 0));

        double t;
        bool flag = false;
        double denom = dot(normal, ray.dir);

        if (denom < 0.0)
        {
            normal = normal * -1.0;
            denom = dot(normal, ray.dir);
        }

        // checking whether the ray intersects with the triangle's plane
        if (abs(denom) < 0.0000001)
            return -1;

        t = (dot(normal, p1 - ray.start)) / denom;
        if (t >= 0)
            flag = true; // yes, intersects the plane

        if (!flag)
            return -1; // no

        // checking whether the point of intersection is inside the triangle
        bool b1, b2, b3;
        Point intersectionPoint = ray.start + ray.dir * t;

        Point N = cross(p2 - p1, p3 - p1);
        double area2 = N.length();

        Point C, edge1, edge2;

        // checking if intersectionPoint is inside the triangle

        edge1 = p2 - p1;
        edge2 = intersectionPoint - p1;
        C = cross(edge1, edge2);
        if (dot(N, C) < 0)
            return -1;

        edge1 = p3 - p2;
        edge2 = intersectionPoint - p2;
        C = cross(edge1, edge2);
        if (dot(N, C) < 0)
            return -1;

        edge1 = p1 - p3;
        edge2 = intersectionPoint - p3;
        C = cross(edge1, edge2);
        if (dot(N, C) < 0)
            return -1;

        return t; // else, yes
    }
};

class GeneralQuadratic : public Object
{

public:
    Point reference_point;
    double height, width, length;
    double A, B, C, D, E, F, G, H, I, J;

    GeneralQuadratic(double coeff[10], Point reference_point, double length, double width, double height)
    {
        this->A = coeff[0];
        this->B = coeff[1];
        this->C = coeff[2];
        this->D = coeff[3];
        this->E = coeff[4];
        this->F = coeff[5];
        this->G = coeff[6];
        this->H = coeff[7];
        this->I = coeff[8];
        this->J = coeff[9];
        this->reference_point = reference_point;
        this->height = height;
        this->width = width;
        this->length = length;
    }

    void draw() {}

    Point getNormal(Point intersection)
    {

        double u = 2 * A * intersection.x + D * intersection.y + F * intersection.z + G;
        double v = 2 * B * intersection.y + D * intersection.x + E * intersection.z + H;
        double z = 2 * C * intersection.z + E * intersection.y + F * intersection.x + I;

        Point normal(u, v, z);
        return normal.normalize();
    }

    double getIntersectionT(Ray &ray)
    {
        double a = A * ray.dir.x * ray.dir.x + B * ray.dir.y * ray.dir.y + C * ray.dir.z * ray.dir.z;
        double b = 2 * (A * ray.start.x * ray.dir.x + B * ray.start.y * ray.dir.y + C * ray.start.z * ray.dir.z);
        double c = A * ray.start.x * ray.start.x + B * ray.start.y * ray.start.y + C * ray.start.z * ray.start.z;

        a += (D * ray.dir.x * ray.dir.y + E * ray.dir.y * ray.dir.z + F * ray.dir.z * ray.dir.x);
        b += (D * (ray.start.x * ray.dir.y + ray.dir.x * ray.start.y) + E * (ray.start.y * ray.dir.z + ray.dir.y * ray.start.z) + F * (ray.start.z * ray.dir.x + ray.dir.z * ray.start.x));
        c += (D * ray.start.x * ray.start.y + E * ray.start.y * ray.start.z + F * ray.start.z * ray.start.x);

        b += (G * ray.dir.x + H * ray.dir.y + I * ray.dir.z);
        c += (G * ray.start.x + H * ray.start.y + I * ray.start.z + J);

        double discriminant = b * b - 4 * a * c;

        if (discriminant < 0)
        {
            return -1;
        }

        double t1 = (-b + sqrt(discriminant)) / (2.0 * a);
        double t2 = (-b - sqrt(discriminant)) / (2.0 * a);

        Point intersectionPoint1 = ray.start + ray.dir * t1;
        Point intersectionPoint2 = ray.start + ray.dir * t2;

        double x_min = reference_point.x;
        double x_max = x_min + length;

        double y_min = reference_point.y;
        double y_max = y_min + width;

        double z_min = reference_point.z;
        double z_max = z_min + height;

        double x1 = intersectionPoint1.x;
        double y1 = intersectionPoint1.y;
        double z1 = intersectionPoint1.z;

        double x2 = intersectionPoint2.x;
        double y2 = intersectionPoint2.y;
        double z2 = intersectionPoint2.z;

        bool flag1 =    length > 0 && (x1 < x_min || x1 > x_max) ||
                        width > 0 && (y1 < y_min || y1 > y_max) ||
                        height > 0 && (z1 < z_min || z1 > z_max);

        bool flag2 =    length > 0 && (x2 < x_min || x2 > x_max) ||
                        width > 0 && (y2 < y_min || y2 > y_max) ||
                        height > 0 && (z2 < z_min || z2 > z_max);

        if (flag1 && flag2) return -1;
        else if (flag1) return t2;
        else if (flag2) return t1;
        else return min(t1, t2);
    }
};

class Floor : public Object
{

public:
    Point origin;
    double floorWidth, tile_size;
    int n_tiles;
    bitmap_image texture;

    Floor(double floorWidth, double tile_size)
    {
        this->floorWidth = floorWidth;
        this->tile_size = tile_size;
        this->origin = Point(-floorWidth / 2.0, -floorWidth / 2.0, 0.0);
        this->n_tiles = floorWidth / tile_size;

        texture = bitmap_image("1605102_bitmap.bmp");
    }

    void draw()
    {
        glBegin(GL_QUADS);
        {
            for (int i = 0; i < n_tiles; i++)
                for (int j = 0; j < n_tiles; j++)
                {
                    bool c = (i + j) % 2;
                    glColor3f(c, c, c);

                    glVertex3f(origin.x + tile_size * i, origin.y + tile_size * j, origin.z);
                    glVertex3f(origin.x + tile_size * (i + 1), origin.y + tile_size * j, origin.z);
                    glVertex3f(origin.x + tile_size * (i + 1), origin.y + tile_size * (j + 1), origin.z);
                    glVertex3f(origin.x + tile_size * i, origin.y + tile_size * (j + 1), origin.z);
                }
        }
        glEnd();
    }

    Point getNormal(Point intersection) { return Point(0, 0, 1); }

    double getIntersectionT(Ray &ray)
    {
        if (ray.dir.z == 0)
            return -1;

        double t = -(ray.start.z / ray.dir.z);
        Point intersectionPoint = ray.start + ray.dir * t;

        double x = intersectionPoint.x - origin.x;
        double y = intersectionPoint.y - origin.y;
        if ((x < 0.0 || x > floorWidth) || (y < 0.0 || y > floorWidth))
            return -1;

        int pixel_x = (intersectionPoint.x - origin.x) / tile_size;
        int pixel_y = (intersectionPoint.y - origin.y) / tile_size;
        if ((pixel_x < 0.0 || pixel_x > n_tiles) || (pixel_y < 0.0 || pixel_y > n_tiles))
            return -1;
        int c = (pixel_x + pixel_y) % 2;

        unsigned char r, g, b;
        int i, j;
        i = (texture.width()-1.0)/1000.0 * x;
        j = (texture.height()-1.0)/1000.0 * y;
        texture.get_pixel(i, texture.height()-j-1, r, g, b);

        double tile_portion = 0.7;
        double texture_portion = 1.0 - tile_portion;

        color[0] = double(c)*tile_portion + (double(r) / 255.0)*texture_portion;
        color[1] = double(c)*tile_portion + (double(g) / 255.0)*texture_portion;
        color[2] = double(c)*tile_portion + (double(b) / 255.0)*texture_portion;

        return t;
    }
};


vector<Object *> objects;
vector<Point> lights;


pair<double, double> get_nearest(Ray &ray)
{
    int nearest = -1;
    double t_min = 9999999;

    for (int k = 0; k < objects.size(); k++)
    {
        double tk = objects[k]->getIntersectionT(ray);

        if (tk <= 0)
        {
            continue;
        }
        else if (tk < t_min)
        {
            t_min = tk;
            nearest = k;
        }
    }

    return make_pair(nearest, t_min);
}

void Object::fill_color(Ray &ray, double t, double current_color[3], int level)
{
    Point intersectionPoint = ray.start + ray.dir * t;
    Point normal = getNormal(intersectionPoint);
    Point reflection = get_reflected_ray_direction(ray, normal);
    Point refraction = get_refracted_ray_direction(ray, normal);

    for (int i = 0; i < lights.size(); i++)
    {
        double ambient = co_efficients[amb], lambert = 0.0, phong = 0.0;

        Point dir = (lights[i] - intersectionPoint).normalize();

        Point start = intersectionPoint + dir * 0.0000001;
        // adding epsilon so that the image itself is not detected as a blockade 
        Ray L(start, dir);
        Point R = (normal * (dot(L.dir, normal) * 2.0) - L.dir).normalize();
        Point V = (intersectionPoint * -1.0).normalize();

        bool flag = false;

        int n_objects = objects.size();
        for (int j = 0; j < n_objects; j++)
        {
            double t = objects[j]->getIntersectionT(L);
            if (t > 0)
            {
                flag = true;
                break;
            }
        }

        if (!flag)
        {
            lambert = source_factor * co_efficients[dif] * dot(L.dir, normal);
            phong = co_efficients[spec] * pow(dot(R, V), shine);

            if (lambert < 0.0) lambert = 0.0;
            if (lambert > 1.0) lambert = 1.0;


            if (phong < 0.0) phong = 0.0;
            if (phong > 1.0) phong = 1.0;
        }
        for (int k = 0; k < 3; k++)
        {
            current_color[k] += ((ambient + lambert + phong) * color[k]);
        }

        if (level < recursion_level)
        {
            Point start = intersectionPoint + reflection * 0.0000001;

            Ray reflectionRay(start, reflection);
            double reflected_color[3] = {0.0, 0.0, 0.0};

            pair<double, double> pair = get_nearest(reflectionRay);
            int nearest = pair.first;
            double t_min = pair.second;

            if (nearest != -1)
            {
                objects[nearest]->fill_color(reflectionRay, t_min, reflected_color, level + 1);

                for (int k = 0; k < 3; k++)
                {
                    current_color[k] += reflected_color[k] * co_efficients[refl];
                }
            }

            start = intersectionPoint + refraction * 0.0000001;

            Ray refractionRay(start, refraction);
            double refracted_color[3] = {0.0, 0.0, 0.0};

            pair = get_nearest(refractionRay);
            nearest = pair.first;
            t_min = pair.second;

            if (nearest != -1)
            {
                objects[nearest]->fill_color(refractionRay, t_min, refracted_color, level + 1);

                for (int k = 0; k < 3; k++)
                {
                    current_color[k] += refracted_color[k] * eta;
                }
            }
        }

        for (int k = 0; k < 3; k++)
        {
            if (current_color[k] < 0.0) current_color[k] = 0.0;
            if (current_color[k] > 1.0) current_color[k] = 1.0;
        }
    }
}

#endif //RAY_TRACING_1605102_CLASSES_H

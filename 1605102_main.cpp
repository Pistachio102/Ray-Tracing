//
// Created by nisa on 25/06/21.
//

#include <bits/stdc++.h>
#include <GL/glut.h>

using namespace std;

#include "bitmap_image.hpp"
#include "1605102_classes.h"

#define pi (2*acos(0.0))
#define unit_dist 4.0


Point eye, l, r, u;
int recursion_level;
int image_width, image_height;

void drawAxes() {
    glBegin(GL_LINES);
    {
        glVertex3f( 100,0,0);
        glVertex3f(-100,0,0);

        glVertex3f(0,-100,0);
        glVertex3f(0, 100,0);

        glVertex3f(0,0, 100);
        glVertex3f(0,0,-100);
    }
    glEnd();
}

void createFloor()
{
    Object *floor = new Floor(1000, 20);
    floor->set_lighting_coefficients(0.3, 0.3, 0.3, 0.1);
    floor->setShine(1.0);
    objects.push_back(floor);
}

void loadData() {

    int n_objects, n_lights;
    string object_type;
    Object *obj;

    freopen("scene.txt", "r", stdin);

    cin >> recursion_level >> image_width >> n_objects;
    image_height = image_width;

    for (int i = 0; i < n_objects; i++)
    {
        cin >> object_type;

        if (object_type == "sphere") {

            double x, y, z; // center
            double radius; // radius
            double red, green, blue;
            double ambient, diffuse, specular, reflection;
            double shine;

            Point center;

            cin >> x >> y >> z >> radius;
            center = Point(x, y, z);
            obj = new Sphere(center, radius);

            cin >> red >> green >> blue;
            obj->setColor(red, green, blue);

            cin >> ambient >> diffuse >> specular >> reflection;
            obj->set_lighting_coefficients(ambient, diffuse, specular, reflection);

            cin >> shine;
            obj->setShine(shine);

            obj->eta = 10.0;

            objects.push_back(obj);
        }

        else if (object_type == "triangle") {

            double x, y, z; // a vertex
            double red, green, blue;
            double ambient, diffuse, specular, reflection;
            double shine;

            Point A, B, C;

            cin >> x >> y >> z;
            A = Point(x, y, z);

            cin >> x >> y >> z;
            B = Point(x, y, z);

            cin >> x >> y >> z;
            C = Point(x, y, z);

            obj = new Triangle(A, B, C);

            cin >> red >> green >> blue;
            obj->setColor(red, green, blue);

            cin >> ambient >> diffuse >> specular >> reflection;
            obj->set_lighting_coefficients(ambient, diffuse, specular, reflection);

            cin >> shine;
            obj->setShine(shine);

            objects.push_back(obj);

        }

        else if (object_type == "general") {

            double coeff[10];
            double x, y, z;
            double length, width, height;
            double red, green, blue;
            double ambient, diffuse, specular, reflection;
            double shine;

            Point base_point;

            for (int c = 0; c < 10; c++) {
                cin>>coeff[c];
            }

            cin >> x >> y >> z;
            base_point = Point(x, y, z);

            cin >> length >> width >> height;
            obj = new GeneralQuadratic(coeff, base_point, length, width, height);

            cin >> red >> green >> blue;
            obj->setColor(red, green, blue);

            cin >> ambient >> diffuse >> specular >> reflection;
            obj->set_lighting_coefficients(ambient, diffuse, specular, reflection);

            cin >> shine;
            obj->setShine(shine);

            objects.push_back(obj);
        }
    }

    cin >> n_lights;
    for (int i = 0; i < n_lights; i++) {
        double x, y, z;
        cin >> x >> y >> z;

        Point light(x, y, z);
        lights.push_back(light);
    }
}

void keyboardListener(unsigned char key, int x, int y)
{
    switch (key)
    {
        case '1':
            l = rotate(l, u, 1);
            r = rotate(r, u, 1);
            break;
        case '2':
            l = rotate(l, u, -1);
            r = rotate(r, u, -1);
            break;
        case '3':
            u = rotate(u, r, 1);
            l = rotate(l, r, 1);
            break;
        case '4':
            u = rotate(u, r, -1);
            l = rotate(l, r, -1);
            break;
        case '5':
            u = rotate(u, l, 1);
            r = rotate(r, l, 1);
            break;
        case '6':
            u = rotate(u, l, -1);
            r = rotate(r, l, -1);
            break;

//        case '0':
//            capture();
//            break;
        default:
            break;
    }
}


void specialKeyListener(int key, int x, int y)
{
    switch(key) {
        case GLUT_KEY_DOWN:
            eye = eye - l * unit_dist;
            break;
        case GLUT_KEY_UP:
            eye = eye + l * unit_dist;
            break;
        case GLUT_KEY_RIGHT:
            eye = eye + r * unit_dist;
            break;
        case GLUT_KEY_LEFT:
            eye = eye - r * unit_dist;
            break;
        case GLUT_KEY_PAGE_UP:
            eye = eye + u * unit_dist;
            break;
        case GLUT_KEY_PAGE_DOWN:
            eye = eye - u * unit_dist;
            break;
        case GLUT_KEY_INSERT:
            break;

        default:
            break;
    }
}


void mouseListener(int button, int state, int x, int y) {	//x, y is the x-y of the screen (2D)
    switch (button) {
        case GLUT_LEFT_BUTTON:
            if (state == GLUT_DOWN) {		// 2 times?? in ONE click? -- solution is checking DOWN or UP
                //drawaxes = 1 - drawaxes;
            }
            break;

        case GLUT_RIGHT_BUTTON:
            //........
            break;

        case GLUT_MIDDLE_BUTTON:
            //........
            break;

        default:
            break;
    }
}


void display() {

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0);    //color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?
    gluLookAt(eye.x, eye.y, eye.z, eye.x + l.x, eye.y + l.y, eye.z + l.z, u.x, u.y, u.z);

    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();


    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}


void animate() {
    //codes for any changes in Models, Camera

    glutPostRedisplay();
}


void init() {
    //codes for initialization
    eye = {0, -200, 30};
    l = { 0.0, 1.0, 0.0 };
    r = { 1.0, 0.0, 0.0 };
    u = { 0.0, 0.0, 1.0 };

    createFloor();
    loadData();

    //clear the screen
    glClearColor(0, 0, 0, 0);


    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(80, 1, 1, 1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}


int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);    //Depth, Double buffer, RGB color

    glutCreateWindow("My OpenGL Program");

    init();

    glEnable(GL_DEPTH_TEST);    //enable Depth Testing

    glutDisplayFunc(display);    //display callback function
    glutIdleFunc(animate);        //what you want to do in the idle time (when no drawing is occurring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();        //The main loop of OpenGL

    return 0;
}
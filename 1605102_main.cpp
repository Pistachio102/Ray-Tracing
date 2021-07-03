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
int recursionLevel;
int imageWidth, imageHeight;
#define WINDOW_WIDTH 700
#define WINDOW_HEIGHT 700
#define VIEW_ANGLE 50
void drawAxes() {
    glBegin(GL_LINES);
    {
        glVertex3f(100, 0, 0);
        glVertex3f(-100, 0, 0);

        glVertex3f(0, -100, 0);
        glVertex3f(0, 100, 0);

        glVertex3f(0, 0, 100);
        glVertex3f(0, 0, -100);
    }
    glEnd();
}

void createFloor() {
    Object *floor = new Floor(1000, 20);
    floor->set_lighting_coefficients(0.3, 0.3, 0.3, 0.1);
    floor->setShine(1.0);
    objects.push_back(floor);
}

void loadData() {

    int totalObjects, totalLights;
    string objectType;
    Object *obj;

    freopen("scene.txt", "r", stdin);

    cin >> recursionLevel >> imageWidth >> totalObjects;
    imageHeight = imageWidth;
    cout << recursionLevel << imageWidth << totalObjects << endl;

    for (int i = 0; i < totalObjects; i++) {
        cin >> objectType;
        cout << objectType << endl;

        if (objectType == "sphere") {

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
        } else if (objectType == "triangle") {

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

        } else if (objectType == "general") {

            double coeff[10];
            double x, y, z;
            double length, width, height;
            double red, green, blue;
            double ambient, diffuse, specular, reflection;
            double shine;

            Point base_point;

            for (int c = 0; c < 10; c++) {
                cin >> coeff[c];
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

    cin >> totalLights;

    for (int i = 0; i < totalLights; i++) {
        double x, y, z, r, g, b;
        cin >> x >> y >> z;
        cin >> r >> g >> b;

        Point pos(x, y, z);
        Light *light = new Light(pos, r, g, b);
        lights.push_back(light);
    }
}

void capture() {

    Color** frameBuffer;
    frameBuffer = new Color* [imageWidth];
    for (int i=0; i<imageWidth; i++) {
        frameBuffer[i] = new Color[imageHeight];
    }

    double plane_distance = (WINDOW_HEIGHT/2)/tan(VIEW_ANGLE*pi/360);
    Point top_left = eye + (l * plane_distance - r * (WINDOW_WIDTH/2) + u * (WINDOW_HEIGHT/2));

    cout << "Eye : "; eye.print();
    cout << "Plane distance : " << plane_distance << endl;
    cout << "top_left : "; top_left.print();
    cout << "Saving...";

    double du = (WINDOW_WIDTH*1.0) / imageWidth;
    double dv = (WINDOW_HEIGHT*1.0) / imageHeight;

    for (int i = 0; i < imageWidth; i++) {
        for (int j = 0; j < imageHeight; j++) {

            Point direction_to_top_left = top_left + r*i*du - u*j*dv;

            Ray ray(eye, direction_to_top_left - eye);
            double dummy_color[3] = {0.0, 0.0, 0.0};

            pair<double, double> pair = get_nearest(ray);
            int nearest = pair.first;
            double t_min = pair.second;

            if(nearest!=-1) {
                objects[nearest]->fill_color(ray, t_min, dummy_color, 1);
            }
            frameBuffer[i][j].r = dummy_color[0];
            frameBuffer[i][j].g = dummy_color[1];
            frameBuffer[i][j].b = dummy_color[2];
        }
    }

    bitmap_image image(imageWidth, imageHeight);

    for (int i=0; i<imageWidth; i++) {
        for (int j=0; j<imageHeight; j++) {
            double r = frameBuffer[i][j].r;
            double g = frameBuffer[i][j].g;
            double b = frameBuffer[i][j].b;
            image.set_pixel(i, j, r*255, g*255, b*255);
        }
    }

    image.save_image("out.bmp");

    cout << "\tSaved\n";
    cout << "\a";
}

void freeMemory() {
    vector<Light*>().swap(lights);
    vector<Object*>().swap(objects);
}
void keyboardListener(unsigned char key, int x, int y) {
    switch (key) {
        case '0':
            capture();
            break;
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


        default:
            break;
    }
}


void specialKeyListener(int key, int x, int y) {
    switch (key) {
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


void mouseListener(int button, int state, int x, int y) {    //x, y is the x-y of the screen (2D)
    switch (button) {
        case GLUT_LEFT_BUTTON:
            if (state == GLUT_DOWN) {        // 2 times?? in ONE click? -- solution is checking DOWN or UP
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

    for (int i = 0; i < objects.size(); i++)
        objects[i]->draw();

    for (int i = 0; i < lights.size(); i++)
        lights[i]->draw();


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
    l = {0.0, 1.0, 0.0};
    r = {1.0, 0.0, 0.0};
    u = {0.0, 0.0, 1.0};

    createFloor();
    loadData();
//    cout << lights[0].lightPos.x << lights[0].lightPos.y << lights[0].lightPos.z << endl;
//    cout << lights[1].lightPos.x << lights[1].lightPos.y << lights[1].lightPos.z << endl;
//    cout << lights[2].lightPos.x << lights[2].lightPos.y << lights[2].lightPos.z << endl;
//    cout << lights[3].lightPos.x << lights[3].lightPos.y << lights[3].lightPos.z << endl;

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
//
// Created by nisa on 25/06/21.
//

#include <bits/stdc++.h>
#include <GL/glut.h>

using namespace std;

#include "bitmap_image.hpp"
#include "1605102_classes.h"

#define pi (2*acos(0.0))
#define unitDist 4.0


Point eye, l, r, u;
int recursionLevel;
int imageWidth, imageHeight;
#define windowWidth 1000
#define windowHeight 1000
#define viewAngle 60


void drawAxes() {
    glBegin(GL_LINES);
    {
        glColor3f(1, 1, 1);
        glVertex3f(100, 0, 0);
        glVertex3f(-100, 0, 0);

        glVertex3f(0, -100, 0);
        glVertex3f(0, 100, 0);

        glVertex3f(0, 0, 100);
        glVertex3f(0, 0, -100);
    }
    glEnd();
}



void loadData() {

    Object *temp = new Floor(1000, 20);
    temp->setCoefficients(0.2, 0.2, 0.2, 0.2);
    temp->setShine(1.5);
    objects.push_back(temp);

    int totalObjects, totalLights;
    string objectType;
    Object *object;

    freopen("scene.txt", "r", stdin);

    cin >> recursionLevel >> imageWidth >> totalObjects;
    imageHeight = imageWidth;
    cout << recursionLevel << imageWidth << totalObjects << endl;

    for (int i = 0; i < totalObjects; i++) {
        cin >> objectType;

        if (objectType == "sphere") {
            Point center;
            double x, y, z, radius;
            double r, g, b;
            double ambient, diffuse, specular, reflection;
            double shine;

            cin >> x >> y >> z >> radius;
            cin >> r >> g >> b;
            cin >> ambient >> diffuse >> specular >> reflection;
            cin >> shine;

            //create sphere
            center = Point(x, y, z);
            object = new Sphere(center, radius);
            //set color
            object->setColor(r, g, b);
            //set coefficients
            object->setCoefficients(ambient, diffuse, specular, reflection);
            //set shine
            object->setShine(shine);
            
            objects.push_back(object);
            
        } else if (objectType == "triangle") {

            double x, y, z;
            double r, g, b;
            double ambient, diffuse, specular, reflection;
            double shine;

            Point a1, a2, a3;

            cin >> x >> y >> z;
            a1 = Point(x, y, z);

            cin >> x >> y >> z;
            a2 = Point(x, y, z);

            cin >> x >> y >> z;
            a3 = Point(x, y, z);

            object = new Triangle(a1, a2, a3);

            cin >> r >> g >> b;
            cin >> ambient >> diffuse >> specular >> reflection;
            cin >> shine;
            
            //set color
            object->setColor(r, g, b);
            //set coefficients
            object->setCoefficients(ambient, diffuse, specular, reflection);
            //set shine
            object->setShine(shine);

            objects.push_back(object);

        } else if (objectType == "general") {

            double coefficients[10];
            Point ori;
            double x, y, z;
            double length, width, height;
            double r, g, b;
            double ambient, diffuse, specular, reflection;
            double shine;


            for (int coeff = 0; coeff < 10; coeff++) {
                cin >> coefficients[coeff];
            }

            cin >> x >> y >> z;
            cin >> length >> width >> height;
            ori = Point(x, y, z);

            object = new GeneralQuadratic(coefficients, ori, width, height, length);

            cin >> r >> g >> b;
            cin >> ambient >> diffuse >> specular >> reflection;
            cin >> shine;

            //set color
            object->setColor(r, g, b);
            //set coefficients
            object->setCoefficients(ambient, diffuse, specular, reflection);
            //set shine
            object->setShine(shine);

            objects.push_back(object);
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
    cout << "Saving..." << endl;
    Color** frameBuffer;
    frameBuffer = new Color* [imageWidth];
    for (int i=0; i<imageWidth; i++) {
        frameBuffer[i] = new Color[imageHeight];
    }

    double planeDistance = (windowHeight/2)/tan(viewAngle*pi/360);
    Point topLeft = eye + (l * planeDistance - r * (windowWidth/2) + u * (windowHeight/2));

    double du = (windowWidth*1.0) / imageWidth;
    double dv = (windowHeight*1.0) / imageHeight;

    for (int i = 0; i < imageWidth; i++) {
        for (int j = 0; j < imageHeight; j++) {

            Point dirTopLeft = topLeft + r*i*du - u*j*dv;

            Ray ray(eye, dirTopLeft - eye);
            double dummyColor[3] = {0.0, 0.0, 0.0};

            pair<double, double> pair = getNearest(ray);
            int nearest = pair.first;
            double t_min = pair.second;

            if(nearest!=-1) {
                objects[nearest]->intersect(ray, t_min, dummyColor, 1);
            }
            frameBuffer[i][j].r = dummyColor[0];
            frameBuffer[i][j].g = dummyColor[1];
            frameBuffer[i][j].b = dummyColor[2];
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
    image.save_image("1605102_out.bmp");
    cout << "Saved\n";
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
            l = afterRotation(l, u, 1);
            r = afterRotation(r, u, 1);
            break;
        case '2':
            l = afterRotation(l, u, -1);
            r = afterRotation(r, u, -1);
            break;
        case '3':
            u = afterRotation(u, r, 1);
            l = afterRotation(l, r, 1);
            break;
        case '4':
            u = afterRotation(u, r, -1);
            l = afterRotation(l, r, -1);
            break;
        case '5':
            u = afterRotation(u, l, 1);
            r = afterRotation(r, l, 1);
            break;
        case '6':
            u = afterRotation(u, l, -1);
            r = afterRotation(r, l, -1);
            break;


        default:
            break;
    }
}


void specialKeyListener(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_DOWN:
            eye = eye - l * unitDist;
            break;
        case GLUT_KEY_UP:
            eye = eye + l * unitDist;
            break;
        case GLUT_KEY_RIGHT:
            eye = eye + r * unitDist;
            break;
        case GLUT_KEY_LEFT:
            eye = eye - r * unitDist;
            break;
        case GLUT_KEY_PAGE_UP:
            eye = eye + u * unitDist;
            break;
        case GLUT_KEY_PAGE_DOWN:
            eye = eye - u * unitDist;
            break;
        case GLUT_KEY_INSERT:
            break;

        default:
            break;
    }
}


void mouseListener(int button, int state, int x, int y) {
    switch (button) {
        case GLUT_LEFT_BUTTON:
            //........
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

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0);    //color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye.x, eye.y, eye.z, eye.x + l.x, eye.y + l.y, eye.z + l.z, u.x, u.y, u.z);
    glMatrixMode(GL_MODELVIEW);

    drawAxes();

    for (int i = 0; i < objects.size(); i++)
        objects[i]->draw();

    for (int i = 0; i < lights.size(); i++)
        lights[i]->draw();


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


    loadData();
//    cout << lights[0].lightPos.x << lights[0].lightPos.y << lights[0].lightPos.z << endl;
//    cout << lights[1].lightPos.x << lights[1].lightPos.y << lights[1].lightPos.z << endl;
//    cout << lights[2].lightPos.x << lights[2].lightPos.y << lights[2].lightPos.z << endl;
//    cout << lights[3].lightPos.x << lights[3].lightPos.y << lights[3].lightPos.z << endl;

    //clear the screen
    glClearColor(0, 0, 0, 0);

    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    gluPerspective(80, 1, 1, 1000.0);

}


int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);    //Depth, Double buffer, RGB color

    glutCreateWindow("1605102: My OpenGL Program");

    init();

    glEnable(GL_DEPTH_TEST);    //enable Depth Testing

    glutDisplayFunc(display);    //display callback function
    glutIdleFunc(animate);        //what you want to do in the idle time (when no drawing is occurring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();        //The main loop of OpenGL
    freeMemory();
    return 0;
}
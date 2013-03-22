/*
	Scene.cpp
*/

#include <iostream>
#include <fstream>
#include <deque>
#include <time.h>

#include <glm/glm.hpp>
typedef glm::dvec3 dvec3; 
typedef glm::dvec4 dvec4; 
typedef glm::dmat4 dmat4; 

#include "Scene.h"
#include "Ray.h"
#include "Triangle.h"
#include "Sphere.h"
#include "Material.h"
#include "Transform.h"

using namespace std;

Scene::Scene(const char* filename) : camera(), screen()
{
    this->outputFilename = "screenshot.png";
    this->attenuation = dvec3(1.0,0.0,0.0);
    this->ambient = dvec3(0.2,0.2,0.2);
    this->samplesPerPixel = 1;
    this->samplerType = REGULAR;
    this->traceType = RAY;
    this->maxDepth = 5;
    this->directLighting = true;
    this->parseFile(filename);
}

Scene::~Scene()
{
    //delete mommaBox;
}

int Scene::renderImage()
{
    time_t begin = time(0);
    //createGrid();
    raytrace();
    time_t end = time(0);
    cout << "Rendered '" << outputFilename << "' in " << (end - begin) << "s.\n";
    return 0;
}

void Scene::raytrace()
{
    camera.setWidthAndHeight(screen.getWidth(), screen.getHeight());
    srand ( time(NULL) );
    int numberOfPixels = screen.getWidth() * screen.getHeight();
    double spp = (double) samplesPerPixel;
    if (traceType == RAY) {
        cout << "Raytracing\n";
    } else if (traceType == PATH) {
        cout << "Pathtracing\n";
    }
    cout << "Samples per pixel: " << samplesPerPixel << "\n";

    cout << "Indirect lighting: ";
    if (indirectLighting == true) {
        cout << "on\n";
    } else {
        cout << "off\n";
    }

    cout << "Direct lighting: ";
    if (directLighting == true) {
        cout << "on\n";
    } else {
        cout << "off\n";
    }


    Raytracer tracer = Raytracer(this);

    #pragma omp parallel for
    for (int i = 0; i < numberOfPixels; i++) {
        Ray ray;
        dvec3 color;
        dvec3 totalColor = dvec3(0.0f);
        int x = i % screen.getWidth();
        int y = i / screen.getWidth();
        Sample sample;
        Sampler sampler = Sampler(x, y, samplesPerPixel, samplerType);
        while (sampler.hasSamples() != false) {
            sample = sampler.getSample();
            camera.generateRay(ray,sample.x,sample.y);
            color = tracer.pathTraceRay(ray, 0, 1.0f, 0);
            totalColor += color;
        }
        totalColor /= spp;
        screen.writePixel(totalColor,x,y);
    }
    screen.saveScreenshot(("testscenes/" + outputFilename).c_str());
}

void Scene::parseFile(const char* filename)
{
    string str, cmd; 
    ifstream in;
    in.open(filename); 
    if (in.is_open()) {
        getline(in, str);

        // Transformation temporary
        transfstack.push(dmat4(1.0));  // identity

        // Geometry temporaries
        vector<dvec3> tempVertices;
        vector<dvec3> tempVerticesN;
        vector<dvec3> tempNormals;

        Material tempMaterial;
        tempMaterial.specular  = dvec3(0.0,0.0,0.0);
        tempMaterial.diffuse   = dvec3(0.0,0.0,0.0);
        tempMaterial.emission  = dvec3(0.0,0.0,0.0);
        tempMaterial.shininess = 0.0;
        tempMaterial.alpha     = 1.0;
        tempMaterial.rindex    = 1.0;
        tempMaterial.type      = LAMBERTIAN;   

        while (in) {
            if ((str.find_first_not_of(" \t\r\n") != string::npos) && (str[0] != '#')) {
                stringstream s(str);
                s >> cmd; 
                int i; 
                double values[15];
                bool validinput;

                // Settings
                if (cmd == "size") {
                    validinput = readvals(s,2,values);
                    if (validinput) {
                        screen.init(values[0], values[1]);
                    }
                } else if (cmd == "maxdepth") {
                    validinput = readvals(s,1,values);
                    if (validinput) {
                        this->maxDepth = values[0];
                    }
                } else if (cmd == "camera") {
                    validinput = readvals(s,10,values);
                    if (validinput) {
                        camera.init(values);
                    }
                } else if (cmd == "output") {
                    s >> outputFilename;
                } else if (cmd == "samplesperpixel") {
                    validinput = readvals(s,1,values);
                    if (validinput) {
                        samplesPerPixel = values[0];
                    }
                } else if (cmd == "tracingtype") {
                    string traceTypes;
                    s >> traceTypes;

                    if (traceTypes.compare("ray") == 0) {
                        this->traceType = RAY;
                    } else if (traceTypes.compare("path") == 0) {
                        this->traceType = PATH;
                    } else {
                        cout << "Unrecognizable tracing type. Defaulting to ray tracing.\n";
                        this->traceType = RAY;
                    }
                } else if (cmd == "samplesperpixel") {
                    validinput = readvals(s,1,values);
                    if (validinput) {
                        samplesPerPixel = values[0];
                    }
                } else if (cmd == "sampling") {
                    string samplingType;
                    s >> samplingType;

                    if (samplingType.compare("regular") == 0) {
                        samplerType = REGULAR;
                    } else if (samplingType.compare("jittered") == 0) {
                        samplerType = JITTERED;
                    } else if (samplingType.compare("random") == 0) {
                        samplerType = RANDOM;
                    } else {
                        cout << "Unrecognizable sampling type. Defaulting to regular\n";
                    }
                } else if (cmd == "direct") {
                    string lightingType;
                    s >> lightingType;

                    if (lightingType.compare("on") == 0) {
                        directLighting = true;
                    } else if (lightingType.compare("off") == 0) {
                        directLighting = false;
                    }
                } else if (cmd == "indirect") {
                    string lightingType;
                    s >> lightingType;

                    if (lightingType.compare("on") == 0) {
                        indirectLighting = true;
                    } else if (lightingType.compare("off") == 0) {
                        indirectLighting = false;
                    }
                }

                // Geometry
                else if (cmd == "sphere") {
                    validinput = readvals(s,4,values);
                    if (validinput) {
                        Sphere* sphere = new Sphere(values[0],values[1],values[2],values[3]);
                        Primitive prim(sphere);
                        prim.setAmbient(ambient);
                        prim.setMaterial(tempMaterial);
                        prim.setTransformation(transfstack.top());
                        primitiveList.push_back(prim);
                    }
                } else if (cmd == "vertex") {
                    validinput = readvals(s,3,values);
                    if (validinput) {
                        dvec3 vert = dvec3(values[0], values[1], values[2]);
                        tempVertices.push_back(vert);
                    }
                } else if (cmd == "vertexnormal") {
                    validinput = readvals(s,6,values);
                    if (validinput) {
                        dvec3 vert = dvec3(values[0], values[1], values[2]);
                        dvec3 norm = dvec3(values[3], values[4], values[5]);
                        tempVerticesN.push_back(vert);
                        tempNormals.push_back(norm);
                    }
                } else if (cmd == "tri") {
                    validinput = readvals(s,3,values);
                    if (validinput) {
                        Triangle* triangle = new Triangle(tempVertices[values[0]],tempVertices[values[1]], tempVertices[values[2]]);
                        Primitive prim(triangle);
                        prim.setAmbient(ambient);
                        prim.setMaterial(tempMaterial);
                        prim.setTransformation(transfstack.top());
                        primitiveList.push_back(prim);
                    }
                } else if (cmd == "trinormal") {
                    validinput = readvals(s,3,values);
                    if (validinput) {
                        Triangle* triangle = new Triangle(tempVerticesN[values[0]],tempVerticesN[values[1]], tempVerticesN[values[2]],
                                                          tempNormals[values[0]], tempNormals[values[1]], tempNormals[values[2]]);
                        Primitive prim(triangle);
                        prim.setAmbient(ambient);
                        prim.setMaterial(tempMaterial);
                        prim.setTransformation(transfstack.top());
                        primitiveList.push_back(prim);
                    }
                } 

                // Transformations
                else if (cmd == "translate") {
                    validinput = readvals(s,3,values); 
                    if (validinput) {
                        dmat4 mat = Transform::translate(values[0],values[1],values[2]);
                        rightMultiply(mat,transfstack);
                    }
                } else if (cmd == "scale") {
                    validinput = readvals(s,3,values); 
                    if (validinput) {
                        dmat4 mat = Transform::scale(values[0],values[1],values[2]);
                        rightMultiply(mat,transfstack);
                    }
                } else if (cmd == "rotate") {
                    validinput = readvals(s,4,values);
                    if (validinput) { 
                        dmat4 mat = Transform::rotate(values[3], dvec3(values[0],values[1],values[2]));
                        rightMultiply(mat,transfstack);
                    }
                } else if (cmd == "pushTransform") {
                    transfstack.push(transfstack.top());
                } else if (cmd == "popTransform") {
                    if (transfstack.size() <= 1) 
                        cerr << "  Stack has no elements.  Cannot Pop\n"; 
                    else {
                        transfstack.pop(); 
                    }
                }

                // Lights
                else if (cmd == "directional") {
                    validinput = readvals(s,6,values);
                    if (validinput) {
                        Light light;
                        light.posdir = dvec4(values[0],values[1],values[2],0);
                        light.color = dvec3(values[3],values[4],values[5]);
                        light.type = DIRECTIONAL;
                        lightList.push_back(light);
                    }
                } else if (cmd == "point") {
                    validinput = readvals(s,6,values);
                    if (validinput) {
                        Light light;
                        light.posdir = dvec4(values[0],values[1],values[2],1);
                        light.color = dvec3(values[3],values[4],values[5]);
                        light.type = POINT;
                        lightList.push_back(light);
                    }
                } else if (cmd == "area") { // This is different from the other light sources
                    string lightType;       // since it will actually create geometry in space
                    s >> lightType;
                    if (lightType.compare("rect") == 0) {
                        validinput = readvals(s,13,values);
                        if (validinput) {
                            /*Light light;
                            light.posdir    = dvec4(values[0],values[1],values[2],1);
                            light.upStep    = dvec4(values[3]/ (double) values[12],values[4]/ (double) values[12],values[5]/ (double) values[12],1.0);
                            light.rightStep = dvec4(values[6]/ (double) values[12],values[7]/ (double) values[12],values[8]/ (double) values[12],1.0);
                            light.color     = dvec3(values[9],values[10],values[11]);

                            light.type = AREA;
                            light.areaType = RECT;
                            light.numSamples = values[12]; // Will actually take square of this!!
                            lightList.push_back(light);*/

                            // Create the light to be used in lighting calculations
                            QuadLight* light = new QuadLight(dvec3(values[0],values[1],values[2]), 
                                                             dvec3(values[3],values[4],values[5]),
                                                             dvec3(values[6],values[7],values[8]),
                                                             dvec3(values[9],values[10],values[11]));
                            areaLightList.push_back(light);

                            // Create the geometry representing the light in the scene
                            dvec3 vert0 = dvec3(values[0],values[1],values[2]);
                            dvec3 vert1 = vert0 + dvec3(values[3],values[4],values[5]);
                            dvec3 vert2 = vert0 + dvec3(values[6],values[7],values[8]);
                            dvec3 vert3 = vert1 + vert2 - vert0;

                            Triangle* triangle = new Triangle(vert0, vert2, vert1);
                            Primitive prim(triangle);
                            prim.setAmbient(ambient);
                            dvec3 emis = dvec3(values[9],values[10],values[11]);
                            dvec3 other = dvec3(0.0);

                            Material matt;
                            matt.specular = other;
                            matt.diffuse = other;
                            matt.emission = emis;
                            matt.alpha = 0;
                            matt.shininess = 1;
                            matt.type = EMISSIVE;

                            prim.setMaterial(matt);
                            prim.setTransformation(transfstack.top());
                            primitiveList.push_back(prim);

                            Triangle* triangle2 = new Triangle(vert2, vert3, vert1);
                            Primitive prim2(triangle2);
                            prim2.setAmbient(ambient);
                            prim2.setMaterial(matt);
                            prim2.setTransformation(transfstack.top());
                            primitiveList.push_back(prim2);
                        }
                    } else {
                        cout << "Unrecognizable light type: " << lightType << ", skipping.";
                    }

                } else if (cmd == "attenuation") {
                    validinput = readvals(s,3,values);
                    if (validinput) {
                        this->attenuation[0] = values[0];
                        this->attenuation[1] = values[1];
                        this->attenuation[2] = values[2];
                    }
                } else if (cmd == "ambient") {
                    validinput = readvals(s,3,values);
                    if (validinput) {
                        ambient = dvec3(values[0], values[1], values[2]);
                    }
                } 

                // Materials
                else if (cmd == "diffuse") {
                    validinput = readvals(s,3,values);
                    if (validinput) {
                        tempMaterial.diffuse = dvec3(values[0], values[1], values[2]);
                    }
                } else if (cmd == "specular") {
                    validinput = readvals(s,3,values);
                    if (validinput) {
                        tempMaterial.specular = dvec3(values[0], values[1], values[2]);
                    }
                } else if (cmd == "emission") {
                    validinput = readvals(s,3,values);
                    if (validinput) {
                        tempMaterial.emission = dvec3(values[0], values[1], values[2]);
                    }
                } else if (cmd == "shininess") {
                    validinput = readvals(s,1,values);
                    if (validinput) {
                        tempMaterial.shininess = values[0];
                    }
                } else if (cmd == "alpha") {
                    validinput = readvals(s,1,values);
                    if (validinput) {
                        tempMaterial.alpha = values[0];
                    }
                } else if (cmd == "rindex") {
                    validinput = readvals(s,1,values);
                    if (validinput) {
                        tempMaterial.rindex = values[0];
                    }
                } else if (cmd == "material") {
                    string materialType;
                    s >> materialType;

                    if (materialType.compare("lambertian") == 0) {
                        tempMaterial.type = LAMBERTIAN;
                    } else if (materialType.compare("glossy") == 0) {
                        tempMaterial.type = GLOSSY;
                    } else if (materialType.compare("transmissive") == 0) {
                        tempMaterial.type = TRANSMISSIVE;
                    } else if (materialType.compare("reflective") == 0) {
                        tempMaterial.type = REFLECTIVE;
                    } else {
                        cout << "Unrecognizable material type. Defaulting to previous material type.\n";
                    }
                }

                // Don't need these two
                else if (cmd == "maxverts") { 
                } else if (cmd == "maxvertnorms") {
                } else {
                    cerr << "Unknown command: " << cmd << " -- skipping \n";
                }
            }
            getline(in, str);
        }
    } else {
        cerr << "Unable to open input data file " << filename << "\n"; 
        throw 2; 
    }
}

bool Scene::readvals(stringstream &s, const int numvals, double* values) 
{
    for (int i = 0 ; i < numvals ; i++) {
        s >> values[i]; 
        if (s.fail()) {
            cout << "Failed reading value " << i << " will skip\n"; 
            return false;
        }
    }
    return true; 
}

void Scene::rightMultiply(const dmat4& M, stack<dmat4> &transfstack) {
    dmat4 &T = transfstack.top(); 
    T = M * T; 
}

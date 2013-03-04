/*
	Scene.cpp
*/

#include <iostream>
#include <fstream>
#include <deque>
#include <time.h>

#include <glm/glm.hpp>
typedef glm::vec3 vec3; 
typedef glm::vec4 vec4; 
typedef glm::mat4 mat4; 

#include "Scene.h"
#include "Ray.h"
#include "Triangle.h"
#include "Sphere.h"
#include "Material.h"
#include "Transform.h"

using namespace std;

int Raytracer::maxDepth = 5;

Scene::Scene(const char* filename) : camera(), screen()
{
    this->outputFilename = "screenshot.png";
    this->attenuation = vec3(1.0,0.0,0.0);
    this->ambient = vec3(0.2,0.2,0.2);
    this->parseFile(filename);
    this->gridStart = vec3(0.0,0.0,0.0);
    this->gridEnd = vec3(10.0,10.0,10.0);
    this->gridSize = 8;
    this->gi = 8;
    this->gidepth = 3;
}

int Scene::renderImage()
{
    time_t begin = time(0);
    createGrid();
    raytrace();
    time_t end = time(0);
    cout << "Rendered '" << outputFilename << "' in " << (end - begin) << "s.\n";
    return 0;
}

void Scene::raytrace()
{
    camera.setWidthAndHeight(screen.getWidth(), screen.getHeight());
    int* prevRow = new int[screen.getWidth()];
    for (int j = 0; j < screen.getHeight(); j++) {
        int thisIndex = -1;
        int prevIndex = -1;
        Ray* ray = new Ray();
        for (int i = 0; i < screen.getWidth(); i++) {
            vec3 totalColor = vec3(0.0);
            vec3 color = vec3(0.0);
            camera.generateRay(ray,i,j);
            thisIndex = Raytracer::traceRay(this, ray, 0, color, 1.0);

            // If we encounter an edge, supersample by 16 and average
            // Make this condition true to make it pure blind supersampling
            if (i != 0 && (prevIndex != thisIndex) || (j != 0 && thisIndex != prevRow[i])) {
                //std::cout << thisIndex << " vs. " << prevRow[i] << " for " << i << "," << j << "\n";
                for (int y = 0; y < 4; y++) {
                    for (int x = 0; x < 4; x++) {
                        float div = 1.0f / (RAND_MAX+1.0f);
                        float randx = (rand()*div+rand()*div*div-0.5f)/4.0f;
                        float randy = (rand()*div+rand()*div*div-0.5f)/4.0f;
                        camera.generateRay(ray, (i-0.5) + (float)x/4.0f,(j-0.5) + (float)y/4.0f);
                        //camera.generateRay(ray, (i-0.5) + (float)x/4.0f+randx,(j-0.5) + (float)y/4.0f+randy);
                        color = vec3(0.0);
                        Raytracer::traceRay(this, ray, 0, color, 1.0);
                        totalColor += color;
                    }
                }
                totalColor /= 16.0f;
            } else {
                totalColor = color;
            }
            prevIndex = thisIndex;
            prevRow[i] = thisIndex;
            screen.writePixel(totalColor,i,j);
        }
    }
    screen.saveScreenshot(("testscenes/" + outputFilename).c_str());
}

// Credit to flipcode for keeping me sane
void Scene::createGrid()
{
    grid.resize(gridSize*gridSize*gridSize);

    vec3 dVec = (gridEnd - gridStart) / (float) gridSize;
    vec3 dVecInv = 1.0f / dVec;
    this->mommaBox = new AABB();
    this->mommaBox->minimum = gridStart;
    this->mommaBox->maximum = gridEnd;
    AABB* gridBox = new AABB();

    for (int i = 0; i < primitiveList.size(); i++) {
        
        int* startVals = new int[3];
        int* endVals   = new int[3];
        for (int k = 0; k < 3; k++) {
            startVals[k] = (int) ((primitiveList[i].getAABB()->minimum[k] - gridStart[k])*dVecInv[k]);
            startVals[k] = startVals[k] < 0 ? 0 : startVals[k];
            endVals[k]   = (int) ((primitiveList[i].getAABB()->maximum[k] - gridStart[k])*dVecInv[k] + 1);
            endVals[k]   = endVals[k] > gridSize - 1 ? gridSize - 1 : endVals[k];
        }

        for (int x = startVals[0]; x < endVals[0]; x++) {
            for (int y = startVals[1]; y < endVals[1]; y++) {
                for (int z = startVals[2]; z < endVals[2]; z++) {
                    int index = x + y * gridSize + z * gridSize * gridSize;
                    vec3 minimum = vec3(gridStart.x + x * dVec.x,
                                        gridStart.y + y * dVec.y,
                                        gridStart.z + z * dVec.z);
                    gridBox->minimum = minimum;
                    gridBox->maximum = minimum + dVec;

                    // Check bounding box of object against gridbox
                    // If intersects, then add the primitive index into grid index list
                    if (primitiveList[i].getAABB()->doesIntersectAABB(gridBox)) {
                        grid[index].push_back(i);
                    }
                }
            }
        }
    }
}

void Scene::parseFile(const char* filename)
{
    string str, cmd; 
    ifstream in;
    in.open(filename); 
    if (in.is_open()) {
        getline(in, str);

        // Transformation temporary
        transfstack.push(mat4(1.0));  // identity

        // Geometry temporaries
        vector<vec3> tempVertices;
        vector<vec3> tempVerticesN;
        vector<vec3> tempNormals;

        Material* tempMaterial = new Material;
        tempMaterial->specular = vec3(0.0,0.0,0.0);
        tempMaterial->diffuse  = vec3(0.0,0.0,0.0);
        tempMaterial->emission = vec3(0.0,0.0,0.0);
        tempMaterial->shininess = 0.0;
        tempMaterial->alpha = 1.0;
        tempMaterial->rindex = 1.0;

        while (in) {
            if ((str.find_first_not_of(" \t\r\n") != string::npos) && (str[0] != '#')) {
                stringstream s(str);
                s >> cmd; 
                int i; 
                float values[10];
                bool validinput;

                // Setup
                if (cmd == "size") {
                    validinput = readvals(s,2,values);
                    if (validinput) {
                        screen.init(values[0], values[1]);
                    }
                } else if (cmd == "maxdepth") {
                    validinput = readvals(s,1,values);
                    if (validinput) {
                        Raytracer::maxDepth = values[0];
                    }
                } else if (cmd == "camera") {
                    validinput = readvals(s,10,values);
                    if (validinput) {
                        camera.init(values);
                    }
                } else if (cmd == "output") {
                    s >> outputFilename;
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
                        vec3 vert = vec3(values[0], values[1], values[2]);
                        tempVertices.push_back(vert);
                    }
                } else if (cmd == "vertexnormal") {
                    validinput = readvals(s,6,values);
                    if (validinput) {
                        vec3 vert = vec3(values[0], values[1], values[2]);
                        vec3 norm = vec3(values[3], values[4], values[5]);
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
                        mat4 mat = Transform::translate(values[0],values[1],values[2]);
                        rightMultiply(mat,transfstack);
                    }
                } else if (cmd == "scale") {
                    validinput = readvals(s,3,values); 
                    if (validinput) {
                        mat4 mat = Transform::scale(values[0],values[1],values[2]);
                        rightMultiply(mat,transfstack);
                    }
                } else if (cmd == "rotate") {
                    validinput = readvals(s,4,values);
                    if (validinput) { 
                        mat4 mat = Transform::rotate(values[3], vec3(values[0],values[1],values[2]));
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
                        light.posdir = vec4(values[0],values[1],values[2],0);
                        light.color = vec3(values[3],values[4],values[5]);
                        light.point = false;
                        lightList.push_back(light);
                    }
                } else if (cmd == "point") {
                    validinput = readvals(s,6,values);
                    if (validinput) {
                        Light light;
                        light.posdir = vec4(values[0],values[1],values[2],1);
                        light.color = vec3(values[3],values[4],values[5]);
                        light.point = true;
                        lightList.push_back(light);
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
                        ambient = vec3(values[0], values[1], values[2]);
                    }
                } 

                // Materials
                else if (cmd == "diffuse") {
                    validinput = readvals(s,3,values);
                    if (validinput) {
                        tempMaterial->diffuse = vec3(values[0], values[1], values[2]);
                    }
                } else if (cmd == "specular") {
                    validinput = readvals(s,3,values);
                    if (validinput) {
                        tempMaterial->specular = vec3(values[0], values[1], values[2]);
                    }
                } else if (cmd == "emission") {
                    validinput = readvals(s,3,values);
                    if (validinput) {
                        tempMaterial->emission = vec3(values[0], values[1], values[2]);
                    }
                } else if (cmd == "shininess") {
                    validinput = readvals(s,1,values);
                    if (validinput) {
                        tempMaterial->shininess = values[0];
                    }
                } else if (cmd == "alpha") {
                    validinput = readvals(s,1,values);
                    if (validinput) {
                        tempMaterial->alpha = values[0];
                    }
                } else if (cmd == "rindex") {
                    validinput = readvals(s,1,values);
                    if (validinput) {
                        tempMaterial->rindex = values[0];
                    }
                } 

                else if (cmd == "gi") {
                    validinput = readvals(s,1,values);
                    if (validinput) {
                        this->gi = values[0];
                    }
                } else if (cmd == "gidepth") {
                    validinput = readvals(s,1,values);
                    if (validinput) {
                        this->gidepth = values[0];
                    }
                } 

                // Grid stuff
                else if (cmd == "gridstart") {
                    validinput = readvals(s,3,values);
                    if (validinput) {
                        this->gridStart = vec3(values[0],values[1],values[2]);
                    }
                }  else if (cmd == "gridend") {
                    validinput = readvals(s,3,values);
                    if (validinput) {
                        this->gridEnd = vec3(values[0],values[1],values[2]);
                    }
                }  else if (cmd == "gridsize") {
                    validinput = readvals(s,1,values);
                    if (validinput) {
                        this->gridSize = values[0];
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

bool Scene::readvals(stringstream &s, const int numvals, float* values) 
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

void Scene::rightMultiply(const mat4 & M, stack<mat4> &transfstack) {
    mat4 &T = transfstack.top(); 
    T = M * T; 
}

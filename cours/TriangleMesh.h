#pragma once

#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include "vector.h"
#include "ray.h"
#include "Object.h"
#include "BoundingBox.h"

#include <list>


#include "stb_image.h"

#define M_PI 3.14159265358979323846


class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};


class TriangleMesh: public Object {

private:
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    std::vector<unsigned char*> textures;
    std::vector<int> Wtex, Htex;
    BoundingBox bb;
    Noeud* BVH;
    Vector up, fwd, right;


    BoundingBox buildBB(int debut, int fin);
    void readOBJ(const char* obj);
    void buildBVH(Noeud* n, int debut, int fin);
    void loadTexture(const char* filename);

    /** 
        Les méthodes translate, inverse, rotate et resize ont des noms explicites, translate relative permet de déplacer l'objet par rapport 
        à son référentiel propre qui a été définit au moment de construire l'objet la première fois. 

        A noter que vu la façon dont son chargés les objets, toutes ces opérations se font par rapport au référentiel de la scene. On pourrait 
        ajouter un attribut position qui définirait la position de l'objet dans la scene et l'ensemble des points seraient alors définit par rapport 
        à cette position de référence.
    */
    void translate(float tx, float ty, float tz);
    void translateRelative(float tx, float ty, float tz);
    void inverse(bool x, bool y, bool z);
    void rotate(float rightRotation, float upRotation, float fwRotation);
    void resize(float k, float l, float m);



public:
   
    ~TriangleMesh() {}
    TriangleMesh(const char* obj, const char* texture, float tx = 0, float ty = 0, float tz = 0, float rx = 0, float ry = 0, float rz = 0, bool mirrorx = false, bool mirrory = false, bool mirrorz = false, float rFactorX = 1, float rFactorY = 1, float rFactorZ = 1) {
        this->albedo = Vector(0,0,0);
        this->isMirror = false;
        this->isTransparent = false;
        readOBJ(obj);
        loadTexture(texture);
        // fwd, right et up permettent de  définir le référentiel propre à l'objet
        fwd = Vector(0, 0, 1);
        right =  Vector(1, 0, 0);
        up = cross(right, fwd);
        transformMesh(tx, ty, tz, rx, ry, rz, mirrorx, mirrory, mirrorz, rFactorX, rFactorY, rFactorZ);
    };

    
    bool intersect(const Ray& r, Vector& P, Vector& normale, double& t, Vector& color);

    /** 
        translateMesh, rotateMesh, inverseMesh, resizeMesh et translateRelativeMesh sont les versions publiques des méthodes définies plus haut. 
        Elles forcent la reconstruction de la BVH à chaque utilisation.

        TransformMesh est une agrégation des méthodes privées translate, inverse, rotate et resize puis reconstruit la BVH.
    */
    void transformMesh(float tx = 0, float ty = 0, float tz = 0, float rx = 0, float ry = 0, float rz = 0, bool mirrorx = false, bool mirrory = false, bool mirrorz = false, float rFactorX = 1, float rFactorY = 1, float rFactorZ = 1);
    void translateMesh(float tx, float ty, float tz);
    void translateRelativeMesh(float tx, float ty, float tz);
    void inverseMesh(bool x, bool y, bool z);
    void rotateMesh(float rightRotation, float upRotation, float fwRotation);
    void resizeMesh(float k, float l, float m);
};



#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
//#include "stb_image.h"

#include "vector.h"
#include "sphere.h"
#include "ray.h"
#include "scene.h"
#include "TriangleMesh.h"

#include <iostream>
#include <sstream>

#define M_PI 3.14159265358979323846

#include <time.h>





int main() {

    clock_t t_debut = clock();
    clock_t t;
    clock_t t_restant;

    /****************************************/
    /*                                      */
    /*      Paramétrisation de l'image      */
    /*                                      */
    /****************************************/
    int W = 1024;
    int H = 1024;
    int nbRays = 300;


    double tailleLentille = 0.000001;
    double focale = 55;

    /******************************************/
    /*                                        */
    /*      Paramétrisation de la caméra      */
    /*                                        */
    /******************************************/
    
    // Position
    float xCam = 20;
    float yCam = 50;
    float zCam = 55;

    // Rotation (en degrés)
    float phiCam = 20; // droite, gauche
    float thetaCam = -35; // haut, bas
    float psiCam = 0; // "torsion"

    // FOV (en degrés)
    double alphaCam = 60;





    Vector C(xCam, yCam, zCam);

    /******************************************/
    /*                                        */
    /*      Paramétrisation de la scene       */
    /*                                        */
    /******************************************/

    float xLight = -10;
    float yLight = 50;
    float zLight = 40;

    Vector lightPos = Vector(xLight, yLight, zLight);

    float IntensiteLight = 5E9;






    Scene scene;
    scene.L = lightPos;
    scene.I = IntensiteLight;

    /******************************************/
    /*                                        */
    /*           Insertion des objets         */
    /*                                        */
    /******************************************/


    /*********** Lumière ************/

    Sphere Source(scene.L, 10, Vector(1, 1, 1), false, false);
    scene.objects.push_back(&Source);




    /*********** Sphères ************/

    Sphere S1(Vector(0, 10, 15), 5, Vector(1, 1, 1), true, false);
    scene.objects.push_back(&S1);

    Sphere S2(Vector(-10, 0, -20), 10, Vector(1, 1, 1), true, false);
    //scene.objects.push_back(&S2);

    Sphere S3(Vector(10, 0, 20), 10, Vector(1, 1, 1), false, true);
    //scene.objects.push_back(&S3);

    Sphere Ssol(Vector(0, -1000, 0), 990, Vector(0., 0., 1));
    scene.objects.push_back(&Ssol);

    Sphere Splafond(Vector(0, 1000, 0), 940, Vector(1, 0., 0.));
    scene.objects.push_back(&Splafond);

    Sphere Smur1(Vector(-1000, 0, 0), 960, Vector(1, 0, 0));
    scene.objects.push_back(&Smur1);

    Sphere Smur2(Vector(1000, 0, 0), 960, Vector(1, 1, 0));
    scene.objects.push_back(&Smur2);

    Sphere Smur3(Vector(0, 0, -1000), 940, Vector(0.18, 0.65, 0.18));
    scene.objects.push_back(&Smur3);

    Sphere Smur4(Vector(0, 0, 1000), 940, Vector(1, 1, 0.5));
    scene.objects.push_back(&Smur4);




    /*********** Meshes ************/
    /*
    TriangleMesh :
    const char* obj;
    const char* texture;

    float translateX = 0;
    float translateY = 0;
    float translateZ = 0;

    // en degré, sens trigo
    float rotationRight = 0;
    float rotationUp = 0;
    float rotationFwd = 0;

    bool mirrorx = false;
    bool mirrory = false;
    bool mirrorz = false;

    float resizeX = 1;
    float resizeY = 1;
    float resizeZ = 1;
    */


    TriangleMesh m("Rayman3model.obj", "rayman.PNG", -10, 0, 0, 0, 0, 0, false, false, false, 1, 1, 1);
    scene.objects.push_back(&m);

    TriangleMesh m1("Rayman3model.obj", "rayman.PNG", 15, 0, 0, 0, 30, 0, false, false, false, 0.5, 0.5, 0.5);
    scene.objects.push_back(&m1);

    TriangleMesh m2("Rayman3model.obj", "rayman.PNG", 0, 0, 30, 0, 180, 0, false, false, false, 0.5, 0.5, 0.5);
    scene.objects.push_back(&m2);

    






    thetaCam *= M_PI / 180;
    phiCam *= M_PI / 180;
    psiCam *= M_PI / 180;
    alphaCam *= M_PI / 180;


    Vector up(0, 1, 0);
    Vector right(1, 0, 0);
    Vector viewDir = cross(up, right);


    right = right.getRotatedVector(up, phiCam);
    viewDir = viewDir.getRotatedVector(up, phiCam);

    up = up.getRotatedVector(right, thetaCam);
    viewDir = viewDir.getRotatedVector(right, thetaCam);

    up = up.getRotatedVector(viewDir, psiCam);
    right = right.getRotatedVector(viewDir, psiCam);

    std::vector<unsigned char> image(W * H * 3, 0);

#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        std::cout << "Pixel hauteur: " << i << std::endl;
        for (int j = 0; j < W; j++) {

            Vector color(0, 0, 0);

            for (int k = 0; k < nbRays; k++) {
                double u1 = uniform(engine);
                double u2 = uniform(engine);
                double dx = 0.25 * cos(2 * M_PI * u1) * sqrt(-2 * log(u2));
                double dy = 0.25 * sin(2 * M_PI * u1) * sqrt(-2 * log(u2));


                double u3 = uniform(engine);
                double u4 = uniform(engine);
                double dx1 = tailleLentille * cos(2 * M_PI * u3) * sqrt(-2 * log(u4));
                double dy1 = tailleLentille * sin(2 * M_PI * u3) * sqrt(-2 * log(u4));

                Vector direction(j - W / 2 + 0.5 + dx, i - H / 2 + 0.5 + dy, W / (2 * tan(alphaCam / 2)));
                direction = direction.getNormalized();
                direction = direction[0] * right + direction[1] * up + direction[2] * viewDir;

                Vector target = C + focale * direction;
                Vector Cprime = C + dx1 * right + dy1 * up;
                Vector directionprime = (target - Cprime).getNormalized();

                Ray r(Cprime, directionprime);

                color = color + scene.getColor(r, 0, false);
            }

            color = color / nbRays;

            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));

        }
        if (i % 100 == 99) {
            t = (clock() - t_debut) / 1000.;
            std::cout << "temps écoulé depuis le début : " << (int)t / 60 << "m" << t % 60 << "s" << std::endl;
            t_restant = t / (i + 1) * (H - i - 1);
            std::cout << "temps restant estimé : " << (int)t_restant / 60 << "m" << t_restant % 60 << "s" << std::endl;
            // Le temps restant estimé ne dépends que de ce qu'il y a au dessus dans l'image donc une image très chargée seulement en haut ou en bas donnera un resultat erroné
            // On pourrait améliorer cela en ne faisant pas les lignes dans l'ordre mais en parcourant une permutation mais ce ne sera pas fait ici
        }
    }

    t = (clock() - t_debut)/1000.; 
    int min = (int)t / 60;

    std::cout << "Entrez le nom de l'image (commencer par _ pour 'test.png') : ";
    char* nomImage = new char;
    std::cin >> nomImage;
   
    if (*nomImage != '_') {
        std::ostringstream out;
        out << nomImage << "-" << W << "x" << H << "_" << nbRays << "Rayons_" << min << "m" << t % 60 << "s.png";
        std::string location = out.str();
        char* cname;
        cname = &location[0];
        stbi_write_png(cname, W, H, 3, &image[0], 0);
    }
    else {
        stbi_write_png("test.png", W, H, 3, &image[0], 0);
    }




    return 0;
}


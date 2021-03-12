#include "scene.h"

/**
    Renvoie un vecteur avec une direction aléatoire proche de la direction de N
*/
Vector random_cos(const Vector& N) {
    double u1 = uniform(engine);
    double u2 = uniform(engine);
    double x = sin(2 * M_PI * u1) * sqrt(1 - u2);
    double y = cos(2 * M_PI * u1) * sqrt(1 - u2);
    double z = sqrt(u2);
    Vector T1;
    // On élimine la composante la plus petite
    if (N[0] < N[1] && N[0] < N[2])
    {
        T1 = Vector(0, -N[2], N[1]);
    }
    else if (N[1] < N[0] && N[1] < N[2])
    {
        T1 = Vector(-N[2], 0, N[0]);
    }
    else
    {
        T1 = Vector(-N[1], N[0], 0);
    }
    T1 = T1.getNormalized();
    Vector T2 = cross(N, T1);
    return (x * T1 + y * T2 + z * N);
}


/*
    Prend en entrée un rayon, trois vecteurs, deux booleens, un double et un entier.
    Renvoie true s'il a une intersection entre un des objets de la scene et un rayon et que le point d'intersection est 
    plus proche du point de départ du rayon que l'intersection la plus proche trouvée précedemment. 
    Si c'est le cas, les entrées sont changées pour stocker le nouveau point d'intersection, la normale à ce point,
    l'albedo de ce point, si l'objet intersecté à des propriété de reflexion ou de transparence, la nouvelle distance entre l'intersection et 
    le point de départ du rayon et l'indice de l'objet intersecté.
*/
bool Scene::intersect(const Ray& r, Vector& P, Vector& N, Vector& albedo, bool& mirror, bool& transparent, double& t, int& objectid ) {
    t = 1E10; //infini
    bool has_inter = false;
    for (int i = 0; i < objects.size(); i++) {
        Vector localP, localN, localAlbedo;
        double localt;
        if (objects[i]->intersect(r, localP, localN, localt, localAlbedo) && localt < t) {
            t = localt;
            has_inter = true;
            // albedo = objects[i]->albedo;
            albedo = localAlbedo;
            P = localP;
            N = localN;
            mirror = objects[i]->isMirror;
            transparent = objects[i]->isTransparent;
            objectid = i;
        }
    }
    return has_inter;
}

/*
    Prend en entrée un rayon et renvoie un vecteur contenaant la couleur vue par ce rayon le paramètre rebond permet de gérer la récursion
    lastdefuse permet de renvoyer la couleur de la source de lumière si le dernier rebond ne s'est pas fait sur une surface diffuse. Ainsi
    elle n'apparaitre pas noire dans un mirroir.
*/

Vector Scene::getColor(const Ray& r, int rebond, bool lastdefuse) {

    double eps = 1E-5; // Variable pour gérer les incerrtitudes dans les floats et empecher les autointersections

    if (rebond > 5)
        return Vector(0, 0, 0);

    Vector P, N, albedo;
    double t;
    bool mirror, transparent;
    int objectid;
    bool inter = intersect(r, P, N, albedo, mirror, transparent, t, objectid);
    Vector color(0, 0, 0);

    if (inter) {
        if (objectid == 0) { // L'objet avec l'idée 0 représente la lumière on considère que ça ne peut qu'être une sphere
            if (rebond == 0 || !lastdefuse) {
                return Vector(I, I, I) / (4 * M_PI * M_PI * sqr(dynamic_cast<Sphere*>(objects[0])->getR()));
            }
            else {
                return Vector(0, 0, 0);
            }
        }


        if (mirror) {
            Vector reflectedDir = r.u - 2 * dot(r.u, N) * N;
            Ray reflectedRay(P + eps * N, reflectedDir);
            return getColor(reflectedRay, rebond + 1, false);
        }
        else if (transparent) {
            double n1 = 1, n2 = 1.4;
            Vector N2 = N;
            if (dot(r.u, N) > 0) { // On sort de la sphere
                std::swap(n1, n2);
                N2 = -N;
            }
            Vector Tt = n1 / n2 * (r.u - dot(r.u, N2) * N2);
            double rad = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.u, N2)));
            if (rad < 0) {
                Vector reflectedDir = r.u - 2 * dot(r.u, N) * N;
                Ray reflectedRay(P + eps * N, reflectedDir);
                return getColor(reflectedRay, rebond + 1, false);
            }
            Vector Tn = -sqrt(rad) * N2;

            Vector refractedDir = Tt + Tn;
            Ray refractedRay(P - eps * N2, refractedDir);
            return getColor(refractedRay, rebond, false);
        }
        else { // ni transparent, ni miroir
            /*
            //eclairage direct
            Vector PL = L - P;
            double d = sqrt(PL.sqrNorm());
            Vector shadowP, shadowN, shadowAlbedo;
            double shadowt;
            Ray shadowRay(P + 0.001 * N, PL / d);
            bool shadowMirror, shadowTransparent;
            bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransparent, shadowt);

            if (shadowInter && shadowt < d) {
                color = Vector(0., 0., 0.);
            }
            else {
                color = I / (4 * M_PI * d * d) * albedo / M_PI * std::max(0., dot(N, PL / d));
            }
            */

            // eclairage direct
            Vector PL = L - P;
            PL = PL.getNormalized();
            Vector w = random_cos(-PL);
            double rayon = dynamic_cast<Sphere*>(objects[0])->getR();
            Vector origine = dynamic_cast<Sphere*>(objects[0])->getO();
            Vector xprime = w * rayon + origine;
            Vector Pxprime = xprime - P;
            double d = sqrt(Pxprime.sqrNorm());
            Pxprime = Pxprime / d;


            Vector shadowP, shadowN, shadowAlbedo;
            double shadowt;
            bool shadowMirror, shadowTransparent;
            int objectid;
            Ray shadowRay(P + eps * N, Pxprime);
            bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransparent, shadowt, objectid);

            if (shadowInter && shadowt < d-3*eps) {
                color = Vector(0., 0., 0.);
            }
            else {
                
                double proba = dot(-PL, w)/(M_PI * rayon * rayon);
                double J = dot(w, -Pxprime)/(d*d);
                color = I / (4 * M_PI * M_PI * sqr(rayon)) * albedo / M_PI * std::max(0., dot(N, Pxprime)) * J / proba;
            }


            // eclairage indirect
            Vector wi = random_cos(N);
            Ray wiRay(P + eps * N, wi.getNormalized());
            color = color + albedo / M_PI * getColor(wiRay, rebond + 1, true);


        }
    }
    return color;
}



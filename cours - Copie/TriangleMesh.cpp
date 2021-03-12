#include "TriangleMesh.h"


void TriangleMesh::readOBJ(const char* obj) {

    char matfile[255];
    char grp[255];

    FILE* f;
    f = fopen(obj, "r");
    int curGroup = -1;
    while (!feof(f)) {
        char line[255];
        if (!fgets(line, 255, f)) break;

        std::string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());

        if (line[0] == 'u' && line[1] == 's') {
            sscanf(line, "usemtl %[^\n]\n", grp);
            curGroup++;
        }

        if (line[0] == 'v' && line[1] == ' ') {
            Vector vec;

            Vector col;
            if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                col[0] = std::min(1., std::max(0., col[0]));
                col[1] = std::min(1., std::max(0., col[1]));
                col[2] = std::min(1., std::max(0., col[2]));

                vertices.push_back(vec);
                vertexcolors.push_back(col);

            }
            else {
                sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                vertices.push_back(vec);
            }
        }
        if (line[0] == 'v' && line[1] == 'n') {
            Vector vec;
            sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
            normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't') {
            Vector vec;
            sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
            uvs.push_back(vec);
        }
        if (line[0] == 'f') {
            TriangleIndices t;
            int i0, i1, i2, i3;
            int j0, j1, j2, j3;
            int k0, k1, k2, k3;
            int nn;
            t.group = curGroup;

            char* consumedline = line + 1;
            int offset;

            nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
            if (nn == 9) {
                if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                indices.push_back(t);
            }
            else {
                nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                if (nn == 6) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    indices.push_back(t);
                }
                else {
                    nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                    if (nn == 3) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        indices.push_back(t);
                    }
                    else {
                        nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                        if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                        if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                        indices.push_back(t);
                    }
                }
            }

            consumedline = consumedline + offset;

            while (true) {
                if (consumedline[0] == '\n') break;
                if (consumedline[0] == '\0') break;
                nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                TriangleIndices t2;
                t2.group = curGroup;
                if (nn == 3) {
                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                    if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                    if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                    if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                    if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                    if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                    if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                    indices.push_back(t2);
                    consumedline = consumedline + offset;
                    i2 = i3;
                    j2 = j3;
                    k2 = k3;
                }
                else {
                    nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                    if (nn == 2) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        indices.push_back(t2);
                    }
                    else {
                        nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                            if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                            if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            k2 = k3;
                            indices.push_back(t2);
                        }
                        else {
                            nn = sscanf(consumedline, "%u%n", &i3, &offset);
                            if (nn == 1) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                indices.push_back(t2);
                            }
                            else {
                                consumedline = consumedline + 1;
                            }
                        }
                    }
                }
            }

        }

    }
    fclose(f);

}

/*
Prend en entrée un rayon, trois vecteurs et un double.
Renvoie true s'il a une intersection entre le mesh et le rayon.
Si c'est le cas, les entrées sont changées pour stocker le point d'intersection, la normale à ce point,
l'albedo de ce point ainsi que la distance e ce point par rapport au point de départ du rayon
*/

bool TriangleMesh::intersect(const Ray& r, Vector& P, Vector& normale, double& t, Vector& color) {
    
    double localt;
    if (!BVH->b.intersect(r , localt))
        return false;

    t = 1E9;
    bool has_inter = false;

    std::list<Noeud*> l;
    l.push_back(BVH);

    while (!l.empty()) {

        Noeud* c = l.front();
        l.pop_front();


        if (c->fg) {
            
            double localtfg;
            double localtfd; 

            bool interfg = c->fg->b.intersect(r, localtfg);
            bool interfd = c->fd->b.intersect(r, localtfd);

            if (interfg && interfd && localtfd < t && localtfg < t) {
                if (localtfg < localtfd) {
                    l.push_front(c->fd);
                    l.push_front(c->fg);
                }
                else {
                    l.push_front(c->fg);
                    l.push_front(c->fd);
                }
            }
            else {
                if (interfg && localtfg < t) {
                    l.push_front(c->fg);
                }
                if (interfd && localtfd < t) {
                    l.push_front(c->fd);
                }
            }

        } else { // C'est une feuille

            for (int i = c->debut; i < c->fin; i++) {
                const Vector& A = vertices[indices[i].vtxi];
                const Vector& B = vertices[indices[i].vtxj];
                const Vector& C = vertices[indices[i].vtxk];

                Vector e1 = B - A;
                Vector e2 = C - A;
                Vector N = cross(e1, e2);
                Vector AO = r.C - A;
                Vector AOu = cross(AO, r.u);

                double invUN = 1. / dot(r.u, N);
                double beta = -dot(e2, AOu) * invUN;
                double gamma = dot(e1, AOu) * invUN;
                double alpha = 1 - beta - gamma;
                localt = -dot(AO, N) * invUN;

                if (beta >= 0 && gamma >= 0 && alpha >= 0 && beta <= 1 && gamma <= 1 && localt > 0) {
                    has_inter = true;
                    if (localt < t) {
                        t = localt;
                        normale = alpha * normals[indices[i].ni] + beta * normals[indices[i].nj] + gamma * normals[indices[i].nk];
                        normale = normale.getNormalized();
                        P = r.C + t * r.u;

                        Vector UV = alpha * uvs[indices[i].uvi] + beta * uvs[indices[i].uvj] + gamma * uvs[indices[i].uvk];
                        int H = Htex[indices[i].group];
                        int W = Wtex[indices[i].group];

                        UV = UV * Vector(W, H, 0);
                        int uvx = UV[0] + 0.5;
                        int uvy = UV[1] + 0.5;

                        uvx = uvx % W;
                        uvy = uvy % H;

                        if(uvx < 0)
                            uvx += W;

                        if (uvy < 0)
                            uvy += H;

                        uvy = H - uvy - 1; // L'origine sur la texture se trouve en bas

                        // puissance 2.2 car le puissance 1/2.2 est deja pris en compte dans la texture
                        color = Vector(
                            std::pow(textures[indices[i].group][(uvy * W + uvx) * 3] / 255., 2.2),
                            std::pow(textures[indices[i].group][(uvy * W + uvx) * 3 + 1] / 255., 2.2),
                            std::pow(textures[indices[i].group][(uvy * W + uvx) * 3 + 2] / 255., 2.2)
                            );  
                    }
                }

            }



        }


    }

    return has_inter;
}


/*
    Construit une boite autour des sommets d'indice debut et fin
*/
BoundingBox TriangleMesh::buildBB(int debut, int fin) {
    BoundingBox bb;

    bb.mini = Vector(1E9, 1E9, 1E9);
    bb.maxi = Vector(-1E9, -1E9, -1E9);

    for (int i = debut; i < fin; i++) {
        for (int j = 0; j < 3; j++) {
            bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxi][j]);
            bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxi][j]);

            bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxj][j]);
            bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxj][j]);

            bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxk][j]);
            bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxk][j]);
        }
    }

    return bb;
}


/*
    Construit une BVH autour des sommets d'indices debut et fin :  Si cela comprend peu de sommet, construit directement une bounding box, 
    sinon sépare en 2 les sommets dans la dimension la plus grande et fais un appel récursif en "triant" au passage avec un quick sort les sommets selon cette dimension
*/
void TriangleMesh::buildBVH(Noeud* n, int debut, int fin) {

    n->debut = debut;
    n->fin = fin;

    n->b = buildBB(n->debut, n->fin);

    Vector diag = n->b.maxi - n->b.mini;

    int dim; 

    if (diag[0] >= diag[1] && diag[0] >= diag[2]) {
        dim = 0;
    }
    else {
        if (diag[1] >= diag[0] && diag[1] >= diag[0])
            dim = 1;
        else
            dim = 2;
    }

    double milieu = (n->b.mini[dim] + n->b.maxi[dim]) / 2;

    int iPivot = n->debut;

    for (int i = n->debut; i < n->fin; i++) {
        double milieuTriangle = (vertices[indices[i].vtxi][dim] + vertices[indices[i].vtxj][dim] + vertices[indices[i].vtxk][dim]) / 3.;
        // double milieuTriangle = vertices[indices[i].vtxi][dim];
        if (milieuTriangle < milieu) {
            std::swap(indices[i], indices[iPivot]);
            iPivot++;
        }
    }
    n->fg = nullptr;
    n->fd = nullptr;

    if (iPivot == debut || iPivot == fin || (fin - debut < 5)) return;

    n->fg = new Noeud;
    n->fd = new Noeud;

    buildBVH(n->fg, n->debut, iPivot);
    buildBVH(n->fd, iPivot, n->fin);
}



void TriangleMesh::loadTexture(const char* filename) {
    int W, H, C;
    unsigned char* texture = stbi_load(filename, &W, &H, &C, 3);
    Wtex.push_back(W);
    Htex.push_back(H);
    textures.push_back(texture);

}


void TriangleMesh::translate(float tx, float ty, float tz) {
    for (size_t i = 0; i < vertices.size(); i++) {
        vertices[i][0] += tx; // > 0 vers la droite
        vertices[i][1] += ty; // > 0 vers le haut
        vertices[i][2] += tz; // > 0 vers l'avant
    }
}

void TriangleMesh::translateRelative(float tx, float ty, float tz) {
    Vector dRight = tx * right;
    Vector dUp = ty * up;
    Vector dFwd = tz* fwd;
    Vector d = tx * right + ty * up + tz * fwd;
    for (size_t i = 0; i < vertices.size(); i++) {
        vertices[i] = vertices[i] + tx * right + ty * up + tz * fwd;

    }
}


void TriangleMesh::inverse(bool x, bool y, bool z) {
    if (x) {
        for (size_t i = 0; i < vertices.size(); i++) {
            vertices[i][0] *= -1;
        }

        for (size_t i = 0; i < normals.size(); i++) {
            normals[i][0] *= -1;
        }
    }

    if (y) {
        for (size_t i = 0; i < vertices.size(); i++) {
            vertices[i][1] *= -1;
        }

        for (size_t i = 0; i < normals.size(); i++) {
            normals[i][1] *= -1;
        }
    }

    if (z) {
        for (size_t i = 0; i < vertices.size(); i++) {
            vertices[i][2] *= -1;
        }

        for (size_t i = 0; i < normals.size(); i++) {
            normals[i][2] *= -1;
        }
    }

}

void TriangleMesh::rotate(float rightRotation, float upRotation, float fwdRotation) {

    float rotationRight = rightRotation * M_PI / 180;
    float rotationUp = upRotation * M_PI / 180;
    float rotationFwd = fwdRotation * M_PI / 180;


    for (size_t i = 0; i < vertices.size(); i++) {
        vertices[i] = vertices[i].getRotatedVector(right, rotationRight);
    }

    for (size_t i = 0; i < normals.size(); i++) {
        normals[i] = normals[i].getRotatedVector(right, rotationRight);
    }

    up = up.getRotatedVector(right, rotationRight);
    fwd = fwd.getRotatedVector(right, rotationRight);



    for (size_t i = 0; i < vertices.size(); i++) {
        vertices[i] = vertices[i].getRotatedVector(up, rotationUp);
    }

    for (size_t i = 0; i < normals.size(); i++) {
        normals[i] = normals[i].getRotatedVector(up, rotationUp);
    }

    up = up.getRotatedVector(up, rotationUp);
    fwd = fwd.getRotatedVector(up, rotationUp);


    for (size_t i = 0; i < vertices.size(); i++) {
        vertices[i] = vertices[i].getRotatedVector(fwd, rotationFwd);
    }

    for (size_t i = 0; i < normals.size(); i++) {
        normals[i] = normals[i].getRotatedVector(fwd, rotationFwd);
    }

    up = up.getRotatedVector(fwd, rotationFwd);
    right = right.getRotatedVector(fwd, rotationFwd);
}


void TriangleMesh::resize(float k, float l, float m) {

    for (size_t i = 0; i < vertices.size(); i++) {
        vertices[i] = Vector(k, l, m) * vertices[i];
    }

    for (size_t i = 0; i < normals.size(); i++) {
        normals[i] = Vector(k, l, m) * normals[i];
    }
}

void TriangleMesh::transformMesh(float tx, float ty, float tz, float rx, float ry, float rz, bool mirrorx, bool mirrory, bool mirrorz, float rFactorX, float rFactorY, float rFactorZ) {
    rotate(rx, ry, rz);
    inverse(mirrorx, mirrory, mirrorz);
    resize(rFactorX, rFactorY, rFactorZ);
    translate(tx, ty, tz);


    BVH = new Noeud;
    buildBVH(BVH, 0, indices.size());
}



void TriangleMesh::translateMesh(float tx, float ty, float tz) {
    translate(tx, ty, tz);
    BVH = new Noeud;
    buildBVH(BVH, 0, indices.size());
}

void TriangleMesh::translateRelativeMesh(float tx, float ty, float tz) {
    translateRelative(tx, ty, tz);
    BVH = new Noeud;
    buildBVH(BVH, 0, indices.size());
}

void TriangleMesh::inverseMesh(bool x, bool y, bool z) {
    inverse(x, y, z);
    BVH = new Noeud;
    buildBVH(BVH, 0, indices.size());
}

void TriangleMesh::rotateMesh(float rightRotation, float upRotation, float fwRotation) {
    rotateMesh(rightRotation, upRotation, fwRotation);
    BVH = new Noeud;
    buildBVH(BVH, 0, indices.size());
}

void TriangleMesh::resizeMesh(float k, float l, float m) {
    resize(k, l, m);
    BVH = new Noeud;
    buildBVH(BVH, 0, indices.size());
}
#pragma once
#include "atom.h"
#include "math/vec3.h"
#include <vector>

using std::vector;

class Cell {
    private:
        vec3 m_size;
        vector<Atom *> m_atomsClose;

    public:
        int nx;
        int ny;
        int nz;
        vec3 position;
        int cellIndex;
        bool calculatedLocally = false; // Avoid computing the forces inside a cell several times.

        Cell();
        ~Cell();
        void addAtom(Atom *atom);
        void clearList();
        bool isInCell(vec3 pos, double rcut);
        void setSize(int nx, int ny, int nz) {m_size = vec3(nx, ny, nz);}
        double calculateLocally(double sigma, double epsilon, double rcut);
        vec3 getSize() {return m_size;}
        vector<Atom *> atomsClose() {return m_atomsClose;}
};

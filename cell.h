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
        vec3 position; // Position of Cell.
        bool calculatedLocally = false;

        Cell();
        ~Cell();
        void addAtom(Atom *atom);
        void clearList();
        bool isInCell(vec3 pos, double rcut);
        void setSize(double length) {m_size = vec3(length, length, length);}
        double calculateLocally(double sigma, double epsilon, double rcut);
        vec3 getSize() {return m_size;}
        vector<Atom *> atomsClose() {return m_atomsClose;}
};

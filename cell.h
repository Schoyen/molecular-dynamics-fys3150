#pragma once
#include "atom.h"
#include "math/vec3.h"
#include <vector>

class CellList;
using std::vector;

class Cell {
    private:
        vec3 size;
        vector<Atom*> m_atomsClose;

    public:
        vec3 position; // Position of Cell.

        Cell();
        ~Cell();
        void addAtom(Atom *atom);
        void clearList();
        void setSize(double length) {size = vec3(length, length, length);}
};

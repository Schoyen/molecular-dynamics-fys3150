#pragma once
#include "atom.h"
#include "math/vec3.h"
#include <vector>

using std::vector;

class Cell {
    private:
        vector<Atom *> m_atomsClose;

    public:

        Cell();
        ~Cell();
        void addAtom(Atom *atom);
        void clearList();
        vector<Atom *> atomsClose() {return m_atomsClose;}
};

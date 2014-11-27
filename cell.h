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
        int cellIndex;

        Cell();
        ~Cell();
        void addAtom(Atom *atom);
        void addIndicesOnAtoms();
        void clearList();
        void setSize(int x, int y, int z) {m_size = vec3(x, y, z);}
        vec3 getSize() {return m_size;}
        vector<Atom *> atomsClose() {return m_atomsClose;}
};

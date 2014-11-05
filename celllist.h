#pragma once
#include "math/vec3.h"
#include "system.h"
#include "atom.h"
#include <vector>

class CellList {
    private:
        System *m_system;
        vector<Atom*> m_atomsClose;
        double m_rcut;

    public:
        CellList(System *system);
        ~CellList();
        void setEmptyList();
        void calculateCellAtoms();

        vector<Atom *> &atomsClose() {return m_atomsClose;}
};

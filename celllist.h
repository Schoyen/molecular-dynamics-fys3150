#pragma once
#include "math/vec3.h"
#include "system.h"
#include "atom.h"
#include "cell.h"
#include <vector>

class System;

class CellList {
    private:
        vector<Cell *> m_listOfCells;
        System *m_system;
        double m_rcut;

    public:
        CellList();
        ~CellList();
        void createCell(vec3 position);
        void calculateCellAtoms();
        double getrcut() {return m_rcut;}
        void setSystem(System *system) {m_system = system;}
        void setrcut(double val) {m_rcut = val;}
};

#pragma once
#include "math/vec3.h"
#include "system.h"
#include "atom.h"
#include "cell.h"
#include <vector>

class CellList {
    private:
        System *m_system;
        vector<Cell *> m_listOfCells;
        double m_rcut;

    public:
        CellList();
        ~CellList();
        void createCells();
        void calculateCellAtoms();
        void setSystem(System *system) {m_system = system;}
        double getrcut() {return m_rcut;}
};

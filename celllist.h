#pragma once
#include "math/vec3.h"
#include "system.h"
#include "atom.h"
#include "cell.h"
#include <vector>

class CellList {
    private:
        vector<Cell *> m_listOfCells;
        double m_rcut;

    public:
        CellList();
        ~CellList();
        void createCell();
        void calculateCellAtoms();
        double getrcut() {return m_rcut;}
        void setrcut(double val) {m_rcut = val;}
};

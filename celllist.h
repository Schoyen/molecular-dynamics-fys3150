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
        int numberOfCellsX;
        int numberOfCellsY;
        int numberOfCellsZ;

    public:
        CellList();
        ~CellList();
        // Storing cells by index.
        void createCell(int i, int j, int k, int ind, int nx, int ny, int nz);
        Cell *getCell(int i, int j, int k);
        void calculateCellAtoms();
        void emptyCells();
        double getrcut() {return m_rcut;}
        void setSystem(System *system) {m_system = system;}
        void setrcut(double val) {m_rcut = val;}
        vector<Cell *> listOfCells() {return m_listOfCells;}
};

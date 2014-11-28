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
        int numberOfCellsX;
        int numberOfCellsY;
        int numberOfCellsZ;

    public:
        CellList();
        ~CellList();
        // Storing cells by index.
        void setNumberOfCells(int nx, int ny, int nz);
        void createCell();
        Cell *getCell(int i, int j, int k);
        void calculateCellAtoms();
        void emptyCells();
        void setSystem(System *system) {m_system = system;}
        vector<Cell *> listOfCells() {return m_listOfCells;}
};

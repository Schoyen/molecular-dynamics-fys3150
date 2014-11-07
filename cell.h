#pragma once
#include "atom.h"
#include "math/vec3.h"
#include <vector>

class CellList;

class Cell {
    private:
        CellList *m_cellList;
        std::vector<Atom*> m_atomsClose;

    public:
        vec3 position;

        Cell();
        ~Cell();
        void addAtom(Atom *atom);
        void setEmptyList();
        void setCellList(CellList *celllist) {m_cellList = celllist;}
        std::vector<Atom *> &atoms() {return m_atomsClose;}
};

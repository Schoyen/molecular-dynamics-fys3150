#include "celllist.h"

CellList::CellList() :
    m_listOfCells(0)
{

}

CellList::~CellList() {
    m_listOfCells.clear();
}

void CellList::createCell() {
    Cell *cell = new Cell();
    cell->setSize(m_rcut);
    m_listOfCells.push_back(cell);
}

// Implement the method.
void CellList::calculateCellAtoms() {
}

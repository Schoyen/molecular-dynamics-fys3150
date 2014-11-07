#include "celllist.h"

CellList::CellList() :
    m_system(0),
    m_listOfCells(0)
{

}

CellList::~CellList() {
    delete m_system;
    m_listOfCells.clear();
}

void createCells() {

}

// Implement the method.
void calculateCellAtoms() {
}

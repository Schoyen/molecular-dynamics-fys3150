#include "celllist.h"

CellList::CellList(System *system) {
    m_system = system;
}

CellList::~CellList() {
    delete m_system;
    m_atomsClose.clear();
}

// Implement the method.
void setEmptyList() {
}

// Implement the method.
void calculateCellAtoms() {
}

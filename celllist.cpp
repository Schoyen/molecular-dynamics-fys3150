#include "celllist.h"
#include <cmath>

CellList::CellList() :
    m_listOfCells(0),
    m_system(0)
{

}

CellList::~CellList() {
    m_listOfCells.clear();
}

void CellList::createCell(vec3 pos) {
    Cell *cell = new Cell();
    cell->setSize(m_rcut);
    cell->position = pos;
    m_listOfCells.push_back(cell);
}

void CellList::calculateCellAtoms() {
    double distance;
    vec3 temp;
    for (int i = 0; i < (int) m_system->atoms().size(); i++) {
        for (int j = 0; j < (int) m_system->atoms().size(); j++) {
            temp = m_listOfCells[j]->position - m_system->atoms()[i]->position;
            distance = temp.lengthSquared();
            if (distance < pow(m_rcut, 2.0) + pow(m_rcut, 2.0) + pow(m_rcut, 2.0)) {
                m_listOfCells[j]->addAtom(m_system->atoms()[i]);
                break;
            }
        }
    }
}

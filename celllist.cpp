#include "celllist.h"
#include <cmath>

CellList::CellList() :
    m_listOfCells(0),
    m_system(0)
{

}

CellList::~CellList()
{
    m_listOfCells.clear();
}

void CellList::createCell(vec3 pos)
{
    Cell *cell = new Cell();
    cell->setSize(m_rcut);
    cell->position = pos;
    m_listOfCells.push_back(cell);
}

void CellList::emptyCells()
{
    for (int i = 0; i < (int) m_listOfCells.size(); i++) {
        m_listOfCells[i]->clearList();
    }
}

/*
 * Calcualting the center of each atom and checking whether or not the positions of the atoms are inside the cell.
 */
void CellList::calculateCellAtoms()
{
    double distance;
    vec3 temp;
    vec3 center = vec3(m_rcut/2.0, m_rcut/2.0, m_rcut/2.0);
    for (int i = 0; i < (int) m_system->atoms().size(); i++) {
        for (int j = 0; j < (int) m_system->atoms().size(); j++) {
            temp = m_listOfCells[j]->position + center;
            temp = m_system->atoms()[i]->position - temp;
            distance = temp.lengthSquared();
            if (distance < m_listOfCells[j]->getSize().lengthSquared()/2.0) {
                m_listOfCells[j]->addAtom(m_system->atoms()[i]);
                break;
            }
        }
    }
}

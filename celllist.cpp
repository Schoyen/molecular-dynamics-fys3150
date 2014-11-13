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
    //double distance;
    // The atoms seems to be uniformly distributed among the cells.
    vec3 temp;
    vec3 center = vec3(m_rcut/2.0, m_rcut/2.0, m_rcut/2.0);
    for (int i = 0; i < (int) m_system->atoms().size(); i++) {
        for (int j = 0; j < (int) m_system->atoms().size(); j++) {
            temp = m_listOfCells[j]->position + center;
            temp = m_system->atoms()[i]->position - temp;
            //distance = temp.lengthSquared();
            if (temp.x() <= center.x() && temp.y() <= center.y() && temp.z() <= center.z()) {
                m_listOfCells[j]->addAtom(m_system->atoms()[i]);
                break;
            }
        }
    }
    int tempt = 0;
    for (int i = 0; i < (int) m_listOfCells.size(); i++) {
        tempt += m_listOfCells[i]->atomsClose().size();
        std::cout << m_listOfCells[i]->atomsClose().size() << std::endl;
    }
    std::cout << tempt << std::endl;
    /*
    temp = m_listOfCells[0]->position + center;
    std::cout << temp << std::endl;
    std::cout << temp.lengthSquared() << std::endl;
    for (int i = 0; i < (int) m_listOfCells[0]->atomsClose().size(); i++) {
        std::cout << i << "\t" << m_listOfCells[0]->atomsClose()[i]->position << std::endl;
    }
    */
}

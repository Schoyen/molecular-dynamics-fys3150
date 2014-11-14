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
 * Calculating the center of each atom and checking whether or not the positions of the atoms are inside the cell.
 * 
 * CHECK: Are the cells overlapping?
 *
 */
void CellList::calculateCellAtoms()
{
    //double distance;
    // The atoms seems to be uniformly distributed among the cells.
    // Negative Houston.
    vec3 temp;
    //vec3 center = vec3(m_rcut/2.0, m_rcut/2.0, m_rcut/2.0);
    vec3 positionOfCenter;
    vec3 cellPos;
    vec3 atomPos;

    for (int i = 0; i < (int) m_system->atoms().size(); i++) {
        // Iterating over all atoms.
        for (int j = 0; j < (int) m_listOfCells.size(); j++) {
            cellPos = m_listOfCells[j]->position;
            atomPos = m_system->atoms()[i]->position;
            if (cellPos.x() <= atomPos.x() && atomPos.x() < cellPos.x() + m_rcut &&
                cellPos.y() <= atomPos.y() && atomPos.y() < cellPos.y() + m_rcut &&
                cellPos.z() <= atomPos.z() && atomPos.z() < cellPos.z() + m_rcut) {
                m_listOfCells[j]->addAtom(m_system->atoms()[i]);
            }
            /*
            // Iterating over all cells in cell list.
            positionOfCenter = m_listOfCells[j]->position + center;
            temp = m_system->atoms()[i]->position - positionOfCenter;
            */
        }
    }



    /*
    for (int i = 0; i < (int) m_system->atoms().size(); i++) {
        for (int j = 0; j < (int) m_system->atoms().size(); j++) {
            temp = m_listOfCells[j]->position + center;
            temp = m_system->atoms()[i]->position - temp;
            //distance = temp.lengthSquared();
            if (temp.x() < center.x() && temp.y() < center.y() && temp.z() < center.z()) {
                m_listOfCells[j]->addAtom(m_system->atoms()[i]);
                break;
            }
        }
    }
    */

    
    /*
    int tempt = 0;
    for (int i = 0; i < (int) m_listOfCells.size(); i++) {
        tempt += m_listOfCells[i]->atomsClose().size();
        temp = m_listOfCells[i]->position;
        std::cout << m_listOfCells[i]->atomsClose().size() << "\t" << temp <<  std::endl;
    }
    std::cout << tempt << std::endl;
    */
    /*
    std::cout << tempt << "\n" << std::endl;
    temp = m_listOfCells[0]->position + center;
    std::cout << temp << std::endl;
    std::cout << temp.lengthSquared() << std::endl;
    for (int i = 0; i < (int) m_listOfCells[0]->atomsClose().size(); i++) {
        std::cout << i << "\t" << m_listOfCells[0]->atomsClose()[i]->position << std::endl;
    }
    */

    /*
    temp = m_listOfCells[0]->position + m_rcut;
    std::cout << m_listOfCells[0]->position << "\t" << temp << std::endl;
    for (int i = 0; i < (int) m_listOfCells[0]->atomsClose().size(); i++) {
        std::cout << m_listOfCells[0]->atomsClose()[i]->position << std::endl;
    }
    temp = m_listOfCells[1]->position + m_rcut;
    std::cout << "\n" << m_listOfCells[1]->position << "\t" << temp << std::endl;
    for (int i = 0; i < (int) m_listOfCells[1]->atomsClose().size(); i++) {
        std::cout << m_listOfCells[1]->atomsClose()[i]->position << std::endl;
    }
    temp = m_listOfCells[3]->position + m_rcut;
    std::cout << "\n" << m_listOfCells[3]->position << "\t" << temp << std::endl;
    for (int i = 0; i < (int) m_listOfCells[3]->atomsClose().size(); i++) {
        std::cout << m_listOfCells[3]->atomsClose()[i]->position << std::endl;
    }
    */
    //std::cout << "Calculating cell atoms" << std::endl;
}

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

/*
 * Method used by system upon creation of the fcc lattices.
 */
void CellList::createCell(int i, int j, int k, int ind, int nx, int ny, int nz)
{
    Cell *cell = new Cell();

    // The size might not be relevant.
    // nx, ny and nz determines the number of cells, not the size.
    cell->setSize(nx, ny, nz);
    
    // Method determining index of cell in cell list.
    cell->positionInList = i * ny * nz + j * nz + k;
    numberOfCellsX = nx;
    numberOfCellsY = ny;
    numberOfCellsZ = nz;
    
    // Value used to determine indices of atoms.
    // This will then again be used to avoid calculating forces between
    // atoms in cells twice. We will then use N3L.
    cell->index = ind;

    m_listOfCells.push_back(cell);
}

Cell *CellList::getCell(int i, int j, int k)
{
    return m_listOfCells[i * numberOfCellsY * numberOfCellsZ + j * numberOfCellsZ + k];
}

/*
 * Removing the atoms from the cells.
 */
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
    vec3 temp;
    //vec3 center = vec3(m_rcut/2.0, m_rcut/2.0, m_rcut/2.0);
    vec3 positionOfCenter;
    vec3 cellPos;
    vec3 atomPos;

    for (int i = 0; i < (int) m_system->atoms().size(); i++) {
        atomPos = m_system->atoms()[i]->position;
        // Iterating over all atoms.
        for (int j = 0; j < (int) m_listOfCells.size(); j++) {
            //cellPos = m_listOfCells[j]->position;
            if (m_listOfCells[j]->isInCell(atomPos, m_rcut)) {
                m_listOfCells[j]->addAtom(m_system->atoms()[i]);
            }
        }
    }


    // A check to see if all the atoms inside the volume of the cell are added.
    /*
    std::cout << m_listOfCells[0]->position << std::endl;
    for (int i = 0; i < (int) m_system->atoms().size(); i++) {
        atomPos = m_system->atoms()[i]->position;
        if (m_listOfCells[0]->isInCell(atomPos, m_rcut)) {
            std::cout << atomPos << std::endl;
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
    //std::cout << "Calculating cell atoms" << std::endl;
    */
}

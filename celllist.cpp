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
void CellList::createCell(int ind, int nx, int ny, int nz)
{
    Cell *cell = new Cell();

    // The size might not be relevant.
    // nx, ny and nz determines the number of cells, not the size.
    cell->setSize(m_system->systemSize().x() / nx, m_system->systemSize().y() / ny, m_system->systemSize().z() / nz);
    
    // Method determining index of cell in cell list.
    numberOfCellsX = nx;
    numberOfCellsY = ny;
    numberOfCellsZ = nz;
    
    // Value used to determine indices of atoms.
    // This will then again be used to avoid calculating forces between
    // atoms in cells twice. We will then use N3L.
    cell->cellIndex = ind;

    // The cells "should" be stored by the index = i * ny * nz + j * nz + k;
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

void CellList::calculateCellAtoms()
{
    int cx, cy, cz; // Cell index in x, y and z direction.
    Cell *c;
    for (int i = 0; i < (int) m_system->atoms().size(); i++) {
        cx = int(m_system->atoms()[i]->position.x() / m_system->systemSize().x() * numberOfCellsX);
        cy = int(m_system->atoms()[i]->position.y() / m_system->systemSize().y() * numberOfCellsY);
        cz = int(m_system->atoms()[i]->position.z() / m_system->systemSize().z() * numberOfCellsZ);
        c = getCell(cx, cy, cz);
        c->addAtom(m_system->atoms()[i]);
    }
}

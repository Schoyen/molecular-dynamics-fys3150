#include "cell.h"
#include <cmath>

Cell::Cell()
{

}

Cell::~Cell() 
{

}

void Cell::clearList() {
    m_atomsClose.clear();
}

void Cell::addAtom(Atom *atom) {
    m_atomsClose.push_back(atom);
}

/*
 * Method used for adding indices on all atoms.
 * We use these indices in potentials/lennardjones.cpp when we are calculating the force.
 * This will make it possible for us to calculate the force using N3L.
 */
void Cell::addIndicesOnAtoms()
{
    for (int i = 0; i < (int) m_atomsClose.size(); i++) {
        m_atomsClose[i]->index = cellIndex + i;
    }
}

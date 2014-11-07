#include "cell.h"

Cell::Cell()
{

}

Cell::~Cell() 
{

}

void Cell::clearList() {
    m_atomsClose.clear();
}

// Implement this method.
void Cell::addAtom(Atom *atom) {
    m_atomsClose.push_back(atom);
}

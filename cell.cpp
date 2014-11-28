#include "cell.h"
#include <cmath>

Cell::Cell() :
    m_atomsClose(0)
{

}

Cell::~Cell() 
{
    m_atomsClose.clear();
}

void Cell::clearList() {
    m_atomsClose.clear();
}

void Cell::addAtom(Atom *atom) {
    m_atomsClose.push_back(atom);
}

#include "MolecularML.hpp"
#include <GraphMol/SmilesParse/SmilesParse.h>

MolecularTranslator::MolecularTranslator() {}
MolecularTranslator::~MolecularTranslator() {}

RDKit::ROMol* MolecularTranslator::parseSmiles(const std::string& smiles)
{
    try {
        return RDKit::SmilesToMol(smiles);
    } catch (...) {
        return NULL;
    }
}

bool MolecularTranslator::validate(const std::string& smiles)
{
    RDKit::ROMol* mol = parseSmiles(smiles);
    if (mol == NULL) return false;
    delete mol;
    return true;
}

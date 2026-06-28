#include "MolecularML.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/RingInfo.h>

FeatureExtractor::FeatureExtractor() {}
FeatureExtractor::~FeatureExtractor() {}

std::vector<std::string> FeatureExtractor::getFeatureNames() const
{
    std::vector<std::string> names;
    names.push_back("Atoms");
    names.push_back("Bonds");
    names.push_back("AvgMW");
    names.push_back("ExactMW");
    names.push_back("Rings");
    names.push_back("HeavyAtoms");
    names.push_back("Donors");
    names.push_back("Acceptors");
    return names;
}

std::vector<double> FeatureExtractor::extract(RDKit::ROMol* mol)
{
    std::vector<double> features;
    if (mol == NULL) return features;
    
    try {
        // Basic properties
        features.push_back(mol->getNumAtoms());
        features.push_back(mol->getNumBonds());
        
        // Molecular weight
        features.push_back(RDKit::Descriptors::calcAMW(*mol));
        features.push_back(RDKit::Descriptors::calcExactMW(*mol));
        
        // Rings
        features.push_back(mol->getRingInfo()->numRings());
        
        // Heavy atoms (non-H)
        int heavy = 0;
        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            if (mol->getAtomWithIdx(i)->getAtomicNum() > 1)
                heavy++;
        }
        features.push_back(heavy);
        
        // H-bond donors and acceptors
        int donors = 0, acceptors = 0;
        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            RDKit::Atom* atom = mol->getAtomWithIdx(i);
            int num = atom->getAtomicNum();
            
            if (num == 7 || num == 8) {
                if (atom->getNumImplicitHs() > 0 || atom->getNumExplicitHs() > 0)
                    donors++;
                if (num == 8 || (num == 7 && atom->getFormalCharge() == 0))
                    acceptors++;
            }
        }
        features.push_back(donors);
        features.push_back(acceptors);
    } catch (...) {
        return features;
    }
    
    return features;
}

std::string FeatureExtractor::formula(RDKit::ROMol* mol)
{
    if (mol == NULL) return "";
    try {
        return RDKit::Descriptors::calcMolFormula(*mol);
    } catch (...) {
        return "";
    }
}

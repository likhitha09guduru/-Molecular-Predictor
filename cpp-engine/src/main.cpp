#include "MolecularML.hpp"
#include <iostream>
#include <iomanip>

void printFeatures(const std::vector<double>& features, 
                   const std::vector<std::string>& names)
{
    std::cout << "\nExtracted Features:" << std::endl;
    std::cout << "---------------------" << std::endl;
    for (size_t i = 0; i < features.size() && i < names.size(); ++i) {
        std::cout << std::setw(12) << std::left << names[i] 
                  << ": " << std::fixed << std::setprecision(3) 
                  << features[i] << std::endl;
    }
}

int main(int argc, char* argv[])
{
    std::string smiles = (argc > 1) ? argv[1] : "CCO"; // Default: ethanol
    
    std::cout << "Molecular ML Engine" << std::endl;
    std::cout << "SMILES: " << smiles << std::endl;
    std::cout << "====================" << std::endl;
    
    // 1. Parse and validate
    MolecularTranslator translator;
    
    if (!translator.validate(smiles)) {
        std::cout << "Invalid SMILES!" << std::endl;
        return 1;
    }
    std::cout << "✓ Valid SMILES" << std::endl;
    
    RDKit::ROMol* mol = translator.parseSmiles(smiles);
    if (mol == NULL) {
        std::cout << "Failed to parse molecule" << std::endl;
        return 1;
    }
    
    // 2. Extract features
    FeatureExtractor extractor;
    std::vector<double> features = extractor.extract(mol);
    printFeatures(features, extractor.getFeatureNames());
    
    // 3. Get formula
    std::cout << "\nMolecular Formula: " << extractor.formula(mol) << std::endl;
    
    // 4. Run benchmarks
    Benchmark bench;
    bench.runAll(smiles);
    
    // Cleanup
    delete mol;
    
    return 0;
}

#include "MolecularML.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

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

// Export features to CSV for Python training
bool exportToCSV(const std::string& smiles, const std::vector<double>& features)
{
    std::ofstream file("features.csv", std::ios::app);
    if (!file.is_open()) return false;
    
    file << smiles;
    for (size_t i = 0; i < features.size(); ++i) {
        file << "," << features[i];
    }
    file << "\n";
    file.close();
    return true;
}

int main(int argc, char* argv[])
{
    std::string smiles = (argc > 1) ? argv[1] : "CCO";
    
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
    
    // 4. Export for Python training
    if (exportToCSV(smiles, features)) {
        std::cout << "\n✓ Features exported to features.csv" << std::endl;
    }
    
    // 5. Run benchmarks
    Benchmark bench;
    bench.printAll(smiles);
    
    // Cleanup
    delete mol;
    
    return 0;
}

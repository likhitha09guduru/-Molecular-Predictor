#include "MolecularML.hpp"
#include <iostream>
#include <iomanip>

Benchmark::Benchmark() {}
Benchmark::~Benchmark() {}

double Benchmark::measureParse(const std::string& smiles)
{
    MolecularTranslator translator;
    const int RUNS = 10;
    double total = 0.0;
    
    for (int i = 0; i < RUNS; ++i) {
        clock_t start = clock();
        RDKit::ROMol* mol = translator.parseSmiles(smiles);
        clock_t end = clock();
        if (mol != NULL) delete mol;
        total += static_cast<double>(end - start) / CLOCKS_PER_SEC;
    }
    return total / RUNS;
}

double Benchmark::measureExtract(RDKit::ROMol* mol)
{
    if (mol == NULL) return 0.0;
    
    FeatureExtractor extractor;
    const int RUNS = 10;
    double total = 0.0;
    
    for (int i = 0; i < RUNS; ++i) {
        clock_t start = clock();
        std::vector<double> features = extractor.extract(mol);
        clock_t end = clock();
        total += static_cast<double>(end - start) / CLOCKS_PER_SEC;
    }
    return total / RUNS;
}

double Benchmark::measureValidate(const std::string& smiles)
{
    MolecularTranslator translator;
    const int RUNS = 10;
    double total = 0.0;
    
    for (int i = 0; i < RUNS; ++i) {
        clock_t start = clock();
        translator.validate(smiles);
        clock_t end = clock();
        total += static_cast<double>(end - start) / CLOCKS_PER_SEC;
    }
    return total / RUNS;
}

double Benchmark::measureFullPipeline(const std::string& smiles)
{
    MolecularTranslator translator;
    FeatureExtractor extractor;
    const int RUNS = 10;
    double total = 0.0;
    
    for (int i = 0; i < RUNS; ++i) {
        clock_t start = clock();
        RDKit::ROMol* mol = translator.parseSmiles(smiles);
        if (mol != NULL) {
            std::vector<double> features = extractor.extract(mol);
            delete mol;
        }
        clock_t end = clock();
        total += static_cast<double>(end - start) / CLOCKS_PER_SEC;
    }
    return total / RUNS;
}

void Benchmark::printAll(const std::string& smiles)
{
    std::cout << "\n=== Benchmark Results (avg of 10 runs) ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Parse SMILES:  " << measureParse(smiles) << " sec" << std::endl;
    std::cout << "Validate:      " << measureValidate(smiles) << " sec" << std::endl;
    std::cout << "Full Pipeline: " << measureFullPipeline(smiles) << " sec" << std::endl;
}

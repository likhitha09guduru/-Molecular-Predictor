#include "MolecularML.hpp"

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

void Benchmark::runAll(const std::string& smiles)
{
    MolecularTranslator translator;
    
    std::cout << "\n=== Benchmark Results ===" << std::endl;
    std::cout << "Parse time: " << measureParse(smiles) << " sec" << std::endl;
    std::cout << "Validate time: " << measureValidate(smiles) << " sec" << std::endl;
    
    RDKit::ROMol* mol = translator.parseSmiles(smiles);
    if (mol != NULL) {
        std::cout << "Extract time: " << measureExtract(mol) << " sec" << std::endl;
        delete mol;
    }
}

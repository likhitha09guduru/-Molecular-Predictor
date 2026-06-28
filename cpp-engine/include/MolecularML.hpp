#ifndef MOLECULAR_ML_HPP
#define MOLECULAR_ML_HPP

#include <string>
#include <vector>
#include <ctime>
#include <GraphMol/ROMol.h>

class MolecularTranslator
{
public:
    MolecularTranslator();
    ~MolecularTranslator();
    
    RDKit::ROMol* parseSmiles(const std::string& smiles);
    bool validate(const std::string& smiles);
};

class FeatureExtractor
{
public:
    FeatureExtractor();
    ~FeatureExtractor();
    
    std::vector<double> extract(RDKit::ROMol* mol);
    std::string formula(RDKit::ROMol* mol);
    std::vector<std::string> getFeatureNames() const;
    int getFeatureCount() const;
};

class Benchmark
{
public:
    Benchmark();
    ~Benchmark();
    double measureParse(const std::string& smiles);
    double measureExtract(RDKit::ROMol* mol);
    double measureValidate(const std::string& smiles);
    double measureFullPipeline(const std::string& smiles);
    void printAll(const std::string& smiles);
};

#endif

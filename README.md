# -Molecular-Predictor
🔬 High-performance molecular property prediction from SMILES using RDKit, Random Forest, C++, and FastAPI.
# Molecular ML - Complete Machine Learning Pipeline

## Overview
A complete C++/Python pipeline for molecular property prediction using SMILES strings.

## Architecture
- **C++ Engine**: Fast feature extraction using RDKit
- **Python ML**: Training and prediction using scikit-learn
- **Integration**: C++ exports features, Python trains models

## Features
- SMILES parsing and validation
- 8 molecular descriptors extraction
- Random Forest regression model
- Batch prediction support
- Performance benchmarking

## Build and Run

### 1. Build C++ Engine
```bash
make

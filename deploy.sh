#!/bin/bash
# Deployment script for molecular-predictor

echo "🚀 Deploying Molecular Predictor..."

# Build C++ engine
echo "Building C++ engine..."
cd cpp-engine
make || echo "C++ build skipped"
cd ..

# Check if model exists
if [ ! -f "model.pkl" ]; then
    echo "⚠️ Model not found. Training..."
    python train.py molecules.csv
fi

# Start API server
echo "🌐 Starting API server..."
python api.py

echo "✅ Deployment complete!"
echo "API running at: http://localhost:8000"
echo "API documentation: http://localhost:8000/docs"

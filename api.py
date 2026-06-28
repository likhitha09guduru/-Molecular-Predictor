#!/usr/bin/env python3
"""
FastAPI server for molecular property prediction
"""

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import pickle
import numpy as np
import subprocess
import os
from typing import List, Optional

app = FastAPI(
    title="Molecular Predictor API",
    description="Predict molecular properties from SMILES strings using Random Forest",
    version="1.0.0"
)

# Load model at startup
MODEL_PATH = "model.pkl"
model = None
scaler = None

class PredictionRequest(BaseModel):
    smiles: str
    features: Optional[List[float]] = None

class PredictionResponse(BaseModel):
    smiles: str
    prediction: float
    features: List[float]
    valid: bool

class BatchPredictionRequest(BaseModel):
    smiles_list: List[str]

@app.on_event("startup")
async def load_model():
    """Load the trained model on startup"""
    global model, scaler
    if os.path.exists(MODEL_PATH):
        with open(MODEL_PATH, 'rb') as f:
            data = pickle.load(f)
            model = data['model']
            scaler = data['scaler']
        print("✓ Model loaded successfully")
    else:
        print("⚠️ No model found. Train first: python train.py molecules.csv")

def extract_features_from_cpp(smiles: str):
    """Call C++ engine to extract features"""
    try:
        # Path to C++ executable
        cpp_engine = "./cpp-engine/molecular_ml"
        
        # If C++ engine doesn't exist, use Python fallback
        if not os.path.exists(cpp_engine):
            return extract_features_python(smiles)
        
        result = subprocess.run(
            [cpp_engine, smiles],
            capture_output=True,
            text=True,
            timeout=5
        )
        
        # Parse features from output
        features = []
        lines = result.stdout.split('\n')
        found_features = False
        
        for line in lines:
            if 'Extracted Features:' in line:
                found_features = True
                continue
            if found_features and ':' in line and '---' not in line:
                parts = line.split(':')
                if len(parts) == 2:
                    try:
                        value = float(parts[1].strip().split()[0])
                        features.append(value)
                    except:
                        pass
        
        # Should have 8 features
        return features if len(features) == 8 else None
        
    except Exception as e:
        print(f"C++ engine error: {e}")
        return None

def extract_features_python(smiles: str):
    """Python fallback for feature extraction"""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        features = [
            float(mol.GetNumAtoms()),
            float(mol.GetNumBonds()),
            Descriptors.MolWt(mol),
            Descriptors.ExactMolWt(mol),
            float(Descriptors.RingCount(mol)),
            float(Descriptors.HeavyAtomCount(mol)),
            float(Descriptors.NumHDonors(mol)),
            float(Descriptors.NumHAcceptors(mol))
        ]
        return features
    except:
        return None

@app.get("/")
async def root():
    return {
        "message": "Molecular Predictor API",
        "endpoints": {
            "/": "This page",
            "/health": "Health check",
            "/predict": "Predict single molecule",
            "/predict/batch": "Predict multiple molecules",
            "/docs": "API documentation"
        }
    }

@app.get("/health")
async def health():
    if model is None:
        return {"status": "error", "message": "Model not loaded"}
    return {"status": "healthy", "model_loaded": True}

@app.post("/predict", response_model=PredictionResponse)
async def predict(request: PredictionRequest):
    """Predict property for a single molecule"""
    if model is None:
        raise HTTPException(status_code=503, detail="Model not loaded")
    
    # Get features
    if request.features:
        features = request.features
    else:
        features = extract_features_from_cpp(request.smiles)
        if features is None:
            features = extract_features_python(request.smiles)
    
    if features is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    
    # Scale and predict
    features_scaled = scaler.transform(np.array(features).reshape(1, -1))
    prediction = model.predict(features_scaled)[0]
    
    return PredictionResponse(
        smiles=request.smiles,
        prediction=float(prediction),
        features=features,
        valid=True
    )

@app.post("/predict/batch")
async def predict_batch(request: BatchPredictionRequest):
    """Predict properties for multiple molecules"""
    if model is None:
        raise HTTPException(status_code=503, detail="Model not loaded")
    
    results = []
    for smiles in request.smiles_list:
        features = extract_features_from_cpp(smiles)
        if features is None:
            features = extract_features_python(smiles)
        
        if features is None:
            results.append({
                "smiles": smiles,
                "error": "Invalid SMILES"
            })
            continue
        
        features_scaled = scaler.transform(np.array(features).reshape(1, -1))
        prediction = model.predict(features_scaled)[0]
        
        results.append({
            "smiles": smiles,
            "prediction": float(prediction),
            "features": features
        })
    
    return {"results": results}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)

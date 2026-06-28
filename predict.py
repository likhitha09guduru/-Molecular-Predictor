#!/usr/bin/env python3
"""
Prediction Script for Molecular Properties
"""

import pandas as pd
import numpy as np
import pickle
import subprocess
import sys
import os

class MolecularMLPredictor:
    def __init__(self, model_file='model.pkl'):
        self.model = None
        self.scaler = None
        self.feature_names = None
        self.load_model(model_file)
    
    def load_model(self, model_file):
        """Load trained model"""
        if not os.path.exists(model_file):
            print(f"Error: {model_file} not found!")
            print("Train a model first: python train.py")
            sys.exit(1)
        
        with open(model_file, 'rb') as f:
            data = pickle.load(f)
            self.model = data['model']
            self.scaler = data['scaler']
            self.feature_names = data['feature_names']
        
        print(f"✓ Model loaded from {model_file}")
        print(f"Features: {', '.join(self.feature_names)}")
    
    def _get_features_from_cpp(self, smiles):
        """Call C++ engine to extract features"""
        try:
            result = subprocess.run(
                ['./molecular_ml', smiles],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if 'Invalid SMILES' in result.stdout:
                print(f"Invalid SMILES: {smiles}")
                return None
            
            features = self._parse_features(result.stdout)
            return features
            
        except Exception as e:
            print(f"Error processing {smiles}: {e}")
            return None
    
    def _parse_features(self, output):
        """Parse features from C++ output"""
        features = []
        lines = output.split('\n')
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
        
        return features if len(features) == 8 else None
    
    def predict_single(self, smiles):
        """Predict property for a single molecule"""
        features = self._get_features_from_cpp(smiles)
        if features is None:
            return None
        
        # Scale and predict
        features_scaled = self.scaler.transform(np.array(features).reshape(1, -1))
        prediction = self.model.predict(features_scaled)[0]
        
        return {
            'smiles': smiles,
            'prediction': prediction,
            'features': features
        }
    
    def predict_batch(self, smiles_list):
        """Predict properties for multiple molecules"""
        results = []
        
        for smiles in smiles_list:
            result = self.predict_single(smiles)
            if result:
                results.append(result)
        
        return results
    
    def print_prediction(self, result):
        """Pretty print prediction result"""
        if result is None:
            return
        
        print("\n" + "="*50)
        print(f"SMILES: {result['smiles']}")
        print("-"*50)
        print(f"Predicted Property: {result['prediction']:.4f}")
        print("\nFeatures:")
        for name, value in zip(self.feature_names, result['features']):
            print(f"  {name:12s}: {value:.3f}")

def main():
    print("=== Molecular ML Predictor ===\n")
    
    # Check if C++ engine exists
    if not os.path.exists('./molecular_ml'):
        print("Error: C++ engine not found! Run 'make' first.")
        return 1
    
    # Check if model exists
    if not os.path.exists('model.pkl'):
        print("Error: model.pkl not found!")
        print("Train a model first: python train.py")
        return 1
    
    predictor = MolecularMLPredictor()
    
    # Single prediction
    if len(sys.argv) == 2:
        smiles = sys.argv[1]
        result = predictor.predict_single(smiles)
        predictor.print_prediction(result)
    
    # Batch prediction from file
    elif len(sys.argv) == 3 and sys.argv[1] == '--file':
        filename = sys.argv[2]
        if not os.path.exists(filename):
            print(f"Error: {filename} not found!")
            return 1
        
        with open(filename, 'r') as f:
            smiles_list = [line.strip() for line in f if line.strip()]
        
        print(f"Processing {len(smiles_list)} molecules...")
        results = predictor.predict_batch(smiles_list)
        
        # Save results
        output_file = 'predictions.csv'
        with open(output_file, 'w') as f:
            f.write("smiles,prediction\n")
            for r in results:
                f.write(f"{r['smiles']},{r['prediction']:.6f}\n")
        
        print(f"\n✓ Predictions saved to {output_file}")
        
        # Print summary
        print("\nSummary:")
        for r in results[:5]:  # Show first 5
            print(f"  {r['smiles']:15s}: {r['prediction']:.4f}")
        if len(results) > 5:
            print(f"  ... and {len(results)-5} more")
    
    else:
        print("Usage:")
        print("  python predict.py <SMILES>          # Single prediction")
        print("  python predict.py --file <file>     # Batch prediction")
        print("\nExample:")
        print("  python predict.py CCO")
        print("  python predict.py --file smiles.txt")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())

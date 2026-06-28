#!/usr/bin/env python3
"""
Machine Learning Training Script for Molecular Properties
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
import pickle
import subprocess
import sys
import os

class MolecularMLTrainer:
    def __init__(self):
        self.model = None
        self.scaler = None
        self.feature_names = None
        
    def load_data(self, csv_file):
        """Load molecular data from CSV file"""
        if not os.path.exists(csv_file):
            print(f"Error: {csv_file} not found!")
            return None, None
        
        df = pd.read_csv(csv_file)
        print(f"Loaded {len(df)} molecules from {csv_file}")
        
        # Separate features and target
        if 'target' not in df.columns:
            print("Error: CSV must contain 'target' column")
            return None, None
        
        # Features are all columns except 'smiles' and 'target'
        feature_cols = [col for col in df.columns if col not in ['smiles', 'target']]
        self.feature_names = feature_cols
        
        X = df[feature_cols].values
        y = df['target'].values
        
        print(f"Features: {', '.join(feature_cols)}")
        print(f"Target range: {y.min():.3f} - {y.max():.3f}")
        
        return X, y
    
    def generate_features_from_smiles(self, smiles_list):
        """Generate features by calling C++ engine"""
        features = []
        
        for smiles in smiles_list:
            try:
                # Call C++ program
                result = subprocess.run(
                    ['./molecular_ml', smiles],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                
                # Parse features from output
                feature_values = self._parse_features(result.stdout)
                if feature_values:
                    features.append(feature_values)
                else:
                    print(f"Warning: Could not extract features for {smiles}")
                    
            except Exception as e:
                print(f"Error processing {smiles}: {e}")
                continue
        
        return np.array(features)
    
    def _parse_features(self, output):
        """Parse feature values from C++ output"""
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
                        # Extract numeric value
                        value = float(parts[1].strip().split()[0])
                        features.append(value)
                    except:
                        pass
        
        # Should have 8 features
        return features if len(features) == 8 else None
    
    def train(self, X, y, test_size=0.2):
        """Train Random Forest model"""
        print("\n=== Training Model ===")
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=42
        )
        print(f"Training set: {len(X_train)} molecules")
        print(f"Test set: {len(X_test)} molecules")
        
        # Scale features
        self.scaler = StandardScaler()
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        # Train Random Forest
        self.model = RandomForestRegressor(
            n_estimators=100,
            max_depth=10,
            min_samples_split=5,
            random_state=42
        )
        
        self.model.fit(X_train_scaled, y_train)
        
        # Evaluate
        train_pred = self.model.predict(X_train_scaled)
        test_pred = self.model.predict(X_test_scaled)
        
        train_r2 = r2_score(y_train, train_pred)
        test_r2 = r2_score(y_test, test_pred)
        test_rmse = np.sqrt(mean_squared_error(y_test, test_pred))
        
        print(f"Training R²: {train_r2:.4f}")
        print(f"Testing R²:  {test_r2:.4f}")
        print(f"Testing RMSE: {test_rmse:.4f}")
        
        # Cross-validation
        cv_scores = cross_val_score(
            self.model, X_train_scaled, y_train, cv=5
        )
        print(f"Cross-validation R²: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
        
        # Feature importance
        self._print_feature_importance()
        
        return self.model
    
    def _print_feature_importance(self):
        """Print feature importance from Random Forest"""
        if self.model is None or self.feature_names is None:
            return
        
        print("\n=== Feature Importance ===")
        importances = self.model.feature_importances_
        indices = np.argsort(importances)[::-1]
        
        for i in indices:
            print(f"{self.feature_names[i]:12s}: {importances[i]:.4f}")
    
    def save_model(self, filename='model.pkl'):
        """Save trained model and scaler"""
        if self.model is None:
            print("No model to save!")
            return
        
        with open(filename, 'wb') as f:
            pickle.dump({
                'model': self.model,
                'scaler': self.scaler,
                'feature_names': self.feature_names
            }, f)
        print(f"\n✓ Model saved to {filename}")
    
    def load_model(self, filename='model.pkl'):
        """Load trained model and scaler"""
        with open(filename, 'rb') as f:
            data = pickle.load(f)
            self.model = data['model']
            self.scaler = data['scaler']
            self.feature_names = data['feature_names']
        print(f"✓ Model loaded from {filename}")

def generate_sample_data():
    """Generate sample molecular data for training"""
    # Common molecules with their properties (logP values)
    data = [
        ('CCO', 0.5),           # Ethanol
        ('c1ccccc1', 2.13),     # Benzene
        ('CC(=O)O', -0.17),     # Acetic acid
        ('C1CCCCC1', 2.43),     # Cyclohexane
        ('CO', -0.77),          # Methanol
        ('CCOC', 0.89),         # Diethyl ether
        ('CCN', 0.38),          # Ethylamine
        ('CCOCC', 0.89),        # Ethoxyethane
        ('c1ccncc1', 0.73),     # Pyridine
        ('CC(C)=O', -0.42),     # Acetone
        ('C1=CC=CC=C1', 2.13),  # Benzene
        ('CCCl', 1.48),         # Chloroethane
        ('CCOC(=O)C', 0.73),    # Ethyl acetate
        ('CC(C)O', 0.05),       # Isopropanol
        ('c1ccccc1O', 1.48),    # Phenol
        ('CCN(CC)CC', 2.31),    # Triethylamine
        ('CC(=O)OC', 0.18),     # Methyl acetate
        ('C1CCCC1', 2.12),      # Cyclopentane
        ('CC(=O)CC', -0.42),    # Butanone
        ('CCOCC', 0.89)         # Diethyl ether
    ]
    
    # Create CSV with features (we'll generate them using C++ engine)
    df = pd.DataFrame(columns=['smiles', 'target'])
    for smiles, target in data:
        df = df.append({'smiles': smiles, 'target': target}, ignore_index=True)
    
    df.to_csv('molecules.csv', index=False)
    print("✓ Sample data generated: molecules.csv")
    return 'molecules.csv'

def main():
    print("=== Molecular ML Training Pipeline ===\n")
    
    # Check if C++ engine exists
    if not os.path.exists('./molecular_ml'):
        print("Error: C++ engine not found! Run 'make' first.")
        return 1
    
    # Generate or load data
    if len(sys.argv) > 1:
        data_file = sys.argv[1]
    else:
        # Generate sample data if no file provided
        data_file = generate_sample_data()
    
    # Load data
    trainer = MolecularMLTrainer()
    
    # Option 1: Load pre-computed features
    if os.path.exists(data_file) and data_file.endswith('.csv'):
        X, y = trainer.load_data(data_file)
        if X is None:
            return 1
        
        # Train model
        trainer.train(X, y)
        trainer.save_model()
    
    # Option 2: Generate features from SMILES
    else:
        print("Loading SMILES from file...")
        df = pd.read_csv(data_file)
        smiles_list = df['smiles'].tolist()
        targets = df['target'].values
        
        print(f"Generating features for {len(smiles_list)} molecules...")
        X = trainer.generate_features_from_smiles(smiles_list)
        
        if len(X) == 0:
            print("No features generated!")
            return 1
        
        print(f"Generated {X.shape[0]} feature vectors with {X.shape[1]} features each")
        trainer.train(X, targets)
        trainer.save_model()
    
    print("\n✓ Training complete!")
    return 0

if __name__ == '__main__':
    sys.exit(main())

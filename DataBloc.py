import sys
import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import torch
from torch import nn
from torch.autograd import Variable
from torch.utils.data import DataLoader

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import PandasTools

# -------------------- Data Processing and preparation -------------------------

# Loading a dataset of 9000 drugs (in smiles format) with their "Binding Affinities"
# to train the model
df = pd.read_csv('dataset.csv')


# From the smiles we can calculate more than 300 DESCRIPTORS or FEATURES 
# of the Molecules or the drugs using RDKIT with generaly a simple and standard 
# line od code. 

# Using more FEATURES gives us the chance to find more correlations 
# between data. 

# In our Case, We have calculated only the most important DESCRIPTORS
# that can affect our target (the Affinity).

# Below are the 37 descriptors calculated and you can use the official docs of RDKIT
# to calculate more descriptors with same way as below. 

# Be aware that the calculated descriptors will be inserted into our dataset
# through the Dataframe (df) :


df['mol'] = df['SMILES sequence'].apply(lambda x: Chem.MolFromSmiles(x)) 
df['mol'] = df['mol'].apply(lambda x: Chem.AddHs(x))
df['num_of_atoms'] = df['mol'].apply(lambda x: x.GetNumAtoms())
df['Num_of_bonds'] = df['mol'].apply(lambda x:x.GetNumBonds())
df['num_of_heavy_atoms'] = df['mol'].apply(lambda x: x.GetNumHeavyAtoms())
df['logP'] = df['mol'].apply(lambda x: Descriptors.MolLogP(x))
df['ssr'] = df['mol'].apply(lambda x: Chem.GetSSSR(x))
df['tpsa'] = df['mol'].apply(lambda x: Descriptors.TPSA(x))
df['mol_exact_w'] = df['mol'].apply(lambda x: Descriptors.ExactMolWt(x))
df['num_valence_electrons'] = df['mol'].apply(lambda x: Descriptors.NumValenceElectrons(x))
df['num_heteroatoms'] = df['mol'].apply(lambda x: Descriptors.NumHeteroatoms(x))
df['rotatable_bond'] = df['mol'].apply(lambda x: Descriptors.NumRotatableBonds(x))
df['mol_mr'] = df['mol'].apply(lambda x: Descriptors.MolMR(x)) 
df['mol_w'] = df['mol'].apply(lambda x: Descriptors.MolWt(x))
df['NumAromaticCarbocycles'] = df['mol'].apply(lambda x: Descriptors.NumAromaticCarbocycles(x))
df['NumAromaticHeterocycles'] = df['mol'].apply(lambda x: Descriptors.NumAromaticHeterocycles(x))
df['NumAromaticRings'] = df['mol'].apply(lambda x: Descriptors.NumAromaticRings(x))
df['NumHAcceptors'] = df['mol'].apply(lambda x: Descriptors.NumHAcceptors(x))
df['NumHDonors'] = df['mol'].apply(lambda x: Descriptors.NumHDonors(x))
df['NumSaturatedRings'] = df['mol'].apply(lambda x: Descriptors.NumSaturatedRings(x))
df['NHOHCount'] = df['mol'].apply(lambda x: Descriptors.NHOHCount(x))
df['NOCount'] = df['mol'].apply(lambda x: Descriptors.NOCount(x))
df['HeavyAtomMolWt'] = df['mol'].apply(lambda x: Descriptors.HeavyAtomMolWt(x))
df['FpDensityMorgan1'] = df['mol'].apply(lambda x: Descriptors.FpDensityMorgan1(x))
df['FpDensityMorgan2'] = df['mol'].apply(lambda x: Descriptors.FpDensityMorgan2(x))
df['FpDensityMorgan3'] = df['mol'].apply(lambda x: Descriptors.FpDensityMorgan3(x))
df['fr_NH1'] = df['mol'].apply(lambda x: Descriptors.fr_NH1(x))
df['fr_NH2'] = df['mol'].apply(lambda x: Descriptors.fr_NH2(x))
df['fr_NH0'] = df['mol'].apply(lambda x: Descriptors.fr_NH0(x))
df['fr_Al_COO'] = df['mol'].apply(lambda x: Descriptors.fr_Al_COO(x))
df['fr_bicyclic'] = df['mol'].apply(lambda x: Descriptors.fr_bicyclic(x))
df['MaxPartialCharge'] = df['mol'].apply(lambda x: Descriptors.MaxPartialCharge(x))
df['HallKierAlpha'] = df['mol'].apply(lambda x: Descriptors.HallKierAlpha(x))
df['RingCount'] = df['mol'].apply(lambda x: Descriptors.RingCount(x))
df['BalabanJ'] = df['mol'].apply(lambda x: Descriptors.BalabanJ(x))

# Here we need a llop to cover all atoms of the molecule
def number_of_atoms(atom_list, df):
    for i in atom_list:
        df['num_of_{}_atoms'.format(i)] = df['mol'].apply(lambda x: len(x.GetSubstructMatches(Chem.MolFromSmiles(i))))
number_of_atoms(['C','O','N','S'], df)

# Now we need to devide the dataset of drugs molecules into Train data ( 80% ) 
# and Test data ( 20% ). The simplest way is to use the function "train_test_split"
# integrated into "sklearn" , but if you want a longer way you can use PyTorch 
# or any other labrary. 
  
from sklearn.model_selection import train_test_split

# Before dividing the dataset, we have to drop (remove) the non numeric features
# such as smiles and mol, as well as our target (Binding Affinity)

X = df.drop(columns=['SMILES sequence', 'mol', 'logP'])
y = df['logP'].values

# This is the spliting function of the library sklearn , 0.2 = 20%
x_train,x_test,y_train,y_test = train_test_split(X, y, test_size=0.2, random_state=0)

# transfering the lists resulted into numpy arrays
xn_train = np.array(x_test)
xn_test = np.array(x_test)
yn_train=np.array(y_train)
yn_test = np.array(y_test)

# transfering the numpy arrays into Torch Tensors so that we can use them 
# for training using Pytorch 

xt_train= torch.from_numpy(xn_train).float()
yt_train= torch.from_numpy(yn_train).float()

xt_test= torch.from_numpy(xn_test).float()
yt_test= torch.from_numpy(yn_test).float()

# REM : 
# To check the size of your arays and tensors use the property shape : xt_train.shape
# To check tables use head() : pd.head()


# -------------------- Regression model ---------------------------------------

# Now you are ready to integrate a good "Regression code" using Pytorch just here !



# -------------------- Model evaluation ---------------------------------------




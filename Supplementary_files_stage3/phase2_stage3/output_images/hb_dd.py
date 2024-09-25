# -*- coding: utf-8 -*-
"""HB_DD.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1kjDkluVPyxqEb9xHvrNuHY7p7x3h5TAP
"""

#HackBio Internship

!pip install chembl-webresource-client
!pip install rdkit-pypi
!pip install mordred

#import libraries
import pandas as pd
from chembl_webresource_client.new_client import new_client
import numpy as np
import rdkit
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from mordred import Calculator, descriptors

#Target search for Histone Deactylase 4
my_target = new_client.target
my_target_query = my_target.search('Histone deacetylase 4')
my_targets = pd.DataFrame.from_dict(my_target_query)
my_targets

#select and retrieve biactivity data for target
selected_target = my_targets.target_chembl_id[1]
selected_target

# Retrieve only bioactivity data for Histone Deactylase 4
activity = new_client.activity
data = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")

df = pd.DataFrame.from_dict(data)

df.head(5)

#Handling missing data
df2 = df[df.standard_value.notna()]
df2 = df2[df.canonical_smiles.notna()]
df2

len(df2.canonical_smiles.unique())

df2_nr = df2.drop_duplicates(['canonical_smiles'])
df2_nr

#Data pre-processing of the bioactivity data

#print the entire column title
columns_list = df2_nr.columns.tolist()

print(columns_list)

#Table of selected bioactivity data of interest
selection = ['molecule_chembl_id','canonical_smiles','standard_value']
df3 = df2_nr[selection]
df3

"""Labeling compunds as either being acive, inactive or intermediate"""

bioactivity_threshold = []
for i in df3.standard_value:
  if float(i) >= 10000:
    bioactivity_threshold.append("inactive")
  elif float(i) <= 1000:
    bioactivity_threshold.append("active")
  else:
    bioactivity_threshold.append("intermediate")

df3.head()

bioactivity_class = pd.Series(bioactivity_threshold, name='class')
df4 = pd.concat([df3, bioactivity_class], axis=1)
df4.head()

"""Calculate Lipinski descriptors

"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import pandas as pd
import numpy as np

# Inspired by: https://codeocean.com/explore/capsules?query=tag:data-curation
# Addition of Lipinski Descriptors
def lipinski(smiles, verbose=False):

    moldata= []
    for elem in smiles:
        # Check if elem is a valid SMILES string before converting
        if isinstance(elem, str):
            mol=Chem.MolFromSmiles(elem)
            moldata.append(mol)

    baseData= np.arange(1,1)
    i=0
    for mol in moldata:

        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])

        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1

    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)

    return descriptors

#Run lipinski using canonical_smiles
df_lipinski = lipinski(df4.canonical_smiles)

#View lipinski Descriptors
df_lipinski.head()

df4.head()

#Combine df4 table with lipinski value table
df_combined = pd.concat([df4,df_lipinski], axis=1)

#View new combined table
df_combined.head()

#Describe combines table with standard value
df_combined.standard_value.describe()

"""Convert IC50 to pIC50"""

#Import library
import numpy as np
#convert IC50 to pIC50
def pIC50(input):
    pIC50 = []

    for i in input['standard_value']:
        try: # this will try to convert i to a float
            i = float(i)
        except ValueError: # if i is not a number, set molar to 0
            molar = 0
        molar = i*(10**-9) # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value', axis = 1)

    return x

df_final = pIC50(df_combined)
df_final

#Describe the pIC50 column
df_final.pIC50.describe()

"""Chemical Space Analysis"""

#Import Library for Data Expoloration
import seaborn as sns
sns.set(style='ticks')
import matplotlib.pyplot as plt

#Plot graph for bioactivity class
plt.figure(figsize=(5.5, 6.5))

sns.countplot(x='class', data=df_final, edgecolor='black', hue = 'class')

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('Frequency', fontsize=14, fontweight='bold')

#plt.savefig('pfht1_plot_bioactivity_class.pdf')

"""Scatter plot of MV versus LogP"""

#Scatter Plot for MV vs LogP
plt.figure(figsize=(5.5, 5.5))

sns.scatterplot(x='MW', y='LogP', data=df_final, hue='class', size='pIC50', edgecolor='black', alpha=0.7)

plt.xlabel('MW', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
#plt.savefig('plot_MW_vs_LogP.pdf')

df_final

#Remove NaN
df_cleaned = df_final.dropna()

df_cleaned

#subset the pIC50
pIC = df_final[("pIC50")]

# There might be one or more valid SMILES that can represent one compound
# Thanks to Pat Walters for this information,checkout his excellent blog: https://www.blogger.com/profile/18223198920629617711
def canonical_smiles(smiles):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    smiles = [Chem.MolToSmiles(mol) for mol in mols]
    return smiles

from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

Chem.MolFromSmiles('C=CCC')

Chem.MolFromSmiles('CCC=C')

"""Calculate descriptors using RDkit"""

def RDkit_descriptors(smiles):
    mols = [Chem.MolFromSmiles(i) for i in smiles]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
    desc_names = calc.GetDescriptorNames()

    Mol_descriptors =[]
    for mol in mols:
        # add hydrogens to molecules
        mol=Chem.AddHs(mol)
        # Calculate all 200 descriptors for each molecule
        descriptors = calc.CalcDescriptors(mol)
        Mol_descriptors.append(descriptors)
    return Mol_descriptors,desc_names

# Function call
Mol_descriptors,desc_names = RDkit_descriptors(df_cleaned["canonical_smiles"])

df_with_200_descriptors = pd.DataFrame(Mol_descriptors,columns=desc_names)
df_with_200_descriptors

fp_pIC = pd.concat([df_with_200_descriptors, df_lipinski, pIC], axis=1)

#end of video
fp_pIC.head()

# the Molecular is now in duplicate, one from Lipinski df (MW) and the other from the molecular descriptor (MolWt), you can drop one
fp_pIC = fp_pIC.drop('MW', axis=1)

fp_pIC

# Drop NaN
fp_pIC = fp_pIC.dropna()





#Library for Model
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

#define X and Y
X= fp_pIC.drop(columns=['pIC50'])
Y=fp_pIC.pIC50

#Slipt the data into trainig(80%) and testing sets(20%)
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

# prompt: remove the column Ipc from X_train data
if 'Ipc' in X_train.columns:
  X_train = X_train.drop('Ipc', axis=1)

# prompt: remove column Ipc from X_test, Y_train and Y_test

if 'Ipc' in X_test.columns:
  X_test = X_test.drop('Ipc', axis=1)

# Y_train and Y_test are series of pIC50 values, they don't have columns to drop
# If you meant to drop rows containing NaN in Y_train and Y_test you can use
# Y_train = Y_train.dropna()
# Y_test = Y_test.dropna()

#Print the shapes of the training and testing sets to confirm
print(f"Training set size: X_train: {X_train.shape}, Y_train: {Y_train.shape}")
print(f"Testing set size: X_test: {X_test.shape}, Y_test: {Y_test.shape}")

Model= RandomForestRegressor(n_estimators=100, random_state=42)
Model.fit(X_train, Y_train)

from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

# Predictions on the test set (effects of the model on the data)
y_pred = Model.predict(X_test)

mse = mean_squared_error(Y_test, y_pred)
mae = mean_absolute_error(Y_test, y_pred)
r2 = r2_score(Y_test, y_pred)

print(f"MSE: {mse}, MAE: {mae}, R-squared: {r2}")

"""Cross-validation (k-fold cross validation,k=5 or k=10) to ensure the model generalises well to unseen data

"""

from sklearn.model_selection import cross_val_score

if 'Ipc' in X_train.columns:
  X_train = X_train.drop('Ipc', axis=1)

if 'Ipc' in X_test.columns:
  X_test = X_test.drop('Ipc', axis=1)

# prompt: remove Ipc from X and Y and Model

if 'Ipc' in X.columns:
  X = X.drop('Ipc', axis=1)

# Y is a series, so we don't need to remove a column.
# If you meant to remove rows with NaN in Y, you can use:
# Y = Y.dropna()

# Refit the model without Ipc
Model = RandomForestRegressor(n_estimators=100, random_state=42)
Model.fit(X_train, Y_train)

# prompt: cross validation of model

cv_scores = cross_val_score(Model, X, Y, cv=5)  # 5-fold cross-validation
print(f"Cross-validation scores: {cv_scores}")
print(f"Mean cross-validation score: {np.mean(cv_scores)}")
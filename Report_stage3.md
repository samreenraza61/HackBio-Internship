# Docking and Pipeline Implementation

**Authors (@slack)**:  
- Samreen Raza (@samRaza)  
- Lwethu Twana (@Lwethu23)  
- Tanvi Thakur (@Jerry)

---

## Phase 1

### Introduction

Ovarian cancer has become one of the main causes of death in women globally, primarily due to its late-stage diagnosis and constrained treatment choices. Although there have been advances in medication, long-term survival rates have not increased considerably. Recent studies emphasize the role of epigenetic regulation, such as histone changes, in cancer growth. Epigenetic alterations, like histone modification and degradation, have a significant impact on gene expression and have been identified as potential targets for therapy. 

Histone Deacetylase 4 (HDAC4) deacetylates lysine residues on core histones (H2A, H2B, H3, and H4) to control transcription, cell cycle progression, and maturation. It forms multiprotein complexes and interacts with MEF2 factors to help in muscle maturation. HDAC4 epigenetically controls ESR1 in breast cancer and deacetylates HSPA1A/B proteins [1]. Elevated levels of HDAC4 in malignancies such as ovarian and non-small cell lung cancer (NSCLC) are connected to active tumor characteristics and poor prognosis, often resulting in decreased survival [2].

Using PubChem [3], we selected 50 phytochemicals from _Camellia sinensis_, a widely consumed beverage recognized for its cancer prevention benefits [4].

---

### Methodology

#### Protein Preparation

- Download the 6fyz structure from PDB [5].
- Process it in Discovery Studio [6] by eliminating heteroatoms and water molecules.
- Integrate polar hydrogens and save in PDB format.

#### Library Creation

- Curate 50 phytochemicals using PubChem.
- Convert them into PDBQT format.

#### Docking Setup

- Load the protein and ligands from the library in PyRx [7].
- Select key residues for active/binding sites.
- Set up a grid box.

#### Docking Process

- Assess the interactions between the receptor protein (6fyz) and the phytochemical library.
- PrankWeb [8] identified additional pockets that aid in understanding binding interactions, detailed in the supplementary file.

| Parameter            | Value                           |
|----------------------|---------------------------------|
| Active site          | 803                             |
| Binding sites        | 667, 669, 675, 751              |
| Receptor             | 6fyz.pdbqt                      |
| Exhaustiveness       | 8                               |
| Grid Box Center (X,Y,Z) | -22.8026, 21.8574, -8.2483    |
| Grid Box Dimensions (X,Y,Z) | 52.8871, 67.2849, 25     |

**Table 1. Docking Parameters and Active Site Information.**

#### Molecular Docking

- Dock the phytochemical library with 6FYZ in PyRx (or AutoDock Vina), noting binding affinities and interactions.
  
#### Results
- Compile all docking-generated tables and images.

## Results

Includes visuals of protein active sites, grid selection, top two ligand-protein interactions, and other input/output files. A phytochemical library listing 50 compounds' structures and action mechanisms is also provided. All these resources are accessible via our GitHub repository for further review. 


| Ligand                               | Binding Affinity | RMSD/ub | RMSD/lb |
|--------------------------------------|------------------|---------|---------|
| 6fyz_5280899_uff_E=704.15 (Zeaxanthin) | -9.2             | 0       | 0       |
| 6fyz_5280899_uff_E=704.15 (Zeaxanthin) | -9.0             | 1.596   | 0.595   |
| 6fyz_5281243_uff_E=655.10 (Leutin)     | -8.9             | 0       | 0       |

**Table 2. Top three ligands based on binding affinity.**

---

### Conclusion

Our docking analysis found that Zeaxanthin and Leutin are the most significant phytochemicals with high binding affinities to 6fyz, particularly at active sites. Zeaxanthin showed the best potential (-9.2 kcal/mol) as an HDAC4 inhibitor in ovarian cancer.


### References

1. Shen YF, Wei AM, Kou Q, Zhu QY, Zhang L. Histone deacetylase 4 increases progressive epithelial ovarian cancer cells via repression of p21 on fibrillar collagen matrices. *Oncol Rep*. 2016 Feb;35(2):948-54. doi: 10.3892/or.2015.4423. PMID: 26572940.
2. Ahn MY, Kang DO, Na YJ, et al. Histone deacetylase inhibitor, apicidin, inhibits human ovarian cancer cell migration via class II histone deacetylase 4 silencing. *Cancer Lett*. 2012 Dec 28;325(2):189-99. doi: 10.1016/j.canlet.2012.06.017. PMID: 22781396.
3. Kim S, Chen J, Cheng T, et al. PubChem 2023 update. *Nucleic Acids Research*, 2023 Jan;51(D1):D1373–D1380. doi: 10.1093/nar/gkac956.
4. Wang Y, Ren N, Rankin GO, Li B, et al. Anti-proliferative effect and cell cycle arrest induced by saponins extracted from tea (Camellia sinensis) flower in human ovarian cancer cells. *J Funct Foods*. 2017;37:310-321. doi: 10.1016/j.jff.2017.08.001.
5. H.M. Berman, J. Westbrook, Z. Feng, et al. The Protein Data Bank (2000) *Nucleic Acids Res*. 28(1):235-242. doi: 10.1093/nar/28.1.235.
6. Discovery Studio, version 24.1.0.23298, Accelrys, 2024.
7. Dallakyan S, Olson AJ. Small-Molecule Library Screening by Docking with PyRx. *Methods Mol Biol.* 2015;1263:243-50. doi: 10.1007/978-1-4939-2269-7_19.
8. Jendele L, Krivak R, Skoda P, Novotny M, Hoksza D. PrankWeb: a web server for ligand binding site prediction and visualization. *Nucleic Acids Res*. 2019 Jul 2;47(W1):W345-W349. doi: 10.1093/nar/gkz424.

---

## Phase 2

### Introduction

This phase involves the use of machine learning and cheminformatics libraries to analyze bioactivity data for the target "Histone Deacetylase 4" (HDAC4), a known therapeutic target for cancer and other diseases.

---
## Methods

### Imports from Libraries

The following libraries were utilized:

- **chembl-webresource-client**: This library was used to obtain bioactivity information from the ChEMBL database.
- **RDKit**: A cheminformatics toolkit, specifically used for the creation and management of molecular descriptors.
- **Mordred**: A molecular descriptor generation tool that computes various molecular features.

---

### Target Search

The script searches for and retrieves bioactivity data for "Histone Deacetylase 4" using the `chembl-webresource-client` library. This allows for fetching relevant bioactivity information, essential for building the dataset.

---

### Descriptor Calculation

The script computes molecular descriptors for each compound using **RDKit** and **Mordred**. These descriptors represent various molecular features and will likely be used as input features for machine learning models.

---

### Bioactivity Data Retrieval

After determining the target (HDAC4), the script retrieves bioactivity data. This data will be used to build a dataset that connects molecular descriptors (features) to bioactivity outcomes (target).

---

### Preprocessing and Model Development

The script contains code to prepare molecular data for machine learning models. This involves several preprocessing techniques, including molecular descriptor calculation and data cleaning, ensuring the dataset is ready for training.

---

### Machine Learning Phases

The machine learning workflow consists of two primary stages:

- **Training**: A portion of the dataset is used to train models that leverage molecular descriptors to predict bioactivity.
- **Testing**: The remaining dataset is used to evaluate the performance of the trained models, ensuring they generalize well to unseen data.

---

### 3. Results

#### Codes and Output of Preprocessing Steps:

The Preprocessing steps begin with searching for the target protein Histone Deacetylase 4 on the ChEMBL database.  
Assign the first entry which corresponds to the target protein to the `selected_target` variable.  
To find the inhibitory concentration, or IC50 values, the bioactivity data for the target "Histone Deacetylase 4" is filtered using the code:

```python
# filter bioactivity data using IC50 column
data = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50") 
```

To guarantee the originality of compounds, duplicate entries are removed using SMILES by running the code:

```python
df2_nr = df2.drop_duplicates(['canonical_smiles'])
```

The data goes through a number of preprocessing stages, such as:

- Standard_value (bioactivity measure) and canonical_smiles (chemical structure) missing value removal:

```python
# Data pre-processing of the bioactivity data
df2 = df[df.standard_value.notna()]
df2 = df2[df.canonical_smiles.notna()]
```

- Combine the 3 columns and bioactivity_class into a DataFrame:

```python
selection = ['molecule_chembl_id','canonical_smiles','standard_value']
df3 = df2_nr[selection]
```

Compounds are classified based on IC50 values into three categories:

- Inactive: IC50 ≥ 10,000 nM.
- Active: IC50 ≤ 1,000 nM.
- Intermediate: 1,000 nM < IC50 < 10,000 nM.

A new column class is added to the dataset to represent these categories using the following code:

```python
# Addition of class column 
bioactivity_threshold = []
for i in df3.standard_value:
  if float(i) >= 10000:
    bioactivity_threshold.append("inactive")
  elif float(i) <= 1000:
    bioactivity_threshold.append("active")
  else:
    bioactivity_threshold.append("intermediate")

bioactivity_class = pd.Series(bioactivity_threshold, name='class')
df4 = pd.concat([df3, bioactivity_class], axis=1)
```

The script below calculates molecular descriptors based on Lipinski’s rule of five, which is used to assess the drug-likeness of a compound:

```python
# Inspired by: https://codeocean.com/explore/capsules?query=tag

:data-curation
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
        row = np.array([desc_MolWt, desc_MolLogP, desc_NumHDonors, desc_NumHAcceptors])
        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1
    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)
    return descriptors
```

Run Lipinski using canonical_smiles:

```python
df_lipinski = lipinski(df4.canonical_smiles)
```

Combine Bioactivity data table and Lipinski values:

```python
# Combine df4 table with lipinski value table
df_combined = pd.concat([df4,df_lipinski], axis=1)
```

Convert IC50 to pIC50: This is to allow the data to be more uniformly distributed. pIC50 is calculated from the previous standard_value column:

```python
# Import library
import numpy as np

# Convert IC50 to pIC50
def pIC50(input):
    pIC50 = []
    for i in input['standard_value']:
        try:
            i = float(i)
        except ValueError:
            molar = 0
        molar = i*(10**-9)  # Converts nM to M
        pIC50.append(-np.log10(molar))
    input['pIC50'] = pIC50
    x = input.drop('standard_value', axis = 1)
    return x
```

Describe the `pIC50` column:

```python
df_final.pIC50.describe()
```

Exploratory Data Analysis via Lipinski descriptors:

```python
# Import Library for Data Exploration
import seaborn as sns
sns.set(style='ticks')
import matplotlib.pyplot as plt

# Plot graph for bioactivity class
plt.figure(figsize=(5.5, 6.5))
sns.countplot(x='class', data=df_final, edgecolor='black', hue = 'class')
plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('Frequency', fontsize=14, fontweight='bold')

# Scatter Plot for MV vs LogP
plt.figure(figsize=(5.5, 5.5))
sns.scatterplot(x='MW', y='LogP', data=df_final, hue='class', size='pIC50', edgecolor='black', alpha=0.7)
plt.xlabel('MW', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
```

Calculate descriptors using RDKit:

```python
# Calculate descriptor using RDKit
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
```

Function call:

```python
Mol_descriptors,desc_names = RDkit_descriptors(df_cleaned["canonical_smiles"])
```

Calculate all 200 descriptors for each molecule:

```python
df_with_200_descriptors = pd.DataFrame(Mol_descriptors,columns=desc_names)
df_with_200_descriptors
```

### Training and Testing Phase

Import library for Random Forest Regressor Model:

```python
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
```

Define `X` and `Y`:

```python
X= fp_pIC.drop(columns=['pIC50'])
Y=fp_pIC.pIC50
```

Split the data into training (80%) and testing sets (20%):

```python
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)
```

Remove the column `Ipc` from `X_train` data:

```python
if 'Ipc' in X_train.columns:
  X_train = X_train.drop('Ipc', axis=1)
```

Train Model to Data:

```python
Model= RandomForestRegressor(n_estimators=100, random_state=42)
Model.fit(X_train, Y_train)
```

### Model Evaluation

After training, we will evaluate the model using MSE, MAE, and R-squared. Also, we will perform cross-validation to assess generalization:

```python
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

# Predictions on the test set
y_pred = Model.predict(X_test)

# Calculate evaluation metrics
mse = mean_squared_error(Y_test, y_pred)
mae = mean_absolute_error(Y_test, y_pred)
r2 = r2_score(Y_test, y_pred)
print(f"MSE: {mse}, MAE: {mae}, R-squared: {r2}")
```

Output:

```plaintext
MSE: 1.434448392089781, MAE: 0.9173858303675583, R-squared: 0.2674666489747807
```

### Cross Validation

```python
from sklearn.model_selection import cross_val_score

# Cross-validation of model
cv_scores = cross_val_score(Model, X, Y, cv=5)  # 5-fold cross-validation
print(f"Cross-validation scores: {cv_scores}")
print(f"Mean cross-validation score: {np.mean(cv_scores)}")
```

Output:

```plaintext
Cross-validation scores: [-0.52565608 -0.15397229 -0.51313287 -3.09446009 -0.24337113]
Mean cross-validation score: -0.906118493479025
```

---

### Importance of Features and Bioactivity Influence

LogP and Molecular Weight were quite predictive, probably because they affect how well a chemical interacts with the target protein. For example, LogP value of 1.4296 indicates that the compound with chemical id `CHEMBL343448` is expected to show better bioactivity as they can balance between solubility and membrane permeability, which is further validated by the bioactivity_class being active.

Better model performance is indicated by a lower MSE. In this instance, there is an average squared error of 1.434 between the model's predictions and the actual results. The mean absolute error (MAE) between the expected and actual values indicates that, on average, the forecasts are 0.917 units off from the genuine values. With an R2 value of 0.267, the model is able to explain 26.7% of the variability in the data, indicating a weak fit.

When utilizing the `scoring='neg_mean_squared_error'` option (as in the previous code), cross-validation typically yields negative values for MSE, which is why negative MSE is utilized. Improved performance is indicated by a smaller value that is closer to zero. A negative value nearer zero denotes a smaller prediction error because the mean square error is often positive.

In terms of fit quality, a score near 0 denotes a better fit. A mean score of -0.906 in this instance indicates that the model is functioning reasonably, if not flawlessly, as it is not too far from zero. We can note a variance between folds which indicates that the model may not generalize equally well across all subsets of our data, as evidenced by the variability between the folds, particularly the -3.09 score.

---

### References

1. Mendez D, Gaulton A, Bento AP, Chambers J, et al. ChEMBL: towards direct deposition of bioassay data. *Nucleic Acids Res*. 2019 Jan;47(D1):D930–D940. doi: 10.1093/nar/gky1075.
2. Landrum, G., et al. RDKit documentation. Available at: https://www.rdkit.org/docs/.



### **Supplementary_files**:  
[Suplementary_files_phase1&phase2_both](https://github.com/samreenraza61/HackBio-Internship/tree/main/Supplementary_files_stage3)


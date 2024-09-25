# Ovarian Cancer HDAC4 Docking and Machine Learning Analysis

**Authors (@slack): Samreen Raza (@samRaza), Lwethu Twana (@Lwethu23), Tanvi Thakur (@Jerry)**

---

## Phase 1: Introduction

Ovarian cancer has become one of the leading causes of death in women globally, primarily due to its late-stage diagnosis and limited treatment options. Although there have been advances in medication, long-term survival rates have not improved significantly. Recent research highlights the role of epigenetic regulation, such as histone modifications, in cancer progression. Epigenetic alterations, like histone modification and degradation, have a profound impact on gene expression and are considered potential therapeutic targets.

Histone Deacetylase 4 (HDAC4) deacetylates lysine residues on core histones (H2A, H2B, H3, and H4) to regulate transcription, cell cycle progression, and differentiation. It forms multiprotein complexes and interacts with MEF2 factors, facilitating muscle differentiation. HDAC4 epigenetically regulates ESR1 in breast cancer and deacetylates HSPA1A/B proteins [1]. Elevated levels of HDAC4 in malignancies, such as ovarian and non-small cell lung cancer (NSCLC), are associated with aggressive tumor traits and a poor prognosis, indicating shorter survival rates [2].

Using PubChem [3], we selected 50 phytochemicals from *Camellia sinensis*, a widely consumed beverage known for its anti-cancer properties [4].

---

## Methodology

### Protein Preparation
The 6fyz structure was downloaded from PDB [5], processed in Discovery Studio [6] by removing heteroatoms and water molecules, adding polar hydrogens, and saved in PDB format.

### Library Creation
We curated 50 phytochemicals using PubChem and converted them into PDBQT format.

### Docking Setup
Both the protein and ligand library were loaded into PyRx [7], key residues were selected for active/binding sites, and a grid box was configured.

### Docking Process
Interactions between the receptor protein (6fyz) and the phytochemical library were evaluated. PrankWeb [8] predicted additional pockets, enhancing the docking analysis. Detailed information is provided in the supplementary file.

| **Parameter**       | **Value**                             |
|---------------------|---------------------------------------|
| Active site         | 803                                   |
| Binding sites       | 667,669,675,751                       |
| Receptor            | 6fyz.pdbqt                            |
| Exhaustiveness      | 8                                     |
| Grid Box Center     | -22.8026, 21.8574, -8.2483            |
| Grid Box Dimensions | 52.8871, 67.2849, 25                  |

*Table 1. Docking Parameters and Active Site Information.*

---

## Results

![](figure1.png)
**Figure 1:** Chain A active site residues highlighted in yellow.

| **Ligand**                                   | **Binding Affinity** | **RMSD/ub** | **RMSD/lb** |
|----------------------------------------------|----------------------|-------------|-------------|
| 6fyz_5280899_uff_E=704.15 (Zeaxanthin)       | -9.2                 | 0           | 0           |
| 6fyz_5280899_uff_E=704.15 (Zeaxanthin)       | -9.0                 | 1.596       | 0.595       |
| 6fyz_5281243_uff_E=655.10 (Leutin)           | -8.9                 | 0           | 0           |

*Table 2. Top three ligands based on binding affinity.*

---

## Conclusion

Our docking analysis identified Zeaxanthin and Leutin as the most promising phytochemicals with high binding affinities to HDAC4, particularly at active sites. Zeaxanthin showed the highest potential (-9.2 kcal/mol) as an HDAC4 inhibitor in ovarian cancer.

---

**Supplementaryfiles**:  
[Suplementary_files_phase1_docking_results](https://github.com/samreenraza61/HackBio-Internship/tree/main/Supplementary_files_stage3)

## References

1. Shen YF, Wei AM, Kou Q, Zhu QY, Zhang L. Histone deacetylase 4 increases progressive epithelial ovarian cancer cells via repression of p21 on fibrillar collagen matrices. Oncol Rep. 2016 Feb;35(2):948-54. doi: 10.3892/or.2015.4423. Epub 2015 Nov 16. PMID: 26572940.
2. Ahn MY, Kang DO, Na YJ, Yoon S, Choi WS, Kang KW, Chung HY, Jung JH, Min do S, Kim HS. Histone deacetylase inhibitor, apicidin, inhibits human ovarian cancer cell migration via class II histone deacetylase 4 silencing. Cancer Lett. 2012 Dec 28;325(2):189-99. doi: 10.1016/j.canlet.2012.06.017. Epub 2012 Jul 7. PMID: 22781396.
3. Kim, S., Chen, J., Cheng, T., Gindulyte, A., He, J., He, S., Li, Q., Shoemaker, B. A., Thiessen, P. A., Yu, B., Zaslavsky, L., Zhang, J., & Bolton, E. E. (2023). PubChem 2023 update. Nucleic acids research, 51(D1), D1373–D1380. https://doi.org/10.1093/nar/gkac956
4. Wang, Y., Ren, N., Rankin, G. O., Li, B., Rojanasakul, Y., Tu, Y., & Chen, Y. C. (2017). Anti-proliferative effect and cell cycle arrest induced by saponins extracted from tea (Camellia sinensis) flower in human ovarian cancer cells. Journal of functional foods, 37, 310–321. https://doi.org/10.1016/j.jff.2017.08.001
5. H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne, The Protein Data Bank (2000) Nucleic Acids Research 28: 235-242 https://doi.org/10.1093/nar/28.1.235.
6. D. Studio, Discovery Studio, version 24.1.0.23298, Accelrys, 2024, p. 420.
7. Dallakyan S, Olson AJ. Small-Molecule Library Screening by Docking with PyRx. Methods Mol Biol. 2015;1263:243-50. https://www.researchgate.net/publication/273954875_Small-Molecule_Library_Screening_by_Docking_with_PyRx
8. Jendele L, Krivak R, Skoda P, Novotny M, Hoksza D. PrankWeb: a web server for ligand binding site prediction and visualization. Nucleic Acids Res. 2019 Jul 2;47(W1):W345-W349. doi: 10.1093/nar/gkz424. PMID: 31114880; PMCID: PMC6602436.

---

# Phase 2

## Introduction
The report involves the use of machine learning and cheminformatics libraries to analyze bioactivity data for the target **"Histone Deacetylase 4" (HDAC4)**, a known therapeutic target for cancer and other diseases.

## Methods

### Imports from Libraries
- **chemical-webresource-client**: Used to obtain bioactivity information from the ChEMBL database.
- **RDKit**: A cheminformatics toolkit used for the creation and management of molecular descriptors.
- **Mordred**: A molecular descriptor generation tool.

### Target Search
The script looks for and retrieves the bioactivity data for "Histone Deacetylase 4" using `chembl-webresource-client`.

### Descriptor Calculation
The script computes each compound's molecular descriptors using **RDKit** and **Mordred**; these features will likely be utilized in machine learning.

### Bioactivity Data Retrieval
The script retrieves bioactivity data after determining the target, most likely to build a dataset that connects molecular characteristics to bioactivity.

### Preprocessing and Model Development
The script contains code to prepare molecular data for machine learning models by preprocessing it using techniques like descriptor creation.

### Machine Learning Phases
Consists of a training and testing stage for creating models that use molecular descriptors to predict bioactivity.

## Results

### Codes and Output of Preprocessing Steps

1. **Target Search**
   The preprocessing steps begin with searching for the target protein **Histone Deacetylase 4** on the ChEMBL database.

2. **Assign Selected Target**
   Assign the first entry corresponding to the target protein to the `selected_target` variable.

3. **Filter by IC50**
   To find the inhibitory concentration (IC50 values), the bioactivity data for the target **"Histone Deacetylase 4"** is filtered using the code:

   ```python
   # Filter bioactivity data using IC50 column
   data = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50") 
4. **Remove Duplicate Entries**
 To guarantee the originality of compounds, duplicate entries are removed using SMILES by running the code:

 ```python
df2_nr = df2.drop_duplicates(['canonical_smiles'])

### **5. Preprocessing Missing Values**
The data goes through a number of preprocessing stages, such as removing missing values for standard_value (bioactivity measure) and canonical_smiles (chemical structure):

```python
# Data pre-processing of the bioactivity data
df2 = df[df.standard_value.notna()]
df2 = df2[df.canonical_smiles.notna()]

6. **Bioactivity Classification**
 Compounds are classified based on IC50 values into three categories:

Inactive: IC50 ≥ 10,000 nM.
Active: IC50 ≤ 1,000 nM.
Intermediate: 1,000 nM < IC50 < 10,000 nM.
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
7. **Lipinski’s Rule of Five**
The script calculates molecular descriptors based on Lipinski’s Rule of Five, which is used to assess the drug-likeness of a compound:

python
def lipinski(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        if isinstance(elem, str):
            mol = Chem.MolFromSmiles(elem)
            moldata.append(mol)
    baseData = np.arange(1,1)
    i = 0
    for mol in moldata:
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
        row = np.array([desc_MolWt, desc_MolLogP, desc_NumHDonors, desc_NumHAcceptors])
        if i == 0:
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i += 1
    columnNames = ["MW", "LogP", "NumHDonors", "NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)
    return descriptors

8. **Combine Bioactivity Data and Lipinski Values**
The bioactivity data table is combined with Lipinski descriptor values:

python
# Combine df4 table with lipinski value table
df_combined = pd.concat([df4, df_lipinski], axis=1)

9. **Convert IC50 to pIC50**
 To allow the data to be more uniformly distributed, pIC50 is calculated from the previous standard_value column:

python
def pIC50(input):
    pIC50 = []
    for i in input['standard_value']:
        try:
            i = float(i)
        except ValueError:
            molar = 0
        molar = i * (10**-9)  # Converts nM to M
        pIC50.append(-np.log10(molar))
    input['pIC50'] = pIC50
    x = input.drop('standard_value', axis=1)
    return x
10. **Exploratory Data Analysis via Lipinski Descriptors**
The script performs exploratory data analysis via visualization of Lipinski descriptors:

python
import seaborn as sns
sns.set(style='ticks')
import matplotlib.pyplot as plt
# Plot graph for bioactivity class
plt.figure(figsize=(5.5, 6.5))
sns.countplot(x='class', data=df_final, edgecolor='black', hue='class')
And also generates a scatter plot for MW vs LogP:

python
sns.scatterplot(x='MW', y='LogP', data=df_final, hue='class', size='pIC50', edgecolor='black', alpha=0.7)
11. **Calculate Descriptors using RDKit**
 RDKit is used to calculate 200 molecular descriptors:

python
def RDkit_descriptors(smiles):
    mols = [Chem.MolFromSmiles(i) for i in smiles]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
    desc_names = calc.GetDescriptorNames()
    Mol_descriptors = []
    for mol in mols:
        mol = Chem.AddHs(mol)
        descriptors = calc.CalcDescriptors(mol)
        Mol_descriptors.append(descriptors)
    return Mol_descriptors, desc_names
12. **Training and Testing Phase**
 The Random Forest Regressor Model is used for training and testing:

python
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
# Define X and Y
X = fp_pIC.drop(columns=['pIC50'])
Y = fp_pIC.pIC50
# Split data into training(80%) and testing sets(20%)
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)
If Ipc column exists, it is removed:

python
if 'Ipc' in X_train.columns:
    X_train = X_train.drop('Ipc', axis=1)
13. **Model Evaluation**
 After training, the model is evaluated using MSE, MAE, and R-squared. Cross-validation is also performed:

python
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
y_pred = Model.predict(X_test)
mse = mean_squared_error(Y_test, y_pred)
mae = mean_absolute_error(Y_test, y_pred)
r2 = r2_score(Y_test, y_pred)
print(f"MSE: {mse}, MAE: {mae}, R-squared: {r2}")
Results:

MSE: 1.434
MAE: 0.917
R-squared: 0.267
Cross-Validation Cross-validation is performed using 5-fold cross-validation:

python
from sklearn.model_selection import cross_val_score
cv_scores = cross_val_score(Model, X, Y, cv=5)
print(f"Cross-validation scores: {cv_scores}")
print(f"Mean cross-validation score: {np.mean(cv_scores)}")
Results:

Cross-validation scores: [-0.52565608, -0.15397229, -0.51313287, -3.09446009, -0.24337113]
Mean cross-validation score: -0.906518.
Conclusion

Importance of features and bioactivity influence 

LogP and Molecular Weight were quite predictive, probably because they affect how well a chemical interacts with the target protein.For example, LogP value of 1,4296 indicates that the compound with chemical id :CHEMBL343448 is expected to show better bioactivity as they can balance between solubility and membrane permeability which is further validates by the bioactivity_class being active.

Better model performance is indicated by a lower MSE. In this instance, there is an average squared error of 1.434 between the model's predictions and the actual results.The mean absolute error (MAE) between the expected and actual values indicates that, on average, the forecasts are 0.917 units off from the genuine values.With an R2 value of 0.267, the model is able to explain 26.7% of the variability in the data, indicating a pretty weak fit.

When utilizing the scoring='neg_mean_squared_error' option (as in the previous code), cross-validation typically yields negative values for MSE, which is why negative MSE is utilized.
Improved performance is indicated by a smaller value that is closer to zero. A negative value nearer zero denotes a smaller prediction error because the mean square error is often positive.
In terms of Fit quality, a score near 0 denotes a better fit. A mean score of -0.906 in this instance indicates that the model is functioning reasonably, if not flawlessly, as it is not too far from zero.
We can note a variance between folds which indicates that the model may not generalize equally well across all subsets of our data, as evidenced by the variability between the folds, particularly the -3.09 score.

**Supplementaryfiles**:  
[Suplementary_files_Phase2_task](https://github.com/samreenraza61/HackBio-Internship/tree/main/Supplementary_files_stage3/phase2_stage3)

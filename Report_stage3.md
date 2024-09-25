**Docking and Pipeline Implementation**

Authors (@slack): Samreen Raza (@samRaza), Lwethu Twana (@Lwethu23),
Tanvi Thakur (@Jerry)

**Phase 1**

1.  **Introduction**

Ovarian cancer has become one of the main causes of death in women
globally, primarily due to its late-stage diagnosis and constrained
treatment choices. Although advances in medication, long-term survival
rates have not increased considerably. The latest study emphasizes the
role of epigenetic regulation, such as histone changes, in cancer
growth. Epigenetic alterations, such as histone modification and
degradation, have an important impact on gene expression and have been
identified as possible targets for therapy.

Histone Deacetylase 4 (HDAC4) deacetylates lysine residues on core
histones (H2A, H2B, H3, and H4) to control transcription, cell cycle
progression, and maturation. It creates multiprotein complexes and acts
with MEF2 factors to help in muscle maturity. HDAC4 epigenetically
controls ESR1 in breast cancer and deacetylates HSPA1A/B proteins
**\[1\]**. Increased levels of HDAC4 in malignancies such as ovarian and
non-small cell lung cancer (NSCLC) are connected to active tumor
characteristics as a negative prognosis, among those exhibiting
decreased longevity **\[2\]**.

Using PubChem **\[3\]** we selected 50 phytochemicals from Camellia
sinensis, a widely used beverage recognized for its cancer prevention
benefits **\[4\]**.

2.  **Methodology**

```{=html}
<!-- -->
```
1.  Protein preparation

-   Download the 6fyz structure from PDB **\[5\]**, process it in
    > Discovery Studio **\[6\]** by eliminating heteroatoms and water
    > molecules, subsequently integrating polar hydrogens, and save it
    > in PDB format

2.  Library creation

-   Curate 50 phytochemicals using PubChem, convert them into PDBQT
    > format

3.  Docking setup

-   Load the protein and ligands from the library in PyRx **\[7\]**,
    > select key residues for active/binding sites, and set up a grid
    > box.

4.  Docking process

-   Assess the interactions within the receptor protein (6fyz) and the
    > phytochemical library, including the docking parameters and
    > settings outlined in Table 1.

-   PrankWeb **\[8\]** discovered more pockets that help to comprehend
    > binding interactions. Details are included in the supplementary
    > file, which improves the docking analysis.

  -----------------------------------------------------------------------
  **Parameter**                         **Value**
  ------------------------------------- ---------------------------------
  Active site                           803

  Binding sites                         667,669,675,751

  Receptor                              6fyz.pdbqt

  Exhaustiveness                        8

  Grid Box Center (X, Y, Z)             -22.8026, 21.8574, -8.2483

  Grid Box Dimensions (X, Y, Z)         52.8871,67.2849,25
  -----------------------------------------------------------------------

Table 1. Docking Parameters and Active Site Information.

5.  Molecular docking

-   Dock the phytochemical library with **6FYZ** in **PyRx** (or
    > **AutoDock Vina**), noting binding affinities and interactions in
    > Table 2.

6.  Results

-   Compile all docking-generated tables and images.

3.  **Results**![](vertopal_9dab5013fede455188ac4efdf1e9d875/media/image5.png){width="4.364583333333333in"
    > height="2.1035258092738407in"}

Figure 1: Chain A active site residues highlighted in yellow.

  -------------------------------------------------------------------------------
  **Ligand**                         **Binding        **RMSD/ub**   **RMSD/lb**
                                     Affinity**                     
  ---------------------------------- ---------------- ------------- -------------
  6fyz_5280899_uff_E=704.15          -9.2             0             0
  (Zeaxanthin)                                                      

  6fyz_5280899_uff_E=704.15          -9.0             1.596         0.595
  (Zeaxanthin)                                                      

  6fyz_5281243_uff_E=655.10 (Leutin) -8.9             0             0
  -------------------------------------------------------------------------------

Table 2. Top three ligands, based on binding affinity.

4.  **Conclusion**

Our docking analysis found Zeaxanthin and Leutin are the most
significant phytochemicals with high binding affinities to 6fyz,
particularly at active sites. Zeaxanthin showed the best potential (-9.2
kcal/mol) as an HDAC4 inhibitor in ovarian cancer.

5.  **References**

```{=html}
<!-- -->
```
1.  Shen YF, Wei AM, Kou Q, Zhu QY, Zhang L. Histone deacetylase 4
    > increases progressive epithelial ovarian cancer cells via
    > repression of p21 on fibrillar collagen matrices. Oncol Rep. 2016
    > Feb;35(2):948-54. doi: 10.3892/or.2015.4423. Epub 2015 Nov 16.
    > PMID: 26572940.

2.  Ahn MY, Kang DO, Na YJ, Yoon S, Choi WS, Kang KW, Chung HY, Jung JH,
    > Min do S, Kim HS. Histone deacetylase inhibitor, apicidin,
    > inhibits human ovarian cancer cell migration via class II histone
    > deacetylase 4 silencing. Cancer Lett. 2012 Dec 28;325(2):189-99.
    > doi: 10.1016/j.canlet.2012.06.017. Epub 2012 Jul 7. PMID:
    > 22781396..

3.  Kim, S., Chen, J., Cheng, T., Gindulyte, A., He, J., He, S., Li, Q.,
    > Shoemaker, B. A., Thiessen, P. A., Yu, B., Zaslavsky, L., Zhang,
    > J., & Bolton, E. E. (2023). PubChem 2023 update. *Nucleic acids
    > research*, *51*(D1), D1373--D1380.
    > [[https://doi.org/10.1093/nar/gkac956]{.underline}](https://doi.org/10.1093/nar/gkac956)

4.  Wang, Y., Ren, N., Rankin, G. O., Li, B., Rojanasakul, Y., Tu, Y., &
    > Chen, Y. C. (2017). Anti-proliferative effect and cell cycle
    > arrest induced by saponins extracted from tea (*Camellia
    > sinensis*) flower in human ovarian cancer cells. *Journal of
    > functional foods*, *37*, 310--321.
    > [[https://doi.org/10.1016/j.jff.2017.08.001]{.underline}](https://doi.org/10.1016/j.jff.2017.08.001)

5.  H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H.
    > Weissig, I.N. Shindyalov, P.E. Bourne, The Protein Data
    > Bank (2000) *Nucleic Acids Research* 28: 235-242
    > <https://doi.org/10.1093/nar/28.1.235>.

6.  D. Studio, Discovery Studio, version 24.1.0.23298, Accelrys,
    > 2024, p. 420.

7.  [Small-Molecule Library Screening by Docking with
    > PyRx.](http://www.ncbi.nlm.nih.gov/pubmed/25618350) Dallakyan S,
    > Olson AJ. *Methods Mol Biol.* 2015;1263:243-50. The full-text is
    > available at
    > <https://www.researchgate.net/publication/273954875_Small-Molecule_Library_Screening_by_Docking_with_PyRx>

8.  Jendele L, Krivak R, Skoda P, Novotny M, Hoksza D. PrankWeb: a web
    > server for ligand binding site prediction and visualization.
    > Nucleic Acids Res. 2019 Jul 2;47(W1):W345-W349. doi:
    > 10.1093/nar/gkz424. PMID: 31114880; PMCID: PMC6602436.

**Phase 2**

1.  **Introduction**

The report involves the use of machine learning and cheminformatics
libraries to analyze bioactivity data for the target \"Histone
Deacetylase 4\" (HDAC4), a known therapeutic target for cancer and other
diseases.

2.  **Methods**

```{=html}
<!-- -->
```
1.  Imports from Libraries: Utilized to obtain bioactivity information
    > from the ChEMBL database is the chemical-web resource-client.

> Rdkit is used as a cheminformatics toolkit, specifically for the
> creation and management of molecular descriptors.
>
> mordred: A molecular descriptor generation tool.

2.  Target Search: The script looks for and retrieves the bioactivity
    > data for \"Histone Deacetylase 4\" using
    > chembl-webresource-client.

3.  Descriptor Calculation: The script computes each compound\'s
    > molecular descriptors using RDKit and Mordred; these features will
    > probably be utilized in machine learning.

4.  Bioactivity Data Retrieval: The script retrieves bioactivity data
    > after determining the target, most likely to build a dataset that
    > connects molecular characteristics to bioactivity.

5.  Preprocessing and Model Development: The script contains code to
    > prepare molecular data for machine learning models by
    > preprocessing it using techniques like descriptor creation.

6.  Machine Learning phases:consists of a training and testing stage for
    > creating models that use molecular descriptors to predict
    > bioactivity.

**3. Results**

1.  Codes and Output of Preprocessing Steps:

The Preprocessing steps begin with searching for the target protein
Histone Deacetylase 4 on the chembl database

Assign the first entry which corresponds to the target protein to the
selected_target variable

To find the inhibitory concentration, or IC50 values, the bioactivity
data for the target \"Histone Deacetylase 4\" is filtered using the
code:

#filter bioactivity data using IC50 column

data =
activity.filter(target_chembl_id=selected_target).filter(standard_type=\"IC50\")

To guarantee the originality of compounds, duplicate entries are removed
using SMILES by running the code:df2_nr =
df2.drop_duplicates(\[\'canonical_smiles\'\])

The data goes through a number of preprocessing stages, such as:
Standard_value (bioactivity measure) and canonical_smiles (chemical
structure) missing value removal

#Data pre-processing of the bioactivity data

df2 = df\[df.standard_value.notna()\]

df2 = df2\[df.canonical_smiles.notna()\]

#Combine the 3 columns and bioactivity_class into a DataFrame

selection =
\[\'molecule_chembl_id\',\'canonical_smiles\',\'standard_value\'\]

df3 = df2_nr\[selection\]

Compounds are classified based on IC50 values into three categories:

-   **Inactive**: IC50 ≥ 10,000 nM.

-   **Active**: IC50 ≤ 1,000 nM.

-   **Intermediate**: 1,000 nM \< IC50 \< 10,000 nM.

A new column class is added to the dataset to represent these categories
using the following codes:

#Addition of class column

bioactivity_threshold = \[\]

for i in df3.standard_value:

if float(i) \>= 10000:

bioactivity_threshold.append(\"inactive\")

elif float(i) \<= 1000:

bioactivity_threshold.append(\"active\")

else:

bioactivity_threshold.append(\"intermediate\")

bioactivity_class = pd.Series(bioactivity_threshold, name=\'class\')

df4 = pd.concat(\[df3, bioactivity_class\], axis=1)

The script below calculates molecular descriptors based on **Lipinski's
rule of five**, which is used to assess the drug-likeness of a compound:

\# Inspired by:
https://codeocean.com/explore/capsules?query=tag:data-curation

\# Addition of Lipinski Descriptors

def lipinski(smiles, verbose=False):

moldata= \[\]

for elem in smiles:

\# Check if elem is a valid SMILES string before converting

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

row = np.array(\[desc_MolWt,

desc_MolLogP,

desc_NumHDonors,

desc_NumHAcceptors\])

if(i==0):

baseData=row

else:

baseData=np.vstack(\[baseData, row\])

i=i+1

columnNames=\[\"MW\",\"LogP\",\"NumHDonors\",\"NumHAcceptors\"\]

descriptors = pd.DataFrame(data=baseData,columns=columnNames)

return descriptors

#Run lipinski using canonical_smiles

df_lipinski =
lipinski(df4.canonical_smiles)![](vertopal_9dab5013fede455188ac4efdf1e9d875/media/image2.png){width="4.84375in"
height="2.6875in"}

Combine Bioactivity data table and Lipinski values:

#Combine df4 table with lipinski value table

df_combined = pd.concat(\[df4,df_lipinski\], axis=1)

Convert IC50 to pIC50:This is to allow the data to be more uniformly
distributed, pIC50 is calculated from the previous standard\_ value
column

#Import library

import numpy as np

#convert IC50 to pIC50

def pIC50(input):

pIC50 = \[\]

for i in input\[\'standard_value\'\]:

try: \# this will try to convert i to a float

i = float(i)

except ValueError: \# if i is not a number, set molar to 0

molar = 0

molar = i\*(10\*\*-9) \# Converts nM to M

pIC50.append(-np.log10(molar))

input\[\'pIC50\'\] = pIC50

x = input.drop(\'standard_value\', axis = 1)

return x

#Describe the pIC50 column

df_final.pIC50.describe()

![](vertopal_9dab5013fede455188ac4efdf1e9d875/media/image1.png){width="2.2604166666666665in"
height="3.34375in"}

Exploratory Data Analysis via lipinski descriptors

#Import Library for Data Expoloration

import seaborn as sns

sns.set(style=\'ticks\')

import matplotlib.pyplot as plt

#Plot graph for bioactivity class

plt.figure(figsize=(5.5, 6.5))

sns.countplot(x=\'class\', data=df_final, edgecolor=\'black\', hue =
\'class\')

plt.xlabel(\'Bioactivity class\', fontsize=14, fontweight=\'bold\')

plt.ylabel(\'Frequency\', fontsize=14, fontweight=\'bold\')

#Scatter Plot for MV vs LogP

plt.figure(figsize=(5.5, 5.5))

sns.scatterplot(x=\'MW\', y=\'LogP\', data=df_final, hue=\'class\',
size=\'pIC50\', edgecolor=\'black\', alpha=0.7)

plt.xlabel(\'MW\', fontsize=14, fontweight=\'bold\')

plt.ylabel(\'LogP\', fontsize=14, fontweight=\'bold\')

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)

![](vertopal_9dab5013fede455188ac4efdf1e9d875/media/image4.png){width="4.5537970253718285in"
height="4.129529746281714in"}

![](vertopal_9dab5013fede455188ac4efdf1e9d875/media/image6.png){width="5.447916666666667in"
height="3.8678083989501313in"}

Calculate descriptors using RDkit:

\# Calculate descriptor using RDkit

def RDkit_descriptors(smiles):

mols = \[Chem.MolFromSmiles(i) for i in smiles\]

calc = MoleculeDescriptors.MolecularDescriptorCalculator(\[x\[0\] for x
in Descriptors.\_descList\])

desc_names = calc.GetDescriptorNames()

Mol_descriptors =\[\]

for mol in mols:

\# add hydrogens to molecules

mol=Chem.AddHs(mol)

\# Calculate all 200 descriptors for each molecule

descriptors = calc.CalcDescriptors(mol)

Mol_descriptors.append(descriptors)

return Mol_descriptors,desc_names

\# Function call

Mol_descriptors,desc_names =
RDkit_descriptors(df_cleaned\[\"canonical_smiles\"\])

#Calculate all 200 descriptors for each molecule

df_with_200_descriptors =
pd.DataFrame(Mol_descriptors,columns=desc_names)

df_with_200_descriptors

2.  Training and Testing Phase

Import library for Random Forest Regressor Model:

from sklearn.ensemble import RandomForestRegressor

from sklearn.model_selection import train_test_split

#define X and Y

X= fp_pIC.drop(columns=\[\'pIC50\'\])

Y=fp_pIC.pIC50

#Split the data into training(80%) and testing sets(20%)

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2,
random_state=42)

\# prompt: remove the column Ipc from X_train data

if \'Ipc\' in X_train.columns:

X_train = X_train.drop(\'Ipc\', axis=1)

The removal of the Ipc table was due to the presence of infinity, which
resulted in an error code when run in the model.

#Train Model to Data

Model= RandomForestRegressor(n_estimators=100, random_state=42)

Model.fit(X_train, Y_train)

![](vertopal_9dab5013fede455188ac4efdf1e9d875/media/image3.png){width="3.2708333333333335in"
height="0.6666666666666666in"}

Model Evaluation :After training, we will evaluate the model using
**MSE**, **MAE**, and **R-squared**. Also, we will perform
**cross-validation** to assess generalization.

from sklearn.metrics import mean_squared_error, mean_absolute_error,
r2_score

\# Predictions on the test set

y_pred = Model.predict(X_test)

#Calculate evaluation metrics

mse = mean_squared_error(Y_test, y_pred)

mae = mean_absolute_error(Y_test, y_pred)

r2 = r2_score(Y_test, y_pred)

print(f\"MSE: {mse}, MAE: {mae}, R-squared: {r2}\")

MSE: 1.434448392089781, MAE: 0.9173858303675583, R-squared:
0.2674666489747807

Cross validation:

from sklearn.model_selection import cross_val_score

\# prompt: cross validation of model

cv_scores = cross_val_score(Model, X, Y, cv=5) \# 5-fold
cross-validation

print(f\"Cross-validation scores: {cv_scores}\")

print(f\"Mean cross-validation score: {np.mean(cv_scores)}\")

Cross-validation scores: \[-0.52565608 -0.15397229 -0.51313287
-3.09446009 -0.24337113\]

Mean cross-validation score: -0.906118493479025

3.  Importance of features and bioactivity influence

LogP and Molecular Weight were quite predictive, probably because they
affect how well a chemical interacts with the target protein.For
example, LogP value of 1,4296 indicates that the compound with chemical
id :CHEMBL343448 is expected to show better bioactivity as they can
balance between solubility and membrane permeability which is further
validates by the bioactivity_class being active.

Better model performance is indicated by a lower MSE. In this instance,
there is an average squared error of 1.434 between the model\'s
predictions and the actual results.The mean absolute error (MAE) between
the expected and actual values indicates that, on average, the forecasts
are 0.917 units off from the genuine values.With an R2 value of 0.267,
the model is able to explain 26.7% of the variability in the data,
indicating a pretty weak fit.

When utilizing the scoring=\'neg_mean_squared_error\' option (as in the
previous code), cross-validation typically yields negative values for
MSE, which is why negative MSE is utilized.

Improved performance is indicated by a smaller value that is closer to
zero. A negative value nearer zero denotes a smaller prediction error
because the mean square error is often positive.

In terms of Fit quality, a score near 0 denotes a better fit. A mean
score of -0.906 in this instance indicates that the model is functioning
reasonably, if not flawlessly, as it is not too far from zero.

We can note a variance between folds which indicates that the model may
not generalize equally well across all subsets of our data, as evidenced
by the variability between the folds, particularly the -3.09 score.

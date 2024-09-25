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

## Phase 2: Machine Learning and Cheminformatics

### Introduction

This phase involves the use of machine learning and cheminformatics libraries to analyze bioactivity data for the target Histone Deacetylase 4 (HDAC4), a known therapeutic target for cancer and other diseases.

---

## Methods

### Library Imports
- **chembl-webresource-client**: Used to retrieve bioactivity information from the ChEMBL database.
- **RDKit**: Cheminformatics toolkit for molecular descriptor creation.
- **Mordred**: A molecular descriptor generation tool.

### Target Search
The target search retrieves bioactivity data for "Histone Deacetylase 4" from ChEMBL.

### Descriptor Calculation
Using RDKit and Mordred, the script computes molecular descriptors, which will be used in machine learning.

---

## Results

### Code Snippets for Preprocessing

- **Bioactivity Data Retrieval**:
    ```python
    data = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
    ```

- **Removing Duplicates**:
    ```python
    df2_nr = df2.drop_duplicates(['canonical_smiles'])
    ```

- **Classifying IC50 Values**:
    ```python
    bioactivity_threshold = []
    for i in df3.standard_value:
        if float(i) >= 10000:
            bioactivity_threshold.append("inactive")
        elif float(i) <= 1000:
            bioactivity_threshold.append("active")
        else:
            bioactivity_threshold.append("intermediate")
    ```

- **Calculating Lipinski Descriptors**:
    ```python
    def lipinski(smiles):
        moldata= []
        for elem in smiles:
            mol=Chem.MolFromSmiles(elem)
            moldata.append(mol)
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
    ```

---

## Conclusion

This phase has integrated molecular descriptors from bioactivity data into the predictive machine learning models, enhancing the understanding of HDAC4 inhibitors' properties.


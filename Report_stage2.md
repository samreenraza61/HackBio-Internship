# Report

**Authors (@slack):** Samreen Raza (@samRaza), Tanvi Thakur (@Jerry)  
**GitHub:** [https://github.com/samreenraza61/HackBio-Internship/blob/main/Report_stage2.md)

## Introduction

Cellular tumor antigen p53 is a transcription factor that responds to DNA damage by halting the cell cycle, allowing time for repair. If repair fails, p53 triggers apoptosis to prevent damaged DNA from being passed on, reducing cancer risk [1]. The human p53 protein consists of 393 amino acids with key domains: transcription activation domains (TAD1: 1-40, TAD2: 41-60), proline-rich region (61-93), DNA-binding domain (102-293), tetramerization domain (323-353), and C-terminal domain (364-393) [2]. Mutations, particularly in key amino acids like R175, R248, and R273 [3]-[4], are found in about 50% of human cancers [5]-[6]. This makes p53 a major focus in cancer research for potential therapies [1].

## Methods

SWISS-MODEL simplifies the process with integrated databases and software, making model construction and validation easier [7]. Template selection, essential for accuracy, was aided by hidden Markov models to identify structurally similar templates [8]. After evaluation, the most suitable template was chosen as shown in the table below.

| Template  | Method      | Coverage | Seq Identity | QMEAN | GMQUE |
|-----------|-------------|----------|--------------|-------|-------|
| 7xzz.1.K  | EM, 4.07Å   | 1.00     | 100          | 0.93  | 0.54  |

AlphaFold2 integrates deep learning into protein structures in three phases: creating MSAs, refining them using Evoformer, and identifying 3D atom positions [9,10]. Our results provide a PAE plot, which indicates general reliability but increased errors in residues 100-200 and 250-350. These locations may require additional validation or signal flexibility.

## Results

- **Modeling:** Swiss-Model findings reveal high quality: 98.47% of residues in the favored region, and the MolProbity score is 0.62. AlphaFold2 has strong confidence in its main structure but less so in residues, indicating flexibility.
- **Domain:** SWISS-MODEL predicts only the DNA-binding domain, while AlphaFold2 predicts all domains.
- **Visualization:** Alignment of predicted structures by PyMOL [11] shows a low RMSD of 0.390 Å for 172 atoms, indicating high structural similarity between the two predictions.
- **Structural Accuracy:**
  - For the apo form, the homology model has an RMSD of 0.271 Å (high similarity) vs. AlphaFold2 (24.494 Å).
  - For the agonist-bound, both methods show high similarity, with 0.323 Å (homology) and 0.340 Å (AlphaFold2).
  - For the antagonist-bound, the Swiss model shows less variation (6.223 Å) than AlphaFold2 (10.855 Å).

## Conclusion

Swiss-Model, which uses verified templates, is reliable for widely studied proteins, providing more accuracy. AlphaFold, on the other hand, employs deep learning and is effective when no templates are available but not as precise. In this study, Homology modeling is more effective for p53, resulting in reduced RMSD as well as greater structural similarity.

## Supplementary Files
[Supplementary Files](https://github.com/samreenraza61/HackBio-Internship/tree/main/Supplementaryfiles-Stage2)

## References

[1] Y. Huang et al., “An overview of the functions of p53 and drugs acting either on wild- or mutant-type p53,” European Journal of Medicinal Chemistry, vol. 265, p. 116121, Feb. 2024, doi: 10.1016/j.ejmech.2024.116121.

[2] A. C. Joerger and A. R. Fersht, “Structural biology of the tumor suppressor p53,” Annu Rev Biochem, vol. 77, pp. 557–582, 2008, doi: 10.1146/annurev.biochem.77.060806.091238.

[3] W. A. Freed-Pastor and C. Prives, “Mutant p53: one name, many proteins,” Genes Dev, vol. 26, no. 12, pp. 1268–1286, Jun. 2012, doi: 10.1101/gad.190678.112.

[4] X. Yue, Y. Zhao, Y. Xu, M. Zheng, Z. Feng, and W. Hu, “Mutant p53 in Cancer: Accumulation, Gain-of-Function, and Therapy,” J Mol Biol, vol. 429, no. 11, pp. 1595–1606, Jun. 2017, doi: 10.1016/j.jmb.2017.03.030.

[5] N. Rivlin, R. Brosh, M. Oren, and V. Rotter, “Mutations in the p53 Tumor Suppressor Gene: Important Milestones at the Various Steps of Tumorigenesis,” Genes Cancer, vol. 2, no. 4, pp. 466–474, Apr. 2011, doi: 10.1177/1947601911408889.

[6] L. Bouaoun et al., “TP53 Variations in Human Cancers: New Lessons from the IARC TP53 Database and Genomics Data,” Hum Mutat, vol. 37, no. 9, pp. 865–876, Sep. 2016, doi: 10.1002/humu.23035.

[7] L. Bordoli, F. Kiefer, K. Arnold, P. Benkert, J. Battey, and T. Schwede, “Protein structure homology modeling using SWISS-MODEL workspace,” Nat Protoc, vol. 4, no. 1, pp. 1–13, 2009, doi: 10.1038/nprot.2008.197.

[8] Dalton, J. A., & Jackson, R. M. (2007). An evaluation of automated homology modelling methods at low target–template sequence similarity. Bioinformatics, 23(15), 1901-1908.

[9] P. Cramer, “AlphaFold2 and the future of structural biology,” Nat Struct Mol Biol, vol. 28, no. 9, pp. 704–705, Sep. 2021, doi: 10.1038/s41594-021-00650-1.

[10] Z. Yang, X. Zeng, Y. Zhao, and R. Chen, “AlphaFold2 and its applications in the fields of biology and medicine,” Sig Transduct Target Ther, vol. 8, no. 1, pp. 1–14, Mar. 2023, doi: 10.1038/s41392-023-01381-z.

[11] The PyMOL Molecular Graphics System, "PyMOL," [Online]. Available: https://www.pymol.org/. [Accessed: Sep. 15, 2024].


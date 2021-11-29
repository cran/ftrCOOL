# News

## ftrCOOL 2.0.0

### Added the following functions:

#### Amino Acid Functions
- AESNN3
- ASDC
- binary_3bit_T1
- binary_3bit_T2
- binary_3bit_T3
- binary_3bit_T4
- binary_3bit_T5
- binary_3bit_T6
- binary_3bit_T7
- binary_5bit_T1
- binary_5bit_T2
- binary_6bit
- DistancePair
- OPF_7bit_T1
- OPF_7bit_T2
- OPF_7bit_T3
- OPF_10bit

#### DNA/RNA Functions

- ASDC_DNA
- ASDC_RNA
- DPCP_DNA
- DPCP_RNA
- KNN_DNA
- KNN_RNA
- Mismatch_DNA
- Mismatch_RNA
- MMI_DNA
- MMI_RNA
- PS2_DNA
- PS2_RNA
- PS3_DNA
- PS3_RNA
- PS4_DNA
- PS4_RNA
- TPCP_DNA
- TPCP_RNA
- Zcurve9bit_DNA
- Zcurve9bit_RNA
- Zcurve12bit_DNA
- Zcurve12bit_RNA
- Zcurve36bit_DNA
- Zcurve36bit_RNA
- Zcurve48bit_DNA
- Zcurve48bit_RNA
- Zcurve144bit_DNA
- Zcurve144bit_RNA
- Zcurve9bit_DNA

### Fixed a bug in 'upto' parameter of CkSNUCpair_DNA and CkSAApair 


## ftrCOOL 1.1.1
### Fixed a bug in the example of function DisorderB


## ftrCOOL 1.1.0

### Fixed bugs in the following functions:
 - AAutoCor()
 - AutoCorDiNUC()
 - AutoCorTriNUC()


### Renamed the following functions:
 - AutoCorDiNUC to AutoCorDiNUC_DNA
 - AutoCorTriNUC to AutoCorTriNUC_DNA
 - CkSNUCpair to CkSNUCpair_DNA
 - CodonUsage to CodonUsage_DNA
 - DiNUC2Binary to DiNUC2Binary_DNA
 - DiNUCindex to DiNUCindex_DNA
 - ENUComposition to ENUComposition_DNA
 - ExpectedValKmerNUC to ExpectedValKmerNUC_DNA
 - G_Ccontent to G_Ccontent_DNA
 - kNUComposition to kNUComposition_DNA
 - LocalPoSpKaaF to LocalPoSpKAAF
 - LocalPoSpKnucF to LocalPoSpKNUCF_DNA
 - maxORFlength to  maxORFlength_DNA
 - NUC2Binary to NUC2Binary_DNA
 - NUCKpartComposition to NUCKpartComposition_DNA
 - PSEkNUCdi to PSEkNUCdi_DNA
 - PSEkNUCTri to PSEkNUCTri_DNA
 - TriNucIndex to TriNucIndex_DNA


### Added the following functions:
 - ANF_RNA()
 - ASA()
 - APkNUCdi_RNA()
 - AutoCorDiNUC_RNA()
 - CkSNUCpair_RNA()
 - CodonUsage_RNA()
 - DiNUC2Binary_RNA()
 - DiNUCindex_RNA()
 - DisorderB()
 - DisorderC()
 - DisorderS()
 - ENUComposition_RNA()
 - ExpectedValKmerNUC_RNA()
 - G_Ccontent_RNA()
 - KNNPeptide()
 - KNNProtein()
 - kNUComposition_RNA()
 - LocalPoSpKNUCF_RNA()
 - maxORF_RNA()
 - maxORFlength_RNA()
 - NCP_DNA()
 - NCP_RNA()
 - NUC2Binary_RNA()
 - NUCKpartComposition_RNA()
 - PCPseDNC()
 - PSEkNUCdi_RNA()
 - PSSM()
 - PSTNPds()
 - PSTNPss_DNA()
 - PSTNPss_RNA()
 - SSEB()
 - SSEC()
 - SSES()
 - TorsionAngle()

### Substituted functions:
- ANF_DNA() with function novel_PseKNC()

### Improving the execution time of most functions 


## ftrCOOL 0.1.0

Initial Release

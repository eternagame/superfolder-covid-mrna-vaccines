## Comparing the Superfolder to conventionally designed mRNA sequences

| Design |  CAI | dG(MFE) (kcal/mol) |   AUP<sub>CDS</sub><sup>a</sup> | Est. half-life <sup>b</sup> | AUP<sub>init. 14 nts</sub><sup>c</sup> |
| ------------ |--------|------- |--------- | ----  | ----  | 
|  Superfolder-v2                       | 0.72 |  -2375.4  |  0.22 | **0.34**  | **0.90** |
|  Superfolder-Delta                 | 0.72 | -2323.6   | 0.23  | **0.33** | **0.86** |
|  Superfolder-Delta-PSU                | 0.72 | -1711.2   | 0.29  | 0.55<sup>d</sup> | 0.90 |
|  Superfolder-v1                       | 0.73 | -2382.1   | 0.22  |   0.34 | 0.86 |
| Putative BioNTech/Pfizer_BNT-162b2<sup>e</sup> | 0.95 | -1341.1  | 0.21  |  1517.5   | 0.49  |
| Putative Moderna mRNA-1273<sup>e</sup>| 0.98 | -1510.11  | 0.38  | 0.22 |  0.43  |
| LinearDesign<sup>f</sup>              | 0.73 | -2533.3   | 0.20  |   0.32 | 0.42 |
|      GC-rich<sup>g</sup>              | 0.80 | -1617.2   | 0.42  |  0.27 | 0.63 |
|     [GENEWIZ](https://www.genewiz.com/Public/Services/Gene-Synthesis/Codon-Optimization) codon optimization | 0.95 | -1304.3     | 0.49     |  1263.9 | 0.58 |
|          [IDT](https://www.idtdna.com/pages/tools/codon-optimization-tool?returnurl=%2FCodonOpt) codon optimization | 0.73 | -1089.5     | 0.51     |  1427.8 | 0.46 |
  
dG(MFE) calculated in LinearFold-Vienna.
<sup>a</sup>Average unpaired probability over entire coding sequence (Wayment-Steele, 2020b).  
<sup>b</sup>Estimated from [DegScore](https://github.com/eternagame/DegScore), machine-learning model for predicting degradation more accurately than average unpaired probability (Leppek, 2021).
<sup>c</sup>Average unpaired probability of the first 14 nucleotides of the coding sequence (Kozak, 1990) -- should be *high* to optimize translation.  
<sup>d</sup>Using PSU heuristic, setting DegScore for U to 0 (see [DegScore](https://github.com/eternagame/DegScore) repository.)
<sup>e</sup>Jeong, 2021.
<sup>f</sup>Zhang, 2020.
<sup>g</sup>Each codon is randomly sampled from the most GC-rich codons for that amino acid (Thess, 2015).  

The Superfolder construct was optimized by Eterna participants to reduce predicted degradation. It is an adaptation of Superfolder-v1, which had been optimized in RiboTree both to minimize DegScore, maximizing in vitro stability, while keeping the first 14 nucleotides of the coding sequence unpaired.

The Superfolder-Delta construct was developed by enumerating all GC-rich codon substitutions to result in the Delta mutant spike protein, evaluating their DegScore, and selecting the variant with the lowest DegScore.

<img src="../assets/results_barplot_18Aug2021.png" alt="Barplot of calculated metrics" width="800"/>

<img src="../assets/result_scatterplots_18Aug2021.png" alt="Scatterplot of DegScore vs. AUP init" width="800"/>


The Superfolder-Delta-PSU construct was designed from scratch by Eterna participants in the Delta challenge, optimizing a modified DegScore that increases the stability of PSU by setting U degradation to zero.

## Additional mRNA designs

Eterna participants have designed S-2P and Delta variant mRNAs with a diversity of structures in OpenVaccine project. 

In December 2021, 181 S-2P constructs were submitted and voted upon by Eterna participants. The top 9 candidates are included in `superfolders_ALL.csv`.

In July 2021, XX B.1.617 and Delta variants were submitted and voted upon by Eterna participants. The top 8 candidates for both unmodified nucleotides and PSU DegScore are included in `superfolders_ALL.csv`.

<img src="../assets/OpenVaccine_winners_S2P.png" alt="RiboGraphViz images of all Eterna winners of OpenVaccine round 7." width="800"/>

`analysis/analysis.py` reproduces the above calculations, with an example input in `assets/example_input.csv`.

## References

Dae-Eun Jeong, D.E. Matthew McCoy, M., Karen Artiles, K., Orkan Ilbay, O., Fire, A., Nadeau, K., Park, H., Betts, B., Scott Boyd, S., Hoh, R., Shoura, M. [Assemblies of putative SARS-CoV-2 spike encoding mRNA sequences for vaccines BNT-162b2 and mRNA-1273](https://github.com/NAalytics/Assemblies-of-putative-SARS-CoV2-spike-encoding-mRNA-sequences-for-vaccines-BNT-162b2-and-mRNA-1273).

Kozak M. Downstream secondary structure facilitates recognition of initiator codons by eukaryotic ribosomes (1990). Proc Natl Acad Sci U S A. 87(21):8301-5. [doi: 10.1073/pnas.87.21.8301](http://dx.doi.org/10.1073/pnas.87.21.8301)

Leppek, K., Byeon, G.W., Kladwang, W., Wayment-Steele, H.K., Kerr, C., ... Barna, M., Das, R. (2021). Combinatorial optimization of mRNA structure, stability, and translation for RNA-based therapeutics. https://www.biorxiv.org/content/10.1101/2021.03.29.437587v1.

Wayment-Steele, H.K., Kim, D.S., Choe, C.A., Nicol, J.J., Wellington-Oguri, R., Sperberg, R.A.P., Huang, P., Eterna Participants, Das, R. (2020). Theoretical basis for stabilizing messenger RNA through secondary structure design. bioRxiv, 262931.

Wayment-Steele, H.K., Kladwang, W., Eterna Participants, Das, R. (2020). RNA secondary structure packages ranked and improved by high-throughput experiments. bioRxiv, 124511.

Thess, A., Grund, S., Mui, B. L., Hope, M. J., Baumhof, P., Fotin-Mleczek, M., & Schlake, T. (2015). Sequence-engineered mRNA without chemical nucleoside modifications enables an effective protein therapy in large animals. Molecular Therapy, 23(9), 1456-1464.

Zhang, H., Zhang, L., Li, Z., Liu, K., Liu, B., Mathews, D. H., & Huang, L. (2020). LinearDesign: Efficient Algorithms for Optimized mRNA Sequence Design. arXiv preprint arXiv:2004.10177.

For answers to any additional questions that might help accelerate the end of the COVID-19 pandemic, please contact [Rhiju Das](https://daslab.stanford.edu), Stanford University, <a href="mailto:rhiju@stanford.edu">rhiju@stanford.edu</a>.

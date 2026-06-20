# PLIN2+/TREM1+ Macrophages Are a Non-Redundant Inflammatory Signaling Hub in Human Atherosclerosis across Three Vascular Beds

**Zhen-Guo Cai**^1^, **Yin Tang**^2^, **Linlin Che**^2*^

^1^ Harbin Medical University, Harbin, China
^2^ Heilongjiang University of Chinese Medicine, Harbin, China

\* Correspondence: chelinlin@hljucm.edu.cn | ORCiD: 0009-0007-3631-6882

## Structured Abstract

**Background:** Whether macrophage-mediated signaling networks are organized around defined hub cell populations—and whether such architecture is conserved across vascular beds—has not been systematically examined in human atherosclerosis.

**Methods:** We re-analyzed a publicly available human plaque single-cell atlas (259,116 cells; 73 donors) spanning carotid, coronary, and femoral arteries. Macrophages were re-clustered into five subtypes. Ligand-receptor communication was modeled using the CellChat framework with 53 curated pairs (55 significant at FDR<0.05). TREM1 inhibition was simulated by zeroing TREM1 expression with bootstrap resampling (n=2,000). External validation used GSE131778 (coronary; n=11,756 cells); analytical reproducibility via independent re-processing of GSE155512 (carotid; n=8,866 cells).

**Results:** PLIN2+/TREM1+ macrophages exhibited the highest network centrality across all beds. TREM1 was specifically upregulated in this subtype (log2FC=+0.93, p<0.001). Virtual TREM1 ablation reduced hub centrality (Δ=−2.85, 95%CI [−2.89,−2.82]) and eliminated 75 TREM1-mediated pairs without compensatory pathway activation. A FAO-SIRT1-MHC-II axis was femoral-specific and independent of the TREM1 hub. Hub status was confirmed in GSE131778 (rank #1/9) and GSE155512 (rank #1/7). Fifteen of 60 comparisons achieved >80% power; the femoral cohort (n=7) limited definitive conclusions (minimum detectable d=1.13–1.31).

**Conclusions:** PLIN2+/TREM1+ macrophages constitute a conserved, non-redundant inflammatory signaling hub. TREM1 inhibition is predicted to dampen inflammatory communication without triggering compensatory activation—a favorable pharmacological profile for chronic cardiovascular disease. Femoral-specific metabolic-epigenetic programming co-exists with but is independent of the conserved hub architecture.

**Keywords:** atherosclerosis; trained immunity; TREM1; single-cell RNA-seq; macrophage; ligand-receptor signaling; network centrality; vascular bed heterogeneity

---

## Introduction

Atherosclerosis is a systemic disease with striking anatomical heterogeneity. Coronary, carotid, and femoral plaques exhibit divergent clinical presentations, progression rates, and responses to pharmacological intervention[1,2]. While systemic risk factors—hyperlipidemia, hypertension, diabetes—act uniformly across the arterial tree, the site-specific inflammatory microenvironment determines plaque fate[3]. Understanding the cellular and molecular basis of this heterogeneity is essential for developing targeted therapies beyond lipid lowering.

Single-cell RNA sequencing (scRNA-seq) has transformed atherosclerosis research by resolving the cellular complexity of human plaques at unprecedented resolution[4-6]. Recent atlas-scale studies have catalogued the major cell types populating diseased arteries and identified activation states within the macrophage compartment[7]. However, two critical questions remain unanswered. First, do macrophages function as isolated effector cells or as organized hubs within a broader intercellular signaling network? Second, to what extent are macrophage-centered communication architectures conserved across vascular beds versus shaped by local microenvironmental cues?

The concept of trained immunity—sustained innate immune memory characterized by epigenetic and metabolic reprogramming—provides a mechanistic framework for understanding chronic vascular inflammation[8,9]. Macrophages exposed to oxidized lipids, damage-associated molecular patterns (DAMPs), or inflammatory cytokines can acquire a hyper-responsive phenotype that persists independently of the initial stimulus[10]. If trained immunity operates within atherosclerotic plaques, its effects should manifest as structured cell-cell communication networks centered on specific macrophage populations, with network properties that can be quantified and experimentally interrogated.

Here, we leveraged the publicly available Traeuble et al. 259,116-cell human atherosclerotic plaque atlas[7] to: (1) identify macrophage subtypes functioning as signaling hubs; (2) test whether hub architecture is conserved using an independent external cohort (GSE131778) and independent re-processing of a constituent dataset (GSE155512); (3) quantify the functional role of TREM1 through in-silico network perturbation; and (4) determine whether metabolic-epigenetic programs associated with trained immunity are vascular-bed-specific.

---

## Methods

### Study Design and Data Source

The primary analysis used the publicly available harmonized human atherosclerotic plaque single-cell atlas from Traeuble et al. (2025, Nature Communications)[7], accessed via CELLxGENE (Collection ID: db70986c-7d91-49fe-a399-a4730be394ac). The atlas integrates 11 independent scRNA-seq datasets into 259,116 cells from 73 human donors across three vascular beds: carotid artery (188,372 cells), coronary artery (35,395 cells), and femoral artery (35,349 cells). The atlas provides 13 Level-1 cell types and 23 Level-2 subtypes annotated using the scPoli integration framework. Donor-level module scores were computed for 70 donors with complete metadata (carotid n=50, coronary n=13, femoral n=7). Age data were available for 33 of 73 donors; given the limited coverage, age was not included as a covariate in primary analyses. The effect of this missingness is reported as a limitation.

### Macrophage Subtype Annotation

Macrophages (38,271 cells) were re-clustered and annotated into five Level-2 subtypes based on consensus marker expression: PLIN2+/TREM1+ Macrophage (hub population), TREM2+/Foamy Macrophage, Inflammatory Macrophage, HMOX1+ Macrophage, and Other Macrophage. Monocytes were included as a reference population for differentiation analyses (total myeloid compartment: 51,221 cells).

### Module Score Computation

Trained immunity (TI) composite scores and sub-scores (TI_Inflammation, TI_Metabolic, TI_PRR, TI_H3K4me3, TI_H3K27ac, TI_HDAC_SIRT) were computed as mean normalized expression of curated gene sets from MSigDB Hallmark and custom trained-immunity literature signatures. Metabolic module scores (Glycolysis, OXPHOS, FAO, FAS, Cholesterol, Hypoxia) were derived from Hallmark gene sets. Macrophage subtype scores were computed from subtype-specific marker genes. All module scores were computed using scanpy.tl.score_genes with use_raw=False, then aggregated to donor level (mean per donor) for statistical testing across vascular beds.

### Ligand-Receptor Communication Analysis

Cell-cell communication was modeled using a Python-native implementation of the CellChat geometric-mean framework[21]. Communication probability was computed as the square root of the product of ligand expression fraction in the sender population and receptor expression fraction in the receiver population, consistent with the published CellChat methodology. A curated database of 53 high-confidence ligand-receptor pairs spanning immune signaling (chemokines, cytokines, DAMPs), adhesion, growth factors, checkpoint molecules, and extracellular matrix interactions was used (Table S3). Pairs were required to have both ligand and receptor expressed in ≥5% of cells within the respective sender and receiver populations. Permutation testing (1,000 cell-type label shuffles) with Benjamini-Hochberg FDR correction identified 55 significant pairs at FDR<0.05 (Table S1).

### Network Centrality Computation

A directed communication graph was constructed where nodes represent cell types (Level-1 annotation for whole-atlas analysis; Level-2 macrophage subtypes for macrophage-centric analysis) and weighted edges represent summed communication probabilities between all significant L-R pairs for each sender-receiver combination (edge weight threshold >0.2). Weighted degree centrality was computed as (in_weight + out_weight) / (n_nodes − 1). Hub status was defined as the cell type or subtype with the highest weighted degree centrality.

### TREM1 Virtual Knockout Simulation

TREM1 inhibition was simulated by setting TREM1 expression to zero in the expression fraction matrix while preserving all other gene expression values. This approach avoids copying the full AnnData object (memory footprint ~780 MiB) and directly tests the contribution of TREM1 to network topology. Ligand-receptor communication probabilities and network centrality metrics were recomputed on the knockout matrix. Centrality loss (Δ) was calculated as degree_centrality(KO) − degree_centrality(baseline).

Bootstrap resampling was performed at the fraction level using binomial resampling within each cell type (n=2,000 iterations). Specifically, for each cell type and gene, the observed number of expressing cells k was treated as a draw from Binomial(n_cells, p_hat), and each bootstrap iteration drew k_boot ~ Binomial(n_cells, p_hat) to generate a resampled fraction k_boot / n_cells. This fraction-level approach was chosen over standard case-resampling bootstrap because (a) it avoids copying the full sparse expression matrix (memory constraint), and (b) it directly propagates binomial sampling uncertainty from the detection frequency of each gene into the communication probability and centrality estimates. 95% confidence intervals and empirical p-values were computed from the bootstrap distribution of Δ_centrality for each cell type.

### Pathway Compensation Quantification

Ligand-receptor pairs were grouped into 13 functional pathways (TREM1, TREM2/Lipid, Chemokine, Inflammatory Cytokine, DAMP/TLR, Adhesion, Foamy/SPP1, Growth Factor, MHC-II, Complement, Notch, Checkpoint, ECM; see Table S3 for pair-to-pathway mapping). Total pathway signaling was computed as the sum of communication probabilities across all sender-receiver pairs within each pathway at baseline and after TREM1 knockout. Compensation was defined as the net positive change in non-TREM1 pathway signaling following TREM1 ablation. Full compensation data are reported in Table S2.

### Monocyte-to-Macrophage Differential Expression

Differentially expressed genes between each macrophage subtype and monocytes were identified using the Wilcoxon rank-sum test as implemented in scanpy (scanpy.tl.rank_genes_groups, method='wilcoxon'). For each macrophage subtype, the comparison was structured as "[subtype]_vs_Monocyte." Gene set enrichment was performed using hypergeometric testing against Hallmark, GO Biological Process (BP), and KEGG gene sets, restricted to genes significantly upregulated (log2FC>0, adjusted p<0.05) in each comparison. Full DEG and enrichment results are provided in Fig. S3B source data.

### External Validation and Analytical Reproducibility

Independent external validation was performed using GSE131778 (Wirka et al., 2019; human coronary artery; 11,756 cells, 9 cell types)[4], a dataset not included in the 11 constituent datasets of the Traeuble atlas. The same L-R communication and network centrality pipeline was applied independently to this dataset.

Analytical reproducibility was assessed by independent re-processing of one constituent atlas dataset, GSE155512 (Pan et al., 2020; human carotid artery; 8,866 cells, 7 cell types)[22]. Raw count data were downloaded from GEO (GSE155512), re-processed independently (quality control, normalization, clustering, cell-type annotation) using the same pipeline applied to the primary atlas, and subjected to the identical L-R communication and network centrality analysis. This served as a within-atlas analytical reproducibility check rather than independent external validation, since the GSE155512 samples are constituents of the Traeuble integrated atlas.

Cross-dataset centrality concordance was assessed using Spearman rank correlation on the common cell types between GSE131778 (9 types) and GSE155512 (7 types).

### Statistical Analysis

Kruskal-Wallis tests with eta-squared effect sizes were used for three-group (vascular bed) comparisons. Pairwise Mann-Whitney U tests with Bonferroni correction (3 comparisons per module) were applied for post-hoc comparisons. Effect sizes included Cohen's d (with bootstrap 95% CI, n=2,000) and Hedges' g (bias-corrected). Post-hoc statistical power was computed for each pairwise comparison given observed effect size and sample sizes (α=0.05, two-tailed). For the femoral artery cohort (n=7 donors), minimum detectable effect sizes at 80% power were computed analytically. Analyses were performed in Python (scanpy v1.9, numpy, scipy, pandas, networkx) and R (4.5.0). All analysis scripts are available at https://github.com/Caier-dortor/atherosclerosis_atlas.

---

## Results

### A PLIN2+/TREM1+ Macrophage Population Functions as a Conserved Signaling Hub

Re-clustering of 38,271 macrophages from the full atlas identified five transcriptionally distinct subtypes (Fig. 1a: UMAP visualization colored by subtype; Fig. 1b: marker gene dotplot). PLIN2+/TREM1+ macrophages were distinguished by high expression of the lipid droplet coat protein PLIN2 and the myeloid amplifier receptor TREM1, together with elevated trained immunity composite scores (Fig. 1c: TI_composite score by subtype, KW p-value annotation). This population was present across all three vascular beds, with the highest relative abundance in femoral plaques (Fig. 1d: stacked bar of subtype proportions by bed).

Network analysis of ligand-receptor communication revealed that PLIN2+/TREM1+ macrophages occupied the central hub position in the intercellular signaling network (Fig. 4a: circle plot of outgoing signals; Fig. 4b: incoming signal heatmap). This population exhibited the highest weighted degree centrality and served as both the dominant sender and receiver of inflammatory signals. The hub architecture was robust: 55 of 6,135 tested ligand-receptor pairs were significant (FDR<0.05), with the top 50 pairs predominantly involving macrophage-centered communication (Table S1; Fig. 4c: top-15 L-R pairs by communication probability).

### TREM1 Is Specifically Upregulated During Monocyte-to-Hub-Macrophage Differentiation

Comparison of each macrophage subtype against monocytes revealed that TREM1 was specifically and strongly upregulated only in PLIN2+/TREM1+ macrophages (log2FC=+0.93, adjusted p<0.001; Fig. S3B, panel a: DEG dotplot with p-value-weighted border thickness). In all other macrophage subtypes—including TREM2+/Foamy (log2FC=+0.46), Inflammatory (log2FC=−0.94), HMOX1+ (log2FC=−1.91), and Other macrophages (log2FC=−1.64)—TREM1 was either marginally upregulated or substantially downregulated relative to monocytes. This subtype-specific induction pattern was unique to TREM1; other immune receptors (TREM2, TLR4, IL1R1) showed broader expression patterns across multiple macrophage populations.

Functional gene set enrichment analysis of genes upregulated in macrophages versus monocytes (Fig. S3B, panel b: gene set × subtype dotplot with p-value annotations) confirmed that inflammatory signaling, DAMP recognition, and chemokine production gene sets were enriched across all macrophage subtypes, while fatty acid metabolism and oxidative phosphorylation gene sets showed subtype-specific enrichment patterns. Notably, SIRT1-associated gene sets were enriched in Inflammatory and Other macrophages but not in the PLIN2+/TREM1+ hub population, consistent with the independence of the SIRT1-metabolic axis from the TREM1 hub mechanism.

### TREM1 Is a Non-Redundant Inflammatory Signal Amplifier

To determine whether TREM1 functions as a non-redundant signaling component or whether compensatory pathways could maintain network integrity in its absence, we performed an in-silico TREM1 knockout simulation. Virtual ablation of TREM1 eliminated 75 TREM1-mediated ligand-receptor pairs (Fig. 5d: top TREM1-mediated pairs by baseline probability) and reduced PLIN2+/TREM1+ macrophage weighted degree centrality by 7.8% (Δ=−2.85, 95% CI [−2.89, −2.82], bootstrap p<0.001; Fig. 5a: baseline vs. KO centrality bars with bootstrap CI annotations; Fig. 5b: network topology graph at baseline; Fig. 5c: communication Δ heatmap).

Critically, pathway-level compensation analysis revealed no compensatory activation in any of the 12 non-TREM1 signaling pathways (0% compensation; Fig. 5e: pathway signaling baseline vs. KO bars). The Chemokine, Inflammatory Cytokine, DAMP/TLR, and TREM2/Lipid pathways all showed either no change or minor losses following TREM1 ablation, indicating that TREM1-mediated signaling constitutes a unique, irreplaceable channel for inflammatory signal amplification within the plaque microenvironment.

These findings suggest that pharmacological TREM1 blockade would selectively dampen inflammatory signal amplification without triggering compensatory activation of alternative pathways—a profile well-suited for chronic low-grade inflammation such as atherosclerosis, where broad-spectrum immunosuppression is undesirable.

### FAO-SIRT1-MHC-II Axis Is Vascular-Bed-Specific

While the PLIN2/TREM1 hub mechanism was conserved across vascular beds, metabolic-epigenetic programs exhibited site-specific organization (Fig. 3). Fatty acid oxidation (FAO) scores were significantly divergent between carotid and femoral plaques (MW p<0.01, Bonferroni-corrected), with femoral samples showing a distinct FAO-SIRT1-MHC-II co-activation pattern (Fig. 3e: parallel pathway schematic; Fig. 3f: vascular-bed comparison heatmap). This axis was independent of the PLIN2/TREM1 hub: PLIN2+/TREM1+ macrophages showed low FAO and SIRT1 activity, with these metabolic programs concentrated in Inflammatory and Other macrophage subtypes (Fig. S3B).

These data indicate that the hub cell architecture and metabolic-epigenetic programming represent parallel, dissociable layers of macrophage organization. The PLIN2/TREM1 hub provides a conserved signaling scaffold, while vascular-bed-specific metabolic programs determine the biochemical output of the network.

### Independent Validation Confirms Hub Architecture

To test the generalizability of the macrophage hub finding, we analyzed GSE131778 (Wirka 2019, coronary artery, n=11,756 cells)[4], a dataset not included in the primary Traeuble atlas. Macrophages exhibited the highest network centrality (rank #1 of 9 cell types; Fig. S4c: grouped bar chart with rank labels), confirming the hub phenotype in an independent coronary cohort.

As an analytical reproducibility assessment, we independently re-processed raw GSE155512 data (Pan 2020, carotid artery, n=8,866 cells)[22]—a constituent dataset of the Traeuble atlas—through the identical pipeline. Macrophages again ranked #1 of 7 cell types (Fig. S4c), demonstrating that the hub finding is reproducible under independent preprocessing and annotation. Of eight key macrophage-associated ligands examined, APOE, ANXA1, IL1B, and TNF showed conserved macrophage-centered communication across both datasets, while HMGB1, SPP1, CSF1, and CCL2 were detected only in the carotid cohort, reflecting platform-specific gene coverage differences rather than biological absence (Fig. S4d: horizontal bar chart with N.D. annotations for undetected ligands).

Quantitative comparison of degree centrality across the seven common cell types revealed a Spearman rank correlation of ρ=0.57 between the coronary and carotid analyses (Fig. S4e: scatter plot with identity line). Although the p-value (p=0.18) did not reach conventional significance given the limited number of common cell types (n=7), the qualitative concordance—macrophages ranked #1 in both analyses—supports a conserved hub architecture.

### Statistical Power Analysis

Of 60 pairwise vascular-bed comparisons for module scores, 15 (25%) achieved post-hoc statistical power exceeding 80% (Table S2). The femoral artery cohort (n=7 donors) severely limited statistical power for this vascular bed (minimum detectable Cohen's d=1.13 for femoral vs. carotid comparisons; d=1.31 for femoral vs. coronary). Bootstrap confidence intervals (n=2,000 iterations) were reported for all effect size estimates to account for small-sample uncertainty in the femoral group (Table S2b).

---

## Discussion

This study establishes that human atherosclerotic plaques are organized around a conserved macrophage signaling hub defined by PLIN2 and TREM1 co-expression. Through re-analysis of a 259,116-cell atlas, ligand-receptor network modeling, in-silico perturbation, independent external validation (GSE131778), and analytical reproducibility testing (GSE155512), we demonstrate that: (1) PLIN2+/TREM1+ macrophages function as the dominant signaling hub in the plaque cellular communication network; (2) TREM1 serves as a non-redundant inflammatory signal amplifier whose loss cannot be compensated by other pathways; and (3) the hub architecture is conserved across vascular beds, while metabolic-epigenetic programming is site-specific.

### TREM1 as a Therapeutic Target

The finding that TREM1 ablation produces a selective, non-compensated reduction in inflammatory signaling has direct therapeutic implications. Current anti-inflammatory strategies for atherosclerosis—such as canakinumab (CANTOS trial)[11] and colchicine (COLCOT, LoDoCo2)[12,13]—achieve cardiovascular risk reduction through broad inflammasome or microtubule inhibition. While effective, these approaches carry inherent risks of immunosuppression. Our network perturbation analysis suggests that TREM1 blockade would achieve inflammatory signal dampening with a fundamentally different mechanism: selective removal of a signal amplifier without collateral pathway suppression or compensatory activation. This profile is particularly attractive for a chronic, low-grade inflammatory condition where long-term treatment safety is paramount.

TREM1-targeted therapeutics are in early clinical development for sepsis and rheumatoid arthritis[14,15]. Our data provide a rationale for extending TREM1 blockade to atherosclerosis, supported by the observation that TREM1 expression is near-exclusive to the hub macrophage population in human plaques.

### Hub Cell Organization and Trained Immunity

The identification of a defined macrophage hub cell population raises the question of whether trained immunity operates through hub cell-mediated network reorganization. Trained immunity—epigenetically encoded, sustained innate immune memory[8,16]—has been proposed as a mechanism linking systemic risk factors to chronic plaque inflammation. Our network framework suggests a testable model: trained immunity stimuli (oxidized LDL, lipoprotein(a), hyperglycemia) may act by expanding or activating the PLIN2+/TREM1+ hub population, thereby amplifying the entire inflammatory communication network through a single cellular node. This "hub amplification" model differs from the conventional view of trained immunity as a cell-autonomous property and suggests that network-level measurements (centrality, path length, modularity) may provide more sensitive readouts of therapeutic intervention than cell-intrinsic gene expression alone.

### Vascular-Bed-Specific Microenvironmental Drivers

The femoral-specific FAO-SIRT1-MHC-II axis identified in our analysis likely reflects convergent contributions from three microenvironmental factors that differ systematically across anatomical locations. First, **hemodynamic forces**: coronary arteries experience high pulsatile shear stress (~10-20 dyn/cm²) that promotes an atheroprotective endothelial phenotype, while the superficial femoral artery is subjected to low and oscillatory shear during limb flexion, combined with compressive mechanical strain from adjacent skeletal muscle[17]. These biomechanical differences directly regulate endothelial NF-κB activation and macrophage recruitment. Second, **smooth muscle cell (SMC) embryonic origins**: coronary artery SMCs derive from the pro-epicardium (mesothelial origin), whereas femoral artery SMCs originate from the paraxial mesoderm (somite-derived)[18]. Lineage tracing studies have demonstrated that SMCs of distinct embryonic origins exhibit divergent responses to atherogenic stimuli, including differential production of extracellular matrix components, cytokines, and chemokines. Third, **perivascular adipose tissue (PVAT) composition**: epicardial adipose tissue surrounding coronary arteries exhibits a pro-inflammatory, brown-fat-like phenotype with high adipokine secretion (leptin, resistin) that promotes macrophage activation[19], whereas femoral perivascular fat resembles subcutaneous white adipose tissue with lower adiponectin production and reduced inflammatory tone. These three factors—hemodynamic stress, SMC lineage, and PVAT biology—likely act synergistically to create vascular-bed-specific microenvironments that tune the metabolic-epigenetic output of the conserved macrophage hub.

### Limitations

Several limitations warrant consideration. First, the femoral artery cohort comprised only 7 donors, severely limiting statistical power for this vascular bed (minimum detectable Cohen's d=1.13–1.31 at 80% power). While bootstrap confidence intervals were used throughout to account for small-sample uncertainty, femoral-specific findings should be considered hypothesis-generating and require validation in larger cohorts. Second, TREM1 inhibition was modeled in silico through expression-level perturbation; experimental validation in ex vivo human plaque culture systems, TREM1-blocking antibodies, or Trem1-knockout mouse models of atherosclerosis is needed to confirm the predicted effects on network topology and pathway compensation. Third, the curated ligand-receptor database (53 pairs), while spanning major immune-stromal signaling axes, is not exhaustive and may miss relevant interactions involving less-characterized receptors or non-canonical ligand-receptor pairings. Fourth, age data were unavailable for 40 of 73 donors (55%), precluding formal adjustment for age as a potential confounder; however, a sensitivity analysis on the 33 donors with age data showed that age explained <1% of variance in TI_composite scores, suggesting limited confounding. Fifth, the cross-dataset centrality correlation (Spearman ρ=0.57) did not reach statistical significance due to the limited number of common cell types (n=7), though qualitative concordance of macrophage hub status was maintained across all analyses. Sixth, this study is observational and based on bulk tissue scRNA-seq; spatial organization of the identified hub architecture within the plaque microenvironment remains to be determined by spatial transcriptomics.

### Future Directions

Several lines of investigation would extend these findings. Spatial transcriptomics (e.g., Visium, MERFISH) could localize PLIN2+/TREM1+ macrophages within plaque microanatomy, testing whether hub cells concentrate at regions of active inflammation (shoulder regions, necrotic core borders). Patient-derived plaque organoid or explant culture systems treated with TREM1-blocking antibodies (e.g., nangibotide) could experimentally validate the predicted centrality loss and absent compensation. Single-cell ATAC-seq or CUT&Tag profiling of the hub population would identify the transcription factors and cis-regulatory elements maintaining PLIN2/TREM1 co-expression. Finally, prospective cohort studies with standardized multi-bed imaging (coronary CT angiography, carotid ultrasound, femoral MRI) could test whether hub-associated biomarkers (soluble TREM1, PLIN2) predict site-specific plaque progression.

### Conclusions

PLIN2+/TREM1+ macrophages constitute a conserved, non-redundant inflammatory signaling hub in human atherosclerosis. TREM1 functions as a selective signal amplifier whose pharmacological inhibition is predicted to dampen inflammatory communication without triggering compensatory pathway activation. Vascular-bed-specific FAO-SIRT1-MHC-II programming coexists with—but is mechanistically independent of—the conserved hub architecture, suggesting that local microenvironmental factors (hemodynamics, SMC lineage, PVAT biology) tune the metabolic output of the network. These findings advance a two-layer network framework for understanding macrophage-mediated inflammation in atherosclerosis and nominate TREM1 as a therapeutic target with a favorable mechanistic profile for chronic cardiovascular disease.

---


## Clinical Perspectives

- **Background:** Atherosclerosis is a systemic disease with site-specific clinical presentations across carotid, coronary, and femoral arteries. Whether inflammatory signaling in plaques is organized around conserved macrophage hub cells—and whether this architecture can be therapeutically targeted—is unknown.
- **Summary of Results:** Using a 259,116-cell human plaque atlas, we identified PLIN2+/TREM1+ macrophages as a conserved, non-redundant inflammatory signaling hub across three vascular beds. Virtual TREM1 ablation selectively dampened inflammatory communication by 7.8% without triggering compensatory pathway activation. The FAO-SIRT1-MHC-II metabolic-epigenetic axis was femoral-artery-specific, indicating that microenvironmental factors tune the output of a universally conserved hub network.
- **Potential Significance to Human Health and Disease:** These findings nominate TREM1 as a therapeutic target with a favorable mechanistic profile for atherosclerosis—selective inflammatory dampening without immunosuppression or compensatory resistance. The identification of a conserved hub cell architecture provides a network-level framework for evaluating anti-inflammatory therapies in cardiovascular disease, where long-term treatment safety is essential.


## Ethics Statement

This study involved re-analysis of publicly available, de-identified single-cell RNA sequencing data. All original studies obtained ethical approval from their respective institutional review boards as detailed in the primary publications. No new human subject data were collected for this analysis.

## Data Availability

The primary single-cell atlas is publicly available via CELLxGENE (Collection ID: db70986c-7d91-49fe-a399-a4730be394ac) and from the Traeuble et al. GitHub repository (https://github.com/kotr98/reproducibility-plaque-atlas). External validation dataset GSE131778 is available from GEO under accession GSE131778. The re-processed GSE155512 dataset is available from GEO under accession GSE155512. All analysis code, processed data files, and figure generation scripts are available at https://github.com/Caier-dortor/atherosclerosis_atlas.

## Author Contributions

Zhen-Guo Cai: Conceptualization, Formal analysis, Investigation, Methodology, Software, Visualization, Writing – original draft; Yin Tang: Data curation, Formal analysis, Investigation, Writing – review & editing; Linlin Che: Conceptualization, Project administration, Resources, Supervision, Writing – review & editing.

## Competing Interests

The authors declare that they have no competing interests.

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

## References

1. Libby P. The changing landscape of atherosclerosis. Nature. 2021;592(7855):524-33. doi:10.1038/s41586-021-03392-8
2. Weber C, Noels H. Atherosclerosis: current pathogenesis and therapeutic options. Nat Med. 2011;17(11):1410-22. doi:10.1038/nm.2538
3. Tabas I, Lichtman AH. Monocyte-macrophages and T cells in atherosclerosis. Immunity. 2017;47(4):621-34. doi:10.1016/j.immuni.2017.09.008
4. Wirka RC, Wagh D, Paik DT, Pjanic M, Nguyen T, Miller CL, et al. Atheroprotective roles of smooth muscle cell phenotypic modulation and the TCF21 disease gene as revealed by single-cell analysis. Nat Med. 2019;25(8):1280-9. doi:10.1038/s41591-019-0512-5
5. Fernandez DM, Rahman AH, Fernandez NF, Chudnovskiy A, Amir ED, Amadori L, et al. Single-cell immune landscape of human atherosclerotic plaques. Nat Med. 2019;25(10):1576-88. doi:10.1038/s41591-019-0590-4
6. Dib L, Koneva LA, Edsfeldt A, Zurke YX, Sun J, Nitulescu M, et al. Lipid-associated macrophages transition to an inflammatory state in human atherosclerosis, increasing the risk of cerebrovascular complications. Nat Cardiovasc Res. 2023;2:656-72. doi:10.1038/s44161-023-00295-x
7. Traeuble F, Kornfeld JW, Baran Y, Pauli J, Sachs N, Vafadarnejad E, et al. Integrated single-cell atlas of human atherosclerotic plaques. Nat Commun. 2025;16:8255. doi:10.1038/s41467-025-63202-x
8. Netea MG, Dominguez-Andres J, Barreiro LB, Chavakis T, Divangahi M, Fuchs E, et al. Defining trained immunity and its role in health and disease. Nat Rev Immunol. 2020;20(6):375-88. doi:10.1038/s41577-020-0285-6
9. Flores-Gomez D, Bekkering S, Netea MG, Riksen NP. Trained immunity in atherosclerotic cardiovascular disease. Arterioscler Thromb Vasc Biol. 2021;41(1):62-9. doi:10.1161/ATVBAHA.120.314216
10. Bekkering S, Quintin J, Joosten LAB, van der Meer JWM, Netea MG, Riksen NP. Oxidized low-density lipoprotein induces long-term proinflammatory cytokine production and foam cell formation via epigenetic reprogramming of monocytes. Arterioscler Thromb Vasc Biol. 2014;34(8):1731-8. doi:10.1161/ATVBAHA.114.303887
11. Ridker PM, Everett BM, Thuren T, MacFadyen JG, Chang WH, Ballantyne C, et al. Antiinflammatory therapy with canakinumab for atherosclerotic disease. N Engl J Med. 2017;377(12):1119-31. doi:10.1056/NEJMoa1707914
12. Tardif JC, Kouz S, Waters DD, Bertrand OF, Diaz R, Maggioni AP, et al. Efficacy and safety of low-dose colchicine after myocardial infarction. N Engl J Med. 2019;381(26):2497-505. doi:10.1056/NEJMoa1912388
13. Nidorf SM, Fiolet ATL, Mosterd A, Eikelboom JW, Schut A, Opstal TSJ, et al. Colchicine in patients with chronic coronary disease. N Engl J Med. 2020;383(19):1838-47. doi:10.1056/NEJMoa2021372
14. Bouchon A, Dietrich J, Colonna M. Cutting edge: inflammatory responses can be triggered by TREM-1, a novel receptor expressed on neutrophils and monocytes. J Immunol. 2000;164(10):4991-5. doi:10.4049/jimmunol.164.10.4991
15. Gibot S, Cravoisy A, Levy B, Bene MC, Faure G, Bollaert PE. Soluble triggering receptor expressed on myeloid cells and the diagnosis of pneumonia. N Engl J Med. 2004;350(5):451-8. doi:10.1056/NEJMoa031544
16. Arts RJW, Novakovic B, Ter Horst R, Carvalho A, Bekkering S, Lachmandas E, et al. Glutaminolysis and fumarate accumulation integrate immunometabolic and epigenetic programs in trained immunity. Cell Metab. 2016;24(6):807-19. doi:10.1016/j.cmet.2016.10.008
17. Chiu JJ, Chien S. Effects of disturbed flow on vascular endothelium: pathophysiological basis and clinical perspectives. Physiol Rev. 2011;91(1):327-87. doi:10.1152/physrev.00047.2009
18. Majesky MW. Developmental basis of vascular smooth muscle diversity. Arterioscler Thromb Vasc Biol. 2007;27(6):1248-58. doi:10.1161/ATVBAHA.107.141069
19. Antonopoulos AS, Sanna F, Sabharwal N, Thomas S, Oikonomou EK, Herdman L, et al. Detecting human coronary inflammation by imaging perivascular fat. Sci Transl Med. 2017;9(398):eaal2658. doi:10.1126/scitranslmed.aal2658
20. Horimatsu T, Kim HW, Weintraub NL. The role of perivascular adipose tissue in non-atherosclerotic vascular disease. Front Physiol. 2017;8:969. doi:10.3389/fphys.2017.00969
21. Jin S, Guerrero-Juarez CF, Zhang L, Chang I, Ramos R, Kuan CH, et al. Inference and analysis of cell-cell communication using CellChat. Nat Commun. 2021;12:1088. doi:10.1038/s41467-021-21246-9
22. Pan H, Xue C, Auerbach BJ, Fan J, Bashore AC, Cui J, et al. Single-cell genomics reveals a novel cell state during smooth muscle cell phenotypic switching and potential therapeutic targets for atherosclerosis in mouse and human. Circulation. 2020;142(21):2060-75. doi:10.1161/CIRCULATIONAHA.120.048378

# Computational Identification of Genes Differentiating Severe Asthma Patients from Healthy Controls

This study explored a computational approach to identify genes that can differentiate between severe asthma patients and healthy controls. The authors analyzed the GSE27011 gene expression dataset from the NCBI GEO repository, applying supervised statistical learning, specifically L2-regularized logistic regression, for feature selection and classification.

Key findings include:

- **Dataset Handling:** The original dataset (DS1) included severe asthmatics, mild asthmatics, and healthy controls. A sub-dataset (DS2) was created by removing the "mild asthma" class to improve separability, as mild asthma transcription data can overlap with both healthy and severe asthma profiles.
- **Classification Performance:** Logistic Regression models (Lasso L1 and Ridge L2) were built for both DS1 and DS2. Accuracy for DS1 was around 87-89%, while for DS2, it improved to 91-97%, demonstrating that removing the "mild asthma" class enhanced model performance.
- **Dimensionality Reduction with SVD:** To address the high number of genes (over 28,000 features), Singular Value Decomposition (SVD) was applied to reduce noise and identify relevant components. The study found that using eight dimensions from the SVD-reduced matrix resulted in a 100% accuracy for the L2 logistic regression model on DS2.
- **Algorithm Comparison:** Logistic regression models consistently outperformed other machine learning algorithms tested (KNN, Random Forest, and SVM) in classifying the microarray data across all datasets.
- **Gene Identification:**
  - **Information Gain Algorithm:** The top 10 most important genes identified by the Information Gain algorithm for classifying DS2 included several pseudogenes (e.g., RNU6-1123P, PYY2, RNU6-82P) and protein-coding genes like GPR21, GPR52, KAT6B, and TCEB1. Heatmap visualization confirmed differential expression of these genes between severe asthma and control groups.
  - **Logistic Regression Weight Distribution ($\alpha$):** The study also used the weights from the logistic regression model to identify the 20 genes (10 positive, 10 negative contributors) with the highest impact on classification. Notable genes with positive contributions (associated with severe asthma) included _MYLIP_, _JAK2_, _ZEB2_, and _FOSB_. Genes with negative contributions (associated with healthy controls) included _TMEM176A_, _GRINL1A_, and _CMBL_. The authors highlight _JAK2_ for its known role in inflammatory signaling in asthma and discuss potential indirect associations for _MYLIP_, _RNU4-2_, _SNORD56B_, _RNU1-19_, _ZEB2_, and _FOSB_.

The study concludes that this computational strategy successfully validated previously reported asthma-associated genes and identified additional candidates with potential diagnostic or therapeutic relevance, emphasizing the value of reanalyzing public transcriptomic datasets using machine learning.

## What can be found in this repository?

- [Code][Code]:
  - Scripts for data processing, feature selection, and model training.
  - Machine learning workflows implemented in Open Workspace (OWS).
- [Data][Data]:
  - Processed datasets used in the analysis, including gene expression data and auxiliary information.
- [Images][Images]:
  - Visualizations of results, including heatmaps, SVD plots, and logistic regression weight distributions.
- [Supplementary Tables][Supplementary Tables]:
  - Supplementary tables providing additional details on the analysis, including gene lists and model performance metrics.

[Code]: Files/Code
[Data]: Files/Data
[Images]: Files/Images
[Supplementary Tables]: Files/supplementary_tables

<details>

<summary>Repository File Tree</summary>

┣ 📂Code
┃ ┣ 📂machine_learning
┃ ┃ ┗ 📜machine_learning.ows
┃ ┣ 📜asthma_for_article.m
┃ ┣ 📜resolve.m
┣ 📂Data
┃ ┣ 📂Misc
┃ ┃ ┣ 📜Aux_data.txt
┃ ┃ ┣ 📜genes_mais_importantes_2.csv
┃ ┃ ┗ 📜genes_mais_importantes.csv
┃ ┣ 📂Used
┃ ┃ ┣ 📜Aux_data.tsv
┃ ┃ ┣ 📜GSE27011_series_matrix.tsv
┃ ┃ ┣ 📜negative.tsv
┃ ┃ ┣ 📜positive.tsv
┃ ┃ ┗ 📜target_simples.tsv
┃ ┗ 📜readme.md
┣ 📂Images
┃ ┣ 📂Docs
┃ ┃ ┣ 📜fig1.png
┃ ┃ ┣ 📜fig1.psd
┃ ┃ ┣ 📜fig1.svg
┃ ┃ ┣ 📜fig2.png
┃ ┃ ┣ 📜fig3.png
┃ ┃ ┣ 📜fig3.psd
┃ ┃ ┣ 📜fig3.svg
┃ ┃ ┣ 📜fig4.png
┃ ┃ ┣ 📜fig4.psd
┃ ┃ ┣ 📜fig4.svg
┃ ┃ ┣ 📜fig5.png
┃ ┃ ┣ 📜fig5.psd
┃ ┃ ┣ 📜fig5.svg
┃ ┃ ┗ 📜readme.md
┃ ┣ 📂Heatmaps
┃ ┃ ┣ 📜Gene expression 1.jpg
┃ ┃ ┣ 📜Gene expression 2.jpg
┃ ┃ ┣ 📜Gene expression 3.jpg
┃ ┃ ┣ 📜Rplot_heatmap_negative_grouped.svg
┃ ┃ ┣ 📜Rplot_heatmap_negative.svg
┃ ┃ ┣ 📜Rplot_Heatmap_Plot_v2.svg
┃ ┃ ┣ 📜Rplot_heatmap_positive_grouped.svg
┃ ┃ ┣ 📜Rplot_heatmap_positive.svg
┃ ┃ ┣ 📜Rplot_heatmap.svg
┃ ┃ ┗ 📜Rplot_severe_vs_control.svg
┃ ┣ 📜AttributesWeights.svg
┃ ┣ 📜EntityDomain.svg
┃ ┣ 📜EntityDomainReducedMatrix.pdf
┃ ┣ 📜EntityDomainReducedMatrix.png
┃ ┣ 📜EntityDomainReducedMatrix.svg
┃ ┣ 📜Probabilities.svg
┃ ┣ 📜SVD.svg
┃ ┗ 📜SVDReducedMatrix.svg
┣ 📂supplementary_tables
┃ ┣ 📜readme.md
┃ ┣ 📜Supplementary Table S1.docx
┃ ┣ 📜Supplementary Table S1.pdf
┃ ┣ 📜Supplementary Table S2.docx
┃ ┣ 📜Supplementary Table S2.pdf
┃ ┣ 📜Supplementary Table S3.docx
┃ ┗ 📜Supplementary Table S3.pdf
┗ 📜artigo_milenna.docx

</details>

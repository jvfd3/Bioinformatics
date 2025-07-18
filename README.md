# Análise de Dados Transcriptômicos de Asma com Regressão Logística e SVD

Os autores analisaram o conjunto de dados GSE27011 do repositório NCBI GEO, aplicando regressão logística regularizada (L1 e L2) para seleção de características e classificação. A metodologia incluiu a redução de dimensionalidade utilizando a técnica de Decomposição de Valor Singular (SVD) no MATLAB.

Os resultados mostraram que a remoção da classe "asma leve" (criando o subconjunto DS2) melhorou a precisão dos modelos de regressão logística, alcançando 97% de acurácia com L2 para o DS2. A aplicação de SVD, utilizando 8 dimensões, resultou em 100% de acurácia para o modelo L2 no DS2, demonstrando a eficácia da redução de ruído.

Comparando com outros algoritmos de aprendizado de máquina (KNN, SVM e Random Forest), a regressão logística (especialmente Ridge L2 com SVD) apresentou o melhor desempenho. O estudo também identificou genes importantes para a classificação usando o algoritmo de Ganho de Informação, além de uma análise de pesos da regressão logística para destacar genes com maior impacto positivo e negativo na classificação. Entre os genes relevantes, JAK2 foi apontado por seu papel central na sinalização inflamatória.

Em conclusão, o estudo validou genes previamente associados à asma e identificou candidatos adicionais com potencial diagnóstico ou terapêutico, enfatizando a utilidade de mineração de dados transcriptômicos públicos com aprendizado de máquina para descobrir assinaturas moleculares de doenças complexas.

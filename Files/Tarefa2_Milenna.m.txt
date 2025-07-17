% ================================================================
% RESUMO DO SCRIPT - ETAPA DE EXPLORAÇÃO DIMENSIONAL (SVD/PCA)
%
% Objetivo: Verificar se os perfis de expressão gênica de pacientes 
% controle e pacientes com asma severa são separáveis, antes de aplicar
% um modelo de classificação.
%
% Etapas realizadas:
% 1. Seleção das amostras:
%    - Usa apenas controles (colunas 1 a 19) e pacientes com asma severa 
%      (colunas 39 ao final), excluindo os com asma moderada.
%
% 2. Preparação da matriz:
%    - Transpõe a matriz para que cada linha seja uma amostra e cada coluna, um gene.
%    - Centraliza os dados subtraindo a média de cada gene (coluna).
%
% 3. Decomposição em valores singulares (SVD):
%    - Aplica SVD para obter os componentes principais (PCs).
%    - Calcula a proporção de variância explicada por cada componente.
%
% 4. Análise da separação dos grupos:
%    - Calcula os centróides dos grupos controle e asma severa no espaço dos 3 primeiros PCs.
%    - Mede a distância entre esses centróides.
%    - Realiza um teste ANOVA para verificar se há diferença significativa entre os grupos no PC1.
%
% 5. Visualização:
%    - Cria um gráfico de dispersão das amostras no espaço PC1 x PC2, destacando os grupos.
%
% Essa análise serve como etapa exploratória para avaliar se os dados
% têm potencial para serem bem classificados por modelos supervisionados.
% ================================================================


% === Primeira Etapa - Decomposição ===
% Escolher um microarray dentre os disponíveis no NCBI. Fazer o download da matriz, caracterizando o problema subjacente ao microarray como um vector space model, 
% onde se deseja resolver problemas de classificação. Descrever as entidades e o tipo de atributos utilizados. Resolver eventuais omissões como os atributos cujos valores não foram informados.

%Seleção das colunas da matriz: aqui estou selecionando apenas os pacientes
%controle (sem asma, do 1 ao 19) e os pacientes com asma severa (39 ao
%final), excluo os pacientes com asma moderada 
Ans = [MatrizNCBI(:,1:19) MatrizNCBI(:,39:end)];

%SVD (decomposição em valores singulares) - usar os componentes para visualizar variância e calcular centróides, 
% é uma etapa de exploração dimensional.

Ans = Ans'

Ans = Ans - mean(Ans, 1);  % Subtrai a média de cada coluna (gene)

[T, S, V] = svd(Ans, "econ");

diagonal_S = diag(S); % Pega a diagonal principal
dist_import_relativa = diagonal_S/sum(diagonal_S);

figure
hold on
    %title('Valores Singulares Relativos da Matriz "Ans" dos Severos') % Define o título
    grid % Liga o grid
    plot(dist_import_relativa, '*') % Marca os pontos
    plot(dist_import_relativa) % Faz a linha
hold off

% A matriz PC_scores representa as entidades (amostras) projetadas no novo espaço dos PCs (35x35). 
PC_scores = T * S;   % PC_scores será 35x35

%Esses "scores" são as representações das suas 35 amostras em um novo espaço de menor dimensão, 
% onde os dados estão reorganizados para capturar a maior variação possível nos primeiros componentes (tipo PCA).

% Labels: 0 para controle, 1 para asma severa
labels = [zeros(19,1); ones(17,1)];

% Índices
controle_idx = (labels == 0);
asma_idx = (labels == 1);

% Cálculo dos centróides nos 3 primeiros PCs - Distância entre centróides
% para verificar se há separação entre os grupos no espaço projetado. 
% Uma distância alta sugere que os grupos são separáveis com base nas expressões gênicas.
centroide_controle = mean(PC_scores(controle_idx, 1:3), 1);
centroide_asma = mean(PC_scores(asma_idx, 1:3), 1);

%Um centróide é o "ponto médio" (média das posições) das amostras de um grupo. 
% obtém um resumo vetorial da posição média dos controles e dos asmáticos no novo espaço reduzido.

% Distância entre centróides
distancia_centroide = norm(centroide_controle - centroide_asma);

fprintf('Distância entre centróides: %.4f\n', distancia_centroide);

group = categorical([repmat("controle",19,1); repmat("asma",17,1)]);
[p, tbl, stats] = anova1(PC_scores(:,1), group, 'off');

fprintf('p-valor da ANOVA no PC1: %.4e\n', p);




% === Etapa 2 ===
% Construção de um Classificador: usar um modelo de regressão logística modificada para proceder à classificação dos problemas associados, considerando a matriz de entidades completa.


Aux = S*V';
x = Aux(1,:);
y = Aux(2,:);
z = Aux(3,:);

matriz_3x36 = [x; y; z];

% Criação do gráfico
figure
hold on
    title('Projeção das Amostras no Espaço PC1 x PC2')
    xlabel('Primeiro Componente Principal (PC1)')
    ylabel('Segundo Componente Principal (PC2)')
    grid on

    % Plotar grupo controle (label == 0) com círculo preenchido
    scatter(x(labels==0), y(labels==0), 50, 'r', 'filled', 'DisplayName','Controle')

    % Plotar grupo asma (label == 1) com círculo não preenchido
    scatter(x(labels==1), y(labels==1), 50, 'r', 'DisplayName','Asma severa')

    legend('show')
hold off

% ================================================================
% RESUMO DO SCRIPT - REGRESSÃO LOGÍSTICA PARA CLASSIFICAÇÃO BINÁRIA
%
% Objetivo: Estimar os pesos (coeficientes) dos atributos (genes)
% que melhor classificam as amostras entre controle e asma severa,
% usando um modelo de regressão logística linear.
%
% Descrição das etapas:
%
% 1. Função 'resolve':
%    - Implementa a solução do sistema linear para encontrar os coeficientes
%      alpha da regressão logística modificada, utilizando um sistema de equações
%      baseado na matriz de dados A (genes x amostras) e no vetor b (rótulos codificados).
%
% 2. Codificação dos rótulos (b):
%    - Define valores log-odds (lgch1 = 12 para controle e lgch0 = -12 para asma severa),
%      que são transformações logarítmicas para a regressão logística.
%
% 3. Chamada da função 'resolve':
%    - Calcula o vetor de coeficientes alpha que representa o peso de cada gene na classificação.
%
% 4. Visualização dos pesos alpha:
%    - Plota os pesos para observar quais genes têm maior influência no modelo.
%
% 5. Cálculo das probabilidades de classificação:
%    - Calcula o valor linear aux = Ans' * alpha.
%    - Aplica a função logística para transformar em probabilidades p entre 0 e 1.
%
% 6. Visualização das probabilidades:
%    - Plota as probabilidades para cada amostra, indicando a confiança do modelo
%      na classificação entre controle e asma severa.
%
% Essa etapa é fundamental para identificar genes importantes e avaliar a capacidade
% do modelo de regressão logística em separar os dois grupos.
% ================================================================

% Seleção dos 10 marcadores mais importantes

function [alpha, x] = resolve(A,b)
    [m,n]=size(A);
    Im = sparse(eye(m));
    In = sparse(eye(n));
    M = sparse(m+n,m+n);
    M = [Im,-A;-A',In];
    nb = zeros(m+n,1);
    nb(1:m) = -b;
    x = M\nb;
    alpha = x(m+1:end);
end

%lgch1 = log(0.999/(1-0.999));
%lgch0 = log(0.001/(1-0.001));

lgch1 = 12;
lgch0 = -12;
[m,n] =  size(Ans);
b = zeros (n,1);
b(1:19) = lgch1;
b(20:end) = lgch0;

alpha = resolve2(Ans', b);

%Esse gráfico mostra o vetor alpha, que contém os pesos 
% (ou coeficientes) calculados pelo seu modelo de regressão logística para cada gene (atributo).

%Eixo X: posição do gene (ex: gene 1, gene 2, ..., gene 36).

%Eixo Y: valor do peso (alpha) atribuído a esse gene no modelo.


figure
hold on
    title('Pesos associados aos atributos')
    % grid
    plot(alpha, '*')
hold off

[~, idx_mais_importantes] = maxk(abs(alpha), 34);

%Verificando em que linha meus dados começam

filename = 'GSE27011_series_matrix.txt';

fid = fopen(filename);
linha_num = 0;
start_line = 0;

while true
    linha = fgetl(fid);
    if ~ischar(linha)
        break;
    end
    linha_num = linha_num + 1;
    
    if contains(linha, 'ID_REF') && contains(linha, 'GSM')  % linha com cabeçalho de dados
        start_line = linha_num;
        break;
    end
end
fclose(fid);

fprintf('A matriz de dados começa na linha %d\n', start_line);

%Identificando os genes
gene_ids = T{:,1};  % Primeira coluna da tabela (IDs das sondas)
idx_mais_importantes = [10 14 17 31 5 33 28 20 12 27 11 22 18 4 24 36 15 26 34 1 8 23 6 32 16 2 29 30 9 13 35 19 25 3];
top10_sondas = gene_ids(idx_mais_importantes)

is_top10 = ismember(Tmap.ID, top10_sondas);  % tira o string()
genes_mais_importantes = Tmap(is_top10, :)

writetable(genes_mais_importantes, 'C:\Users\milen\OneDrive\Documentos\Doutorado\Disciplinas\Algoritmos para Bioinformática I\genes_mais_importantes.csv');


aux = Ans'* alpha;
num = exp(aux);
p = num./(1+num);


figure
hold on
    title('Classificação por Regressão Logística - Todos os Atributos')
    plot(p, '*')
    xlim([min(p)-0.01, max(p)+0.01])
    ylim([min(p)-0.01, max(p)+0.01])  % margem para ficar proporcional
    xlabel('Genes')
    ylabel('Valores de \alpha_i', 'Interpreter', 'tex')
    axis tight  % força ajuste mais justo
    axis square % força proporção 1:1 (opcional, só se fizer sentido)

hold off

%Outra forma (IA):
figure
hold on
    title('Classificação por Regressão Logística - Todos os Atributos')
    xlabel('Índice da amostra')
    ylabel('Probabilidade de controle')
    grid on

    % Plotar controles em azul
    scatter(find(labels==0), p(labels==0), 50, 'b', 'filled', 'DisplayName','Controle')


%validar os marcadores
[valores, pos] = sort(alpha); %pos é onde estão os valores negativos e positivos 
escolha = [pos(1:7)' pos(end-6:end)'];
Matrizreduzida = [Ans(escolha, :)];

% whos Matrizreduzida

[T, S, V] = svd(Matrizreduzida, "econ");

diagonal_S = diag(S); % Pega a diagonal principal
dist_import_relativa = diagonal_S/sum(diagonal_S);

figure
hold on
    title('Valores Singulares Relativos da Matriz Reduzida') % Define o título
    grid % Liga o grid
    plot(dist_import_relativa, '*') % Marca os pontos
    plot(dist_import_relativa) % Faz a linha
hold off

Aux = S*V';
x = Aux(1,:);
y = Aux(2,:);
z = Aux(3,:);

figure
hold on
    title('Visualização do domínio das entidades')
    grid
    plot3(x, y, z, 'or')
    %hold on
    plot3(x(1:20), y(1:20), z(1:20),'*r')
hold off


%P(x) da matriz reduzida

alphanovo = Matrizreduzida'\b;
title('Pesos associados aos atributos escolhidos')
plot(alphanovo, '*')

%Ver se o gráfico P(x) está longe do zero

aux = Matrizreduzida' * alphanovo;
num = exp(aux);
p = num ./ (1+num);

figure
hold on
    title('Classificação por Regressão Logística - Atributos do modelo reduzido')
    plot(p, '*')
hold off

d = pdist(Aux(1:3,:)');
L = linkage(d);
dendrogram(L)



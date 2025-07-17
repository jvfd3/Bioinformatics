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

%Visualização antes da regressão logística modificada

Ans_Teste = [MatrizNCBI(:,1:19), MatrizNCBI(:,39:end)]';
Ans_Teste = Ans_Teste - mean(Ans, 1);  % Centraliza

[U, S, V] = svd(Ans_Teste, 'econ');

PC_scores_Teste = U * S;  % Essa é a projeção com escala real

xt = PC_scores_Teste(:,1);
yt = PC_scores_Teste(:,2);

labels_Teste = [zeros(19,1); ones(17,1)];

figure
hold on
    title('Projeção das Amostras no Espaço PC1 x PC2 (via SVD)')
    xlabel('PC1')
    ylabel('PC2')
    grid on

    scatter(xt(labels==0), yt(labels==0), 60, 'r', 'filled', 'DisplayName','Controle');
    scatter(xt(labels==1), yt(labels==1), 60, 'k', 'o', 'LineWidth', 1.2, 'DisplayName','Asma severa');

    legend('Location','best')
hold off






% ======
% Selection of the 10 most important markers

function [alpha, x] = solve_system(A, b)
    [m, n] = size(A);
    Im = sparse(eye(m));
    In = sparse(eye(n));
    M = sparse(m+n, m+n);
    M = [Im, -A; -A', In];
    nb = zeros(m+n,1);
    nb(1:m) = -b;
    x = M \ nb;
    alpha = x(m+1:end);
end

% Modified logistic regression thresholds
logit_ch1 = 12;
logit_ch0 = -12;

[m, n] = size(DataMatrix);
b = zeros(n,1);
b(1:19) = logit_ch1;
b(20:end) = logit_ch0;

alpha = solve_system(DataMatrix', b);

figure
hold on
    title('Weights Associated with Attributes')
    plot(alpha, '*')
hold off

aux = DataMatrix' * alpha;
num = exp(aux);
p = num ./ (1 + num);

figure
hold on
    title('Logistic Regression Classification - All Attributes')
    plot(p, '*')
hold off

% Validate the markers
[values, positions] = sort(alpha); % 'positions' points to negative and positive values
selected = [positions(1:7)' positions(end-6:end)'];
ReducedMatrix = [DataMatrix(selected, :)];

% whos ReducedMatrix

[T, S, V] = svd(ReducedMatrix, "econ");

singular_values = diag(S); % Extract the main diagonal
relative_importance = singular_values / sum(singular_values);

figure
hold on
    title('Relative Singular Values of the Reduced Matrix')
    grid on
    plot(relative_importance, '*')
    plot(relative_importance)
hold off

AuxiliaryMatrix = S*V';
x = AuxiliaryMatrix(1,:);
y = AuxiliaryMatrix(2,:);
z = AuxiliaryMatrix(3,:);

figure
hold on
    title('Entity Domain Visualization')
    grid on
    plot3(x, y, z, 'or')
    plot3(x(1:20), y(1:20), z(1:20), '*r')
hold off

% P(x) for the reduced matrix

new_alpha = ReducedMatrix' \ b;
title('Weights Associated with Selected Attributes')
plot(new_alpha, '*')

% Check if P(x) plot is far from zero

aux = ReducedMatrix' * new_alpha;
num = exp(aux);
p = num ./ (1 + num);

figure
hold on
    title('Logistic Regression Classification - Selected Attributes')
    plot(p, '*')
hold off

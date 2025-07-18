% ================================================================
% RESUMO DO SCRIPT - ETAPA DE EXPLORAÇÃO DIMENSIONAL (SVD)
%
% Objetivo: Verificar se os perfis de expressão gênica de pacientes 
% controle e pacientes com asma severa são separáveis, antes de aplicar
% um modelo de classificação.
%
% Etapas realizadas:
% 1. Seleção das amostras: Usa apenas controles (colunas 1 a 19) e pacientes com asma severa (colunas 39 ao final), excluindo os com asma moderada.
% 2. Decomposição em valores singulares (SVD);
% 3. Visualização: Cria um gráfico de dispersão das amostras no espaço, destacando os grupos.
%
% Essa análise serve como etapa exploratória para avaliar se os dados
% têm potencial para serem bem classificados por modelos supervisionados.
% ================================================================


% === Primeira Etapa - Decomposição ===
% Seleção das colunas da matriz: aqui estou selecionando apenas os pacientes
% controle (sem asma, do 1 ao 18) e os pacientes com asma severa (38 ao
% final)
Ans = [MatrizNCBI(:,1:18) MatrizNCBI(:,38:end)];

% SVD (decomposição em valores singulares) 

[T, S, V] = svd(Ans, "econ"); 

diagonal_S = diag(S); % Pega a diagonal principal
dist_import_relativa = diagonal_S/sum(diagonal_S);

figure
hold on
    title('Valores Singulares Relativos da Matriz "Ans" dos Severos') % Define o título
    grid % Liga o grid
    plot(dist_import_relativa, '*') % Marca os pontos
    plot(dist_import_relativa) % Faz a linha
hold off

% === Etapa 2 ===
% Visualização da dispersão dos grupos
Aux = S*V';
x = Aux(1,:);
y = Aux(2,:);
z = Aux(3,:);

% Criação do gráfico
figure
hold on
    title('Domínio das Entidades')
    xlabel('Produto do Primeiro Padrão')
    ylabel('Produto do Segundo Padrão')
    plot3(x,y,z,'or')
    hold on
    plot3(x(1:18), y(1:18), z(1:18),'+r') % Alterar 18 para o nº do seu primeiro grupo
    grid on
hold off

% ================================================================
% RESUMO DO SCRIPT - REGRESSÃO LOGÍSTICA MODIFICADA PARA CLASSIFICAÇÃO
%
% Objetivo: Estimar os pesos (coeficientes) dos atributos (genes)
% que melhor classificam as amostras entre controle e asma severa,
% usando um modelo de regressão logística linear.
%
% Descrição das etapas:
%
% 1. Função 'resolve':
%    - Implementa a solução do sistema linear para encontrar os coeficientes alpha da regressão logística modificada, utilizando um sistema de equações baseado na matriz de dados A (genes x amostras) e no vetor b (rótulos codificados).
%
% 2. Codificação dos rótulos (b):
%    - Define valores log-odds (lgch1 para controle e lgch0 para asma severa),
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
%    - Plota as probabilidades para cada amostra, indicando a confiança do modelo na classificação entre controle e asma severa.
%
% Essa etapa é fundamental para identificar genes importantes e avaliar a capacidade
% do modelo de regressão logística em separar os dois grupos.
% ================================================================

% lgch1 = log(0.999/(1-0.999));
% lgch0 = log(0.001/(1-0.001));

lgch1 = 12;
lgch0 = -12;
[m,n] =  size(Ans);
b = zeros (n,1);
b(1:18) = lgch1;
b(19:end) = lgch0;

alpha = resolve(Ans', b);

% Criação do Gráfico

figure
hold on
    title('Pesos associados aos atributos')
    plot(alpha, '*')
hold off

aux = Ans'* alpha;
num = exp(aux);
p = num./(1+num);

% Criação do gráfico 
figure
plot(p,'*')
hold off

% Seleção dos 10 marcadores mais importantes
[valores, pos] = sort(alpha); % pos é onde estão os valores negativos e positivos 
escolha = [pos(1:10)' pos(end-10:end)'];
Matrizreduzida = [Ans(escolha, :)];

% Decomposição da Matriz Reduzida
[T1, S2, V3] = svd(Matrizreduzida, "econ");

diagonal_S = diag(S); % Pega a diagonal principal
dist_import_relativa = diagonal_S/sum(diagonal_S);

% Criação do gráfico de decomposição da matriz reduzida
figure
hold on
    title('Valores Singulares Relativos da Matriz Reduzida') % Define o título
    grid % Liga o grid
    plot(dist_import_relativa, '*') % Marca os pontos
    plot(dist_import_relativa) % Faz a linha
hold off

Auxreduzido = S2*V3';
x = Auxreduzido(1,:);
y = Auxreduzido(2,:);
z = Auxreduzido(3,:);

figure
hold on
    title('Visualização do domínio das entidades')
    grid
    plot3(x, y, z, 'or')
    % hold on
    plot3(x(1:18), y(1:18), z(1:18),'*r')
hold off


% P(x) da matriz reduzida

alphanovo = Matrizreduzida'\b;
figure
    title('Pesos associados aos atributos escolhidos')
    plot(alphanovo, '*')
grid on 
% Ver se o gráfico P(x) está longe do zero

aux3 = Matrizreduzida' * alphanovo;
num = exp(aux3);
p = num ./ (1+num);

figure
hold on
    title('Classificação por Regressão Logística - Atributos do modelo reduzido')
    plot(p, '*')
hold off

d = pdist(Aux(1:3,:)');
L = linkage(d);
dendrogram(L)

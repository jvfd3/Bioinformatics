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

[~, idx_mais_importantes] = maxk(abs(alpha), 10);

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
idx_mais_importantes = [10 14 17 31 5 33 28 20 12 27];
top10_sondas = gene_ids(idx_mais_importantes)

is_top10 = ismember(Tmap.ID, top10_sondas);  % tira o string()
genes_mais_importantes = Tmap(is_top10, :)
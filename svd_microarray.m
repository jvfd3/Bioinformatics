% Choose a microarray from those available in NCBI. Download the matrix, characterizing 
% the underlying problem as a vector space model aimed at solving classification problems.
% Describe the entities and the type of attributes used. Handle missing data if any.

DataMatrix = [MatrizNCBI(:,1:19) MatrizNCBI(:,39:end)];

[T, S, V] = svd(DataMatrix, "econ");

singular_values = diag(S); % Extract the main diagonal
relative_importance = singular_values / sum(singular_values);

figure
hold on
    title('Relative Singular Values of Matrix A')
    grid on
    plot(relative_importance, '*')
    plot(relative_importance)
hold off

PC_scores = T * S;   % PC_scores will be 35x35

% Labels: 0 for control, 1 for severe asthma
labels = [zeros(19,1); ones(16,1)];

% Indices
control_idx = (labels == 0);
asthma_idx = (labels == 1);

% Calculate centroids in the first 3 PCs
centroid_control = mean(PC_scores(control_idx, 1:3), 1);
centroid_asthma = mean(PC_scores(asthma_idx, 1:3), 1);

% Distance between centroids
centroid_distance = norm(centroid_control - centroid_asthma);

fprintf('Distance between centroids: %.4f\n', centroid_distance);

% ======
% Classifier Construction: Use a modified logistic regression model to classify 
% the associated problems considering the complete entity matrix.

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

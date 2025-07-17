df = cleandata
[T, S, V] = svd(df, "econ")
s = diag(S)
s = s*(1/sum(s));
figure
plot(s, '*')
hold on
plot(s)
title('Valores Singulares Relativos da Matriz A')
grid on
hold off


Aux = S*V;
x = Aux(1,:);
y = Aux(2,:);
z = Aux(3,:);
figure
title('Visualização do dompinio das entidades')
grid on
plot3(x,y,z,'or')
hold on
plot3(x(1:47), y(1:47), z(1:47),'+r') %Alterar 47 para o n do seu primeiro grupo
hold off
lgch1 = log(0.999/(1-.0999))
lgch0 = log(0.001/(1-0.001))
[m,n] =  size(df)
b = zeros (n,1);
b(1:47) = lgch1 %Alterar 47 para o n do seu primeiro grupo
b(47:end) = lgch0 %Alterar 47 para o n do seu primeiro grupo

alpha = resolve(df',b) %Erro

figure
plot(alpha, '*')
title('Pesos associados aos atributos')
hold off
aux = A'* alpha;
num = exp(aux)
p = num./1+num;
figure
plot(p,'*')
hold off
[valores,pos] = sort(alpha);
Ar = [df([pos(1:7),pos(end-6:end)],:)];
alphar = Ar'\b;
aux = Ar'*alphar;
num = exp(aux);
pr = num./(1+num);
figure
title('P após seleção de atributos')

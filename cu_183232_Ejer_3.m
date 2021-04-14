% Ejercicio 3
% Santiago Fernández Gutiérrez Zamora
% 183232

% Tchebysheff

n = 10;
a = -1;
b = 1;
eps = 0.01;
hf = @(x) sin(pi.*x)./(pi.*x);

k = 0:n-1;
k = fliplr(k);
w = pi/(2*n);
theta_k = (2*k + 1)*w;
uk = cos(theta_k) % los puntos (3.1)
        
MT = ones(n);
MT(:,2) = uk;
for j = 3:n
    for r = 1:n
    	MT(r,j) = 2*(uk(r)*MT(r,j-1)) - MT(r,j-2); 
    end
end

MT % los vectores (3.1)

% aproximación de sinc(x) = sin(pi*x) / pi*x
        
xk = ((b-a)/2)*uk + (a+b)/2;
y = hf(pi*xk); % es por pi*x
coefs = y*0; 
mult = 1;
for l = 1:n
    coefs(l) = (sum(MT(:,l).*y')/n)*mult;
    if l == 1
        mult = 2;
    end
end

x = a:eps:b; % el intervalo
l = length(x);

u = (2*x-a-b)/(b-a);

MT1 = ones(l,n);
MT1(:,2) = u;
for j = 3:n
    for r = 1:l
        MT1(r,j) = 2*(u(r)*MT1(r,j-1)) - MT1(r,j-2); 
    end
end

y = x.*0;
for t = 1:l
   y(t) = sum(coefs(1:n).*MT1(t,:)); 
end

yr = hf(pi*x);

plot(x,y);
hold on
plot(x,yr);
legend({'Aproximación','Función'},'Location','southwest')
hold off
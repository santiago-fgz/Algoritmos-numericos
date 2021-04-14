function A = aprox_sistEc(datos)
% Entrega los coeficientes del polinomio que pasa por los puntos dados
%
% input:
% -datos: matriz de nx2 con los puntos (xi,yi)
%
% output:
% -A: vector con n elementos, los coeficientes del polinomio:
%    p(x) = A0 + A1x + A2(x^2) + ... + A(n-1)(x^n)
%
% Santiago F.G. Zamora

    mn = size(datos);
    n = mn(1);
    
    if n>1
        x = zeros(n);
        y = zeros(n,1);
        for i = 1:n
            xi = datos(i,1);
            y(i) = datos(i,2);
            for j = 1:n
                x(i,j) = xi^(j-1);
            end
        end
        A = x\y;
    else
        A = 0;
    end
end
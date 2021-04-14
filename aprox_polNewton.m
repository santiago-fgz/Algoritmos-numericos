function A = aprox_polNewton(datos)
% Entrega los coeficientes (Ai) para formar el polinomio de newton que
% aproxima a los datos dados.
%
% inputs:
% -datos: matriz de nx2 con los puntos (xi,yi)
%
% outputs:
% -A: Coeficientes (Ai) del polinomio
%
% Ejemplo:
% datos = [x1 y1;x2 y2;x3 y3;x4 y4;...;xn yn]
% El grado del polinomio es siempre tres o menos
% p(x) = A1 + A2*(x-x1) + A3*(x-x1)*(x-x2) + A4*(x-x1)*(x-x2)*(x-x3)
%
% Santiago F.G. Zamora

    mn = size(datos);
    n = mn(1);
    
    if n>1
       f = zeros(n);
       for i = 1:n
          f(1,i) = datos(i,2);
          if  i-1 > 0
              f(2,i) = pend(datos(i-1,:),datos(i,:));
              if i-2 > 0
                 f(3,i) =  pend([datos(i-2,1),f(2,i-1)],[datos(i,1),f(2,i)]);
                 if i-3 > 0
                    f(4,i) = pend([datos(i-3,1),f(3,i-1)],[datos(i,1),f(3,i)]);
                 end
              end
          end
       end
       A = diag(f);
    end
    
end

function m = pend(a,b)
% pendiente
    if a(1) ~= b(1)
       m = (b(2)-a(2))/(b(1)-a(1));
    end
end

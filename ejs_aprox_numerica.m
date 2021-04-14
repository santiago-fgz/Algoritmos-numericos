%==========================================================================
%                 EJERCICIOS DE APROXIMACIÓN NUMÉRICA
%
% Santiago F.G. Zamora                                              183232
%==========================================================================

%% 1) Encontrar el polinomio de menor grado que pase por:

d1 = [0 1;1 3;5 7]; % [x1 y1;x2 y2;x3 y3]

A = aprox_sistEc(d1);
A0 = A(1,1);
A1 = A(2,1);
A2 = A(3,1);
plnm1 = @(x) A0 + A1.*x + A2.*(x.^2);

Ap = aprox_polNewton(d1);
plnm2 = @(x) Ap(1) + Ap(2).*(x-d1(1,1)) + Ap(3).*(x-d1(1,1)).*(x-d1(1,2));


x = 0:0.1:10;
plot(x,plnm1(x))
plot(x,plnm2(x))

hold on
plot(d1(1,1),d1(1,2),'r*')
plot(d1(2,1),d1(2,2),'r*')
plot(d1(3,1),d1(3,2),'r*')
hold off

%% ========================================================================
% 2) Obtener el polinomio y = p(x) de menor grado que cumpla con:
%                 [x y  y']
d2 = [-1 -1 0;2 1 0;5 2 -1];
% ver span
% splines cúbicos
%% ========================================================================
% 3) Valores A y B para la que y(x) = A*exp(B*x) pase por:
 d3 = [-1 1;1 3;4 16];
 
 AB = aprox_logaritmica(d3);
 A = AB(1);
 B = AB(2);
 
 x = -2:0.1:5;
 y = A.*exp(x.*B);
 plot(x,y)
 hold on
 [n,m] = size(d3);
 for j = 1:n
     plot(d3(j,1),d3(j,2),'r*')
 end
 hold off

 %% =======================================================================
 % 4) Valores del polinomio p(x) tq y = exp(p(x)) pase por:
 d4 = [-1 2;0 8;1 1;2 2];
 as = aprox_logaritmica_polinomios(d4);
 
 x = -1:0.1:3;
 y = exp(as(1) + as(2).*x + as(3).*(x.^2) + as(4).*(x.^3));
 plot(x,y)
 hold on
 n = size(d4);
 for j = 1:n(1)
     plot(d4(j,1),d4(j,2),'r*')
 end
 hold off
 
 %% =======================================================================
 % 5) Obtener funcionalente una curva en 2d que pase en orden por:
 d5 = [1 4;2 2;4 3;3 5];
 
 % trayectorias x(t), y(t)
 %% =======================================================================
 % 6) Obtener la ecuación del plano que pasa por los puntos:
 d6 = [1 1;3 -1;2 2];
 
 x = 0:0.1:4;
 hold on
 
 % Por sistemas de ecuaciones
 c01 = aprox_sistEc(d6);
 y01 = c01(1) + c01(2).*x + c01(3).*(x.^2);
 plot(x,y01)
 
 % Polinomio de Newton
 c02 = aprox_polNewton(d6);
 y02 = c02(1) + c02(2).*(x-d6(1,1)) + c02(3).*(x-d6(1,1)).*(x-d6(2,1));
 plot(x,y02)
 
 % Aproximación logaritmica
 c03 = aprox_logaritmica(d6);
 y03 = c03(1).*exp(x.*c03(2));
 plot(x,real(y03))
 
 % exponencial de polinomio
 c04 = aprox_logaritmica_polinomios(d6);
 y04 = exp(c04(1) + c04(2).*x + c04(3).*(x.^2));
 plot(x,real(y04))
 
 legend({'Sistema de eq','Pol. de Newton','Logaritmica','exp de polinomio'},'Location','southwest')
 n = size(d6);
 for j = 1:n(1)
     plot(d6(j,1),d6(j,2),'r*')
 end
 hold off
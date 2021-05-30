% 
% ANALISIS DE ECUACIONES DIFERENCIALES
% SANTIAGO F.G. ZAMORA
% 1/5/2021
%==========================================================================

% PARAMETROS

% Intervalo de "tiempo"
T = 10;
% Número de subintervalos
n = 1000;

% SISTEMAS DE ECUACIONES DIFERENCIALES
% x' = Ax + Ac

% MATRIZ A
A = [8 -8;24 -6];
% Constantes
Ac = [0;0];
% CONDICION INICIAL
x0 = [1;0];

eps = 5.5511e-16; % para comparar con cero
%==========================================================================

lnA = length(A);                 % tamaño
[eig_vectA,eig_valA] = eig(A);   % Eigenvectores y Eigenvalores
detA = det(A);                   % Determinante
traA = trace(A);                 % Traza

% CLASIFICACION DE PUNTOS FIJOS
clas = '';
if detA ~= 0 && eigenv_distintos(eig_valA,lnA) && eigenv_complejos_pteReal_dist_cero(eig_valA,lnA)
    % NO DEGENERADO
    clas = strcat(clas,'no degenerado');
    b = zeros(lnA,1);
    pfs = A\b; % punto fijo único
    [est,caso] = estabilidad(eig_valA,lnA);
    if est
       clas = strcat(clas,', estable');
    else
        clas = strcat(clas,', inestable');
        switch caso
            case 1
                clas = strcat(clas,', respulsor');
            case 2
                clas = strcat(clas,', punto silla');
            case 3
                clas = strcat(clas,', atractor');
            otherwise
                clas = strcat(clas,', nada');
        end
    end
else
    % DEGENERADO
    clas = strcat(clas,'degenerado');
    pfs = zeros(lnA,1);
end
pfs
clas

% CLASIFICACION DE PUNTOS FIJOS 2
c_traA = (1/4)*traA^2;
dif_det_tra = detA - c_traA;
if igual(detA,0,eps) || igual(traA,0,eps) || igual(detA,c_traA,eps)
    clas2 = 'degenerado';
else
    if detA < 0
       clas2 = 'punto silla';
    else
        if traA > 0
            if dif_det_tra > 0
                clas2 = 'espiral respulsora';
            else
                clas2 = 'repulsor';
            end
        else
            if dif_det_tra > 0
                clas2 = 'espiral atractora';
            else
                clas2 = 'atractor';
            end
        end
    end
end
clas2

% Solución aproximada dado el punto inicial
dt = T/n;
x = zeros(lnA,n+1); % vector(es) solución
x(:,1) = x0;
for k = 2:n+1
   x(:,k) = x(:,k-1) + (A*x(:,k-1))*dt+Ac;
end

% GRAFICAS
% Solución
% figure('Name','Solución','NumberTitle','off')
% tiempo = dt*(0:n);
% plot(tiempo,x)
% hold on
% title('Soluciones por coordenada')
% xlabel('tiempo') 
% ylabel('x')
% yline(0);
% legend(labels(lnA),'Location','southwest')
% hold off

% Diagrama de fase
if lnA == 2
    figure('Name','Diagrama de fase','NumberTitle','off')
    title('Diagrama de fase')
    hold on
    xlabel('x') 
    ylabel('y')
    c = max(abs(x(1,:))) + 1;
    df = (2*c)/n;
    xx = -c:df:c;
    [yv,yh] = isoclinas(A,Ac,xx);
    if length(yv) == length(xx)
       plot(xx,yv,'r')
    else
        xline(yv,'r');
    end
    if length(yh) == length(xx)
       plot(xx,yh,'g')
    else
        xline(yh,'g');
    end
    plot(x(1,:),x(2,:),'b')
    plot(x0(1),x0(2),'o')
    plot(pfs(1),pfs(2),'*')
    xline(0);
    yline(0);
    legend({'Solución','xp = 0 (verticales)','yp = 0 (horizontales)','condición incial','punto fijo','eje x','eje y'},'Location','southwest')
    hold off
end

%                         FUNCIONES
%==========================================================================
function vpd = eigenv_distintos(eig_val,ln)
    % Verificar si hay eigenvalores iguales
    vpd = true; % valores propios distintos
    j = 1;
    while j<ln && vpd
        k = j+1;
        while k<=ln && vpd
            if eig_val(j,j) == eig_val(k,k)
               vpd = false; 
            end
            k = k+1;
        end
        j = j+1;
    end
end

function vpcd0 = eigenv_complejos_pteReal_dist_cero(eig_val,ln)
    % Verificar si hay valores propios comlejos con parte real igual a cero
    vpcd0 = true; % valores complejos con parte real distinta de cero [true]
    j = 1;
    while j <= ln && vpcd0
       if igual(real(eig_val(j,j)),0,eps) && abs(imag(eig_val(j,j))) >= eps
          vpcd0 = false; 
       end
       j = j+1;
    end
end

function [vpcrn,caso] = estabilidad(eig_val,ln)
    % Verificar si todos los eigenvalores tienen valor real negativo
    vpcrn = true; % valores propios complejos con parte real negativa
    cont = 0;
    for j = 1:ln
        if real(eig_val(j,j)) > 0
            if vpcrn
                vpcrn = false;
            end
            cont = cont+1;
        end
    end
    % CASOS
    if cont == ln % todos los valores propios tienen parte real positiva
       caso = 1; % repulsor
    else
        if cont == 0 % todos los valores propios tienen parte real negativa
            caso = 3; % atractor
        else 
            caso = 2; % punto silla
        end
    end
end

function res = igual(a,b,eps)
    res = false;
    if  abs(a-b)<=eps
        res = true;
    end
end

function str = labels(ln)
    str = string.empty;
    for j = 1:ln
        num = string(j);
        str(j) = strcat('x',num);
        if j == ln
            str(j+1) = 'eje x';
        end
    end
end

function [v,h] = isoclinas(mat,matc,tiempo)
    if mat(1,2) ~= 0
        v = -(mat(1,1).*tiempo + matc(1))./mat(1,2); 
    else
        v = -matc(1)/mat(1,1);
    end
    
    if mat(2,2) ~= 0
        h = -(mat(2,1).*tiempo + matc(2))./mat(2,2); 
    else
        h = -matc(2)/mat(2,1);
    end
end

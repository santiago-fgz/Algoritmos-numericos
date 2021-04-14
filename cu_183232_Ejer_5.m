% Ejercicio 5
% Santiago Fernández Gutiérrez Zamora
% 183232
%==========================================================================

% Aproximación de Newton

n = 4; % el número de datos a analizar por ciclo

datos = readtable('Ejer_5_datos.csv');
datos = table2array(datos);
m = length(datos); 

x = zeros(m,1);
y = zeros(m,1);
for j = 1:m-n
    datosActx = [datos(j:j+n-1,1),datos(j:j+n-1,2)]; % tomar cuatro datos [t,x]
    datosActy = [datos(j:j+n-1,1),datos(j:j+n-1,3)]; % tomar cuatro datos [t,y]
    cfsx = newton(datosActx,n); % coeficientes del polinomio
    cfsy = newton(datosActy,n);
    tpred = j+n; % tiempo para el que se predice el valor
    % polinomios:
    x(j+n) = cfsx(1) + cfsx(2)*(tpred-datosActx(1,1)) + cfsx(3)*(tpred-datosActx(1,1))*(tpred-datosActx(2,1)) + cfsx(4)*(tpred-datosActx(1,1))*(tpred-datosActx(2,1))*(tpred-datosActx(3,1));
    y(j+n) = cfsy(1) + cfsy(2)*(tpred-datosActy(1,1)) + cfsy(3)*(tpred-datosActy(1,1))*(tpred-datosActy(2,1)) + cfsy(4)*(tpred-datosActy(1,1))*(tpred-datosActy(2,1))*(tpred-datosActy(3,1));
end
errorx = abs(datos(:,2)-x(:)); % la diferencia entre los valores predecidos 
errory = abs(datos(:,3)-y(:)); % y los reales
errorx(1:n) = 0; % los primeros cuatro valores no se consideran
errory(1:n) = 0;
todo = [datos x y errorx errory];
titulos = ["t","x","y","predicción de x","predicción de y","error de x (absoluto)","error de y (absoluto)"];
respuesta = array2table(todo,'VariableNames',titulos)
meanErry = mean(errory(n+1:m)) % promedio del error
meanErrx = mean(errorx(n+1:m))

function m = pend(a,b)
% pendiente
    if a(1) ~= b(1)
       m = (b(2)-a(2))/(b(1)-a(1));
    end
end

function A = newton(dat,n)
    % entrega los coeficientes para formar el polinomio de newton
    % dat tiene que ser una matriz de n*2
    if n>1
       f = zeros(n);
       for i = 1:n
          f(1,i) = dat(i,2);
          if  i-1 > 0
              f(2,i) = pend(dat(i-1,:),dat(i,:));
              if i-2 > 0
                 f(3,i) =  pend([dat(i-2,1),f(2,i-1)],[dat(i,1),f(2,i)]);
                 if i-3 > 0
                    f(4,i) = pend([dat(i-3,1),f(3,i-1)],[dat(i,1),f(3,i)]);
                 end
              end
          end
       end
       A = diag(f);
    end
end
% Ejercicio 6
% Santiago Fernández Gutiérrez Zamora
% 183232
%==========================================================================

% Regresión lineal

tabla = readtable('Ejercicio_6.csv');

% la regresión:
rl = fitlm(tabla,"VentasMensuales~Hogs_1_AB+Hogs_2_Cp+Hogs_3_C+Hogs_4_Cm+Hogs_5_Dp+Hogs_6_D+Hogs_7_E")
rcuad1 = rl.Rsquared % error

% coeficientes de correlación
datos = table2array(tabla);
[n,m] = size(datos);
datos = datos(:,2:m-1); % la primera columna es solo el número de dato, y la última no interesa al análisis de colineal
m = m-2;
pMax = 0.075;

colineal = 1; % condición de salida
colscol = zeros(1,m); % indica las columnas que son colineales
while colineal==1
    if m>=2
        R = corrcoef(datos);
        mMax = 2;
        max = abs(R(1,mMax));
        for i=1:n
            for j=1:m
                if j>i
                    if abs(R(i,j))>max
                       mMax = j;
                       max = abs(R(i,j));
                    end
                end
            end
        end
        if max>pMax
            colscol(mMax) = 1;
            datos(:,mMax) = []; % elminiar la columna colineal
            m = m-1;
        else
            colineal = 0;
        end
    else
        colineal = 0;
    end
end

% una mejor regresión

str = 'VentasMensuales~';
if colscol(1)==0
   str = append(str,'Hogs_1_AB+'); 
end
if colscol(2)==0
   str = append(str,'Hogs_2_Cp+'); 
end
if colscol(3)==0
   str = append(str,'Hogs_3_C+'); 
end
if colscol(4)==0
   str = append(str,'Hogs_4_Cm+'); 
end
if colscol(5)==0
   str = append(str,'Hogs_5_Dp+'); 
end
if colscol(6)==0
   str = append(str,'Hogs_6_D+'); 
end
if colscol(7)==0
   str = append(str,'Hogs_7_E+'); 
end
str = str(1:end-1);
rl2 = fitlm(tabla,str) % la regresión sin colineales

rcuad2 = rl2.Rsquared

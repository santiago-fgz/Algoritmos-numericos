% Splines

% El spline S1 sale (t=0) de P1 y en tal punto es tangente a P0.
% S1 termina en P2 (t=1) y en tal punto es tangente a la línea LP2P3.
% El spline S2 sale de P2 (para t=1) y en tal punto es tangente a la línea LP2P3,
% S2 termina en P4 (para t=2) y en tal punto es tangente a la línea recta LP4P5

P0 = [1,5];
P1 = [2,6];
P2 = [3,4];
P3 = [2,3];
P4 = [4,2];
P5 = [6,3];

a = 0;
b = 6;
brk = 3;
eps = 0.01;

% spline 1

pasa_por = [P1;P2];
refs = [P0;P3];

s1 = spline_cub_trayect_r2(pasa_por,refs,0,1);
tt1 = 0:eps:1;
xx1 = s1(1,1) + s1(2,1).*tt1 + s1(3,1).*tt1.^2 + s1(4,1).*tt1.^3;
yy1 = s1(1,2) + s1(2,2).*tt1 + s1(3,2).*tt1.^2 + s1(4,2).*tt1.^3;

plot(xx1,yy1)
hold on

pasa_por = [P2;P4];
refs = [P3;P5];

s2 = spline_cub_trayect_r2(pasa_por,refs,1,2);
tt2 = 1:eps:2;
xx2 = s2(1,1) + s2(2,1).*tt2 + s2(3,1).*tt2.^2 + s2(4,1).*tt2.^3;
yy2 = s2(1,2) + s2(2,2).*tt2 + s2(3,2).*tt2.^2 + s2(4,2).*tt2.^3;
plot(xx2,yy2)

hold off

function spline_cs = spline_cub_trayect_r2(puntos,refs,a,n)
    spline_cs = zeros(4,2); % [cx,cy] constantes para x en la primera col
    if length(puntos) == 2 && length(refs) == 2
        b = a+1;
        % valores y sus derivadas
        T = [1, a, a^2, a^3;    ...
             0, 1, 2*a, 3*(a^2);...
             1, b, b^2, b^3;    ...
             0, 1, 2*b, 3*(b^2)];
        vx = [puntos(1,1);puntos(1,1)-refs(1,1);puntos(2,1);refs(2,1)-puntos(2,1)];
        vy = [puntos(1,2);puntos(1,2)-refs(1,2);puntos(2,2);refs(2,2)-puntos(2,2)];
        if n == 2
           vx(2) = vx(2)*-1; 
           vy(2) = vy(2)*-1;
        end
        cx = T\vx;
        cy = T\vy;
        spline_cs(:,1) = cx;
        spline_cs(:,2) = cy;
    end
end
% Circuito RLC 

R = 10; % Ohms
L = 2.5; % Henr
C = 0.5e-3; % Farads
Vs = 5; % Volts

v0 = 0;
vp0 = 0;

T = 3; % segundos
n = 10000; % veces

h = T/n; % intervalo
tiempo = (0:n)*h; 

% Solución numérica

unoentreLC = 1/(L*C);
RentreL = R/L;
A = [0,1;-unoentreLC,-RentreL];
x = zeros(2,n+1);
vext = [0;Vs];
x(:,1) = [v0,vp0];
for t = 2:n+1
   x(:,t) = x(:,t-1) + (A*x(:,t-1)+vext)*h; 
end

voltaje = x(1,:);

% Solución teórica

alfa = R/(2*L);
w0 = 1/sqrt(L*C);
if  w0>alfa % Subamortiguada
    wd = sqrt(w0^2 - alfa^2);
    A1 = v0-Vs;
    A2 = -vp0/(alfa*A1);
    volt_teo = Vs + (A1*cos(wd*tiempo) + A2*sin(wd*tiempo)).*exp(-alfa*tiempo);
    %volt_teo = Vs + (A1*cos(wd*tiempo)).*exp(-alfa*tiempo);
end

%plot(tiempo,volt_teo,tiempo,voltaje)
%hold on
plot(tiempo,voltaje)
%hold off
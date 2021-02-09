B = 3;
tt = 30;
cmax = 0.25;
cmin = (1/4)*cmax;
res = zeros;
t = 1:tt;
%
cont = 1;
A = (cmax*exp(1))/B;
%
c1 = @(x,a,b) a.*x.*exp(-(x./b));
ct1 = @(x) c1(x,A,B);
%
tmax = B;
tf1 = fsecante(ct1,cmin,B,tt,0.01,100);
%
res(1,1:4) = [cont,A,0,tmax];
%
cont = cont+1;
%
c2 = @(x,a,b,x2,a2) c1(x,a,b) + c1(x-x2,a2,b);
cta2 = @(x,a2) c2(x,A,B,tf1,a2);
taprox = tf1+B;
aaprox = (2/3)*A;

while abs(cta2(taprox,aaprox)-cmax) > 0.01
   %
   ca2 = @(x) cta2(taprox,x);
   aaprox = fsecante(ca2,cmax,A/2,A,0.01,100);
   %
   ct2 = @(x) cta2(x,aaprox);
   ct2s = sym(ct2);
   ct2ps = diff(ct2s);
   ct2p = matlabFunction(ct2ps);
   taprox = fsecante(ct2p,0,tf1,tt,0.001,100);
end
A2 = aaprox;
tmax2 = taprox;
%
res(2,1:4) = [cont,A2,tf1,taprox];
%
c2l = @(x,a,b,x2,a2) c1(x,a,b) + (x>x2)*c1(x-x2,a2,b);
ct2l = @(x) c2l(x,A,B,tf1,A2);
csim = zeros;
for i = 1:tt
    csim(i) = ct2l(i);
end
plot(t,csim(t))
xlabel('tiempo desde la primera dosis [hrs]')
ylabel('concentración en la sangre [mg/mm]')
yline(cmin,'--','concentración mínima')
yline(cmax,'--','concentración máxima')
xline(tf1,'--','2da dosis')
title('Concentración de la droga en la sangre')

titulos = ["dosis N°","dosis [mg]","tiempo de aplicación [hrs]","tiempo de máximo efecto [hrs]"];
T = array2table(res,'VariableNames',titulos)




function T = frecuencias(L,E,n)
    %
    % Determinar las dos frecuencias más bajas de una viga en centilever
    % 
    % Ejemplo
    % res = frecuencias(0.9,200,2)
    %
    % T: frecuencias [Hz]
    %
    % L: longitud [m]
    % E: módulo de elasticidad [GPa]
    % n: número de frecuencias [entero]
    % rango: 
    %
    T = zeros;
    if L > 0 && E > 0
        K = (48*(pi^2)*L)/E;
        beta = @(x) ((x.*x)*K)^(1/4);
        fct = @(x) cosh(beta(x)).*cos(beta(x)) - 1;
        T = zeros; 
        rs = 0;
        i = 1;
        while rs < n 
            if fct(i)*fct(i+1) <= 0
               r = fSecante(fct,0,i,i+1,0.01,100);
               rs = rs + 1;
               T(rs) = r;
            end
            i = i+1;
        end
    end
end

function x = fSecante(f,alfa,a,b,eps,mxv)
  fa = f(a) - alfa;
  fb = f(b) - alfa;
  if  fa * fb > 0 || eps < 0 || mxv < 1
      x = nan;
  else
      i = 1;
      while i <= mxv
        xm = -(fa*((b-a)/(fb-fa)))+a;
        fxm = f(xm) - alfa;
        if abs(fxm) < eps
            i = mxv;
        else
            if fxm*fb > 0
                b = xm;
                fb = fxm;
            else
                a = xm;
                fa = fxm;
            end
        end
        i = i + 1;
      end
    x = xm; 
  end
end

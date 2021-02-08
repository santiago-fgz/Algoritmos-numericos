function T = anestesia(w,n,cmax,cmin,b)
    %  
    % 
    % determinar la dosis de una droga asi como los tiempos en los que se
    % debe aplicar otra(s) dosis.
    %
    % Ejemplo:
    % res = anestesia(70,60,0.25,0.0625,3)
    %
    % T: [ cont , dosis , tap , tmax]
    %   cont: número de dosis
    %   dosis: cantidad de droga por dosis [mg]
    %   tap: tiempo en el que  se debe aplicar la dosis [hrs]
    %   tmax: tiempo en el que se alcanza el mayor efecto [hrs]
    %
    % w: peso del paciente [kgs]
    % n: tiempo que el paciente debe permanecer bajo efecto [hrs]
    % cmax: concentración máxima [mg/cc]
    % cmin: concentración mínima [mg/cc]
    % b: tiempo que tarda la droga en llegar a su máximo efecto despues de
    %    la primera inyección [hrs]
    %
    % Santiago F.G. Zamora
    
    titulos = ["dosis N°","dosis [mg]","tiempo de aplicación [hrs]","tiempo de máximo efecto [hrs]"];
    if n > 0 && b > 0 && cmax > 0 && cmin > 0 && cmin < cmax
        dos = (cmax*exp(1))/b;
        tmax = b;
        cont = 1;
        tac = 0;
        A(1,1:4) = [cont,dos,tac,tmax];
        cc = @(t,D) t.*exp(-(t./b)).*D;
        x = 0:n;
        plot(x,cc(x,dos))
        xlabel('tiempo desde la primera dosis [hrs]')
        ylabel('concentración en la sangre [mg/mm]')
        yline(cmin,'--','concentración mínima')
        yline(cmax,'--','concentración máxima')
        title('Concentración de la droga en la sangre') 
        hold on
        tac = funcSecanter2(cc,cmin,b,n,0.01,1000,dos);
        while tac < n
            dos = ((cmax-cmin)*exp(1))/b;
            tmax = tac + b;
            tsig = funcSecanter2(cc,cmin,b,n,0.01,1000,dos);
            extra = cc(tsig + tac,A(cont,2));
            cont = cont + 1;
            if tac < n
                A(cont,1:4) = [cont,dos,tac,tmax];
                plot(x+tac,cc(x,dos))
            end
            tac = funcSecanter2(cc,cmin-extra,b,n,0.01,1000,dos) + tac;
        end
        hold off
        T = array2table(A,'VariableNames',titulos);
    end
end

function x = funcSecanter2(f,alfa,a,b,eps,mxv,c)
  %
  % Obtiene raices de funciones mediante el método de la secante
  %
  % f: handler de la función
  % alfa: valor deaseado de f
  % a y b : límites
  % eps: error del resultado
  % mxv: máximo de vueltas
  
  fa = f(a,c) - alfa;
  fb = f(b,c) - alfa;
  if  fa * fb > 0 || eps < 0 || mxv < 1
      x = nan;
  else
      i = 1;
      while i <= mxv
        xm = -(fa*((b-a)/(fb-fa)))+a;
        fxm = f(xm,c) - alfa;
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
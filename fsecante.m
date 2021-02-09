function x = fsecante(f,alfa,a,b,eps,mxv)
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

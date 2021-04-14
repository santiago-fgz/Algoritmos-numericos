function coefs = Tch(hf,n,a,b)
    if n>1 && a<b
        k = 0:n-1;
        k = fliplr(k);
        w = pi/(2*n);
        theta_k = (2*k + 1)*w;
        uk = cos(theta_k);
        
        MT = ones(n);
        MT(:,2) = uk;
        for j = 3:n
            for r = 1:n
                MT(r,j) = 2*(uk(r)*MT(r,j-1)) - MT(r,j-2); 
            end
        end
        
        xk = ((b-a)/2)*uk + (a+b)/2;
        y = hf(xk);
        coefs = y*0;
        mult = 1;
        for l = 1:n
            coefs(l) = (sum(MT(:,l).*y')/n)*mult;
            if l == 1
                mult = 2;
            end
        end
        
    end
end
function res = aprox_logaritmica_polinomios(d)
    [n,m] = size(d);

    if n>1 && m==2
        x = zeros(n);
        y = zeros(n,1);
        for i = 1:n
            xi = d(i,1);
            y(i) = log(d(i,2));
            for j = 1:n
                x(i,j) = xi^(j-1);
            end
        end
        res = x\y;
    else
        res = 0;
    end
end
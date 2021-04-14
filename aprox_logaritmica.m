function res = aprox_logaritmica(d)
    % A*exp(B*x) = y
    % ln(A)+(B*x) = ln(y)
    % Ap + B*x = yp
    % (1,x)*(Ap;B) = yp
    % X * (Ap;B) = yp
    % A = exp(Ap)
    % y = exp(yp)
    [n,m]=size(d);
    
    as = zeros(n-1,1);
    bs = zeros(n-1,1);
    for j = 1:n-1
        X = ones(2);
        X(:,2) = [d(j);d(j+1)];
    
        y = [d(j,2);d(j+1,2)];
        yp = log(y);
        
        apb = X\yp;
        as(j) = exp(apb(1));
        bs(j) = apb(2);
    end
    res = [mean(as),mean(bs)];
end

function an = determinaAncho(a,b,h)
    %
    % Ejercicio 1
    %
    % Determinar el ancho de dos paredes paralelas dados
    % los tamaños de dos escaleras recargadas sobre las 
    % paredes y que se intersectan a una altura dada sobre
    % suelo
    %
    % an: distancia entre ambas paredes
    % 
    % a: tamaño de la escalera 1
    % b: tamaño de la escalera 2
    % h: altura a la que se intersectan ambas escaleras
    %
    % Santiago F.G. Zamora
    
    if a>0 && b>0 && h>0
        fun = @(x) (h*((1/sqrt(a^2 - x.*x))+(1/sqrt(b^2 - x.*x))))-1;
        %an = raizPorNewton(fun,min(a,b)-0.01,0.01,100);
        an = fsecante(fun,0,min(a,b)*(0.9),min(a,b)-0.01,0.01,100);
    else
        an = -1;
    end
end
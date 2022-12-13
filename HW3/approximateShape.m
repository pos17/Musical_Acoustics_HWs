function [coneX,coneY] =  approximateShape(a0,m,nPoints,delta)
    coneX = zeros(nPoints,1); 
    coneY = zeros(nPoints,1);
    coneX(1) = 0;
    coneY(1) = a0;
    for n=2:nPoints
        coneX(n) = (n-1)*delta;
        coneY(n) = a0*exp(m*coneX(n));
    end
end 
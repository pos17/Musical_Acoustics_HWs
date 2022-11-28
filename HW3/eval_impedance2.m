function Zin = eval_impedance2(Npoints, L, a0, rho, c, k)
    Zin = zeros(1,length(k));
    x = linspace(0,L,Npoints);
    m = 4.2;
    a = a0.*exp(m*x);

    a=a0*exp(m*L);
    w = k*c;
    ZL0 = 0.25*(w.^2*rho)/(c*pi) + 0.61*1i*(rho*w)/(pi*a);
    Sp = a^2*pi;
    
    ii = length(x);
    
    if ii>1 
        l = x(2)-x(1);
        while ii>1
            a2 = a(ii);
            a1 = a(ii-1);
            S2 = (a2)^2*pi;
            S1 = (a1)^2*pi;
        
            phi = (a2-a1)/l;
        
            x20 = a2/phi;
            x10 = a1/phi;
            
            theta2 = atan(k.*x20);
            theta1 = atan(k.*x10);

            if ii==length(x)
                theta_pr = atan(a2/x20);
                Ss = 2*Sp/(1+cos(theta_pr));
                ZL = ZL0.*Sp./Ss;
                Z = ZL;
            end
        
            num = 1i*Z.*(sin(k*l-theta2)./sin(theta2)) + rho*c/S2 .* sin(k*l);
            den = Z.*(sin(k*l+theta1-theta2)./(sin(theta1).*sin(theta2))) - 1i*rho*c/S2*(sin(k*l+theta1)./sin(theta1));
            
            Zin = rho*c/S1 * (num./den);
        
            Z = Zin;
            ii=ii-1;
        end
    else
        l = x;
        a2 = a0*exp(m*L);
        a1 = a0;
        phi = (a2-a1)/L;
        x10 = a1/phi;
        theta1 = atan(k.*x10);
        S1 = a0^2*pi;
        num = sin(k*L).*sin(theta1);
        den = sin(k*l+theta1);
        Zin = 1i*rho*c/S1*num./den;
    end

end
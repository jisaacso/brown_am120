function MonteCarlo

N = 1000;    %number of spatial points to evaluate at
M = 1000;    %number of monte carlo sample points per spatial point
x = linspace(0,2*pi,N);
u2 = zeros(1,N);
for(i=1:N)
    integrand = @(y)(exp(-abs(x(i)-y).^2./(2*pi)).*sin(y));
    u(i) = 1/(pi*sqrt(2)) * quadgk(integrand,-inf,inf);
    for(j=1:M)
        r = sqrt(pi).*randn(1,1);
        u2(i) = u2(i)+sin(x(i)+r)/M;
    end
end


h2 = plot(x,u2,'r');
hold on;
h1 = plot(x,u,'k');
set(h1,'LineWidth',5)
legend('Monte Carlo','Quadrature');
drawnow();
end
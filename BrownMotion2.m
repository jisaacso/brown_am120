function BrownMotion2

N = 1000;
M = 800;

t = linspace(0,2*pi,N);
Z = randn(M+1,1);
w = zeros(size(t));
w = Z(1)*(2*pi)^(-1/2).*t;

for(i = 1:M)
    w = w+2*pi^(-1/2)*i^(-1)*Z(i).*sin((i/2).*t);
end

plot(t,w,'k-.');

title(['Brownian Motion Trajectories M = ',num2str(M)]);
xlabel('t');
ylabel('W(t)');
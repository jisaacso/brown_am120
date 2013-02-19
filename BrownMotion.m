function BrownMotion

N = 1000;   %Number of time steps between 0 and 2pi
M = 800;    %Number of RW Steps

hold on;
for(numIts = 1:100)
    t = linspace(0,2*pi,N);
    Z = randn(M+1,1);
    w = zeros(size(t));
    w = Z(1)*(2*pi)^(-1/2).*t;

    for(i = 1:M)
        w = w+2*pi^(-1/2)*i^(-1)*Z(i).*sin((i/2).*t);
    end

    plot(t,w,'k-.');
end

%Construct 2nd order polynomial p with
%asymptotics equivalent to law of iterated logs
x = [0,t([2:9,N-7:N])];
y(1) = 0;
y(2:9) = sqrt(2*t(2:9).*log(log(1./t(2:9))));
y(10:17) = sqrt(2*t(N-7:N).*log(log(t(N-7:N))));
p = polyfit(x,y,2);

h1=plot(t,p(3)+p(2)*t+p(1)*t.*t,'b');
h2=plot(t,-p(3)-p(2)*t-p(1)*t.*t,'b');
set(h1,'LineWidth',8);
set(h2,'LineWidth',8);

title(['Brownian Motion Trajectories M = ',num2str(M)]);
xlabel('t');
ylabel('W(t)');
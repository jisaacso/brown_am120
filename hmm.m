function piN=hmm

data = [2.2569 2.3975 2.9830 -2.5518 2.8655 -3.3622 -2.4127 0.3327 -1.7862 2.4159 3.1662 ...
2.8880 2.6574 3.2255 2.6949 -3.3070 -2.3132 2.1371 2.8722 2.6642 3.0503 -1.9723 ...
-2.1519 3.4754 2.6399 -2.5641 2.1840 2.0443 2.7239 0.9758 -2.2464 0.6123 2.0927...
-2.5792 -3.1641 -3.2582 -2.0708 2.5494 -3.0232 -3.0280 3.2116 -2.7659 -0.1542 2.3244...
-2.0709 0.0271 2.5650 -2.8149 2.4351 3.1217 -2.6108];
size(data)

P = [.5,.2,.1,.1,.1;.1,.6,.1,.1,.1;.1,.2,.4,.2,.1;.1,.2,.1,.5,.1;.1,.1,.1,.3,.4];

pi0 = [.5,.3,.1,.1,0];

psi = zeros(51,5);
phi = zeros(size(psi));
sumphi = zeros(size(phi));
piN = zeros(size(phi));
for(i=1:51)
    for(j=1:5)
        psi(i,j) = 1/sqrt(2*pi)*exp(-(data(i)-3*sin(j))^2/2);
    end
end

for(i=1:5)
    integrand = @(y)(exp(-(y-3.*sin(i)).^2./(2)));
    phi(1,i) = pi0(i)/sqrt(2*pi) * quadgk(integrand,0,inf);
end



for(i=1:50)
    for(j=1:5)
        for(k=1:5)
            phi(i+1,j) = phi(i+1,j)+psi(i+1,j)*phi(i,k)*P(k,j);
        end
    end
end

for(i=1:51)
    for(j=1:5)
        for(k=1:5)
            sumphi(i,j) = sumphi(i,j)+phi(i,k);
        end
        piN(i,j) = phi(i,j)/sumphi(i,j);
    end
end
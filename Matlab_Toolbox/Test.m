clc; clear
Tm=100;
d = zeros(1, 100);
for k = 1:Tm
    if k<= 1/3*Tm
        d(k) = 0.1;
    else
        d(k) = 0.1 + 0.4*tanh((18*k-6*Tm)/Tm);
    end
end

figure(1)
plot(d)
ylim([0, 0.6])

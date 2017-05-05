N=[10 20 40 60 80 100 120 140 160 200 320];

NormU=[288 328 413 485 548 605 657 705 749.5 832 1.0417e3];

n=length(N);
for i=2:n-1
    NormDiff(i) = NormU(i+1) - NormU(i);
end

plot(N(1:10),NormDiff,'*')
figure(1)
xlabel('Nodes')
ylabel('Norm of U difference')
for i=1:10
Error(i)=abs(NormDiff(i)-NormDiff(8))/(NormDiff(8));
end

figure(2)
plot(N(2:10),Error(2:10),'*')

xlabel('Nodes')
ylabel('Error')
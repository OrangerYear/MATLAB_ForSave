% kaczmarz_main

clc
clear all

%n=6,K=100;
n=100;K=1000;
[A1,b1,x1]=phillips(n);
x0(1:n,1)=0;
normx(1:K)=0;
% use the subroutine kaczmarz method
[x_solution1,error]=kaczmarz(A1,b1,x0,K,x1);

nerror(1:K)=0;
for i=1:100
     [x_solution1,error1]=kaczmarz1(A1,b1,x0,K,x1);
     nerror=nerror+error1;
end
error1=nerror/100;
plot(error,'b');
hold on
plot(error1,'r');
legend('kaczmarz','randomized kaczmarz')
set(gca,'yscale','log')



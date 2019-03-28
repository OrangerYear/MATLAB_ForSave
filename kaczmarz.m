function [x,error]=kaczmarz(A,b,x0,Kmax,x_true)
   %  kaczmarz is the program of standard kaczmarz method
   %  usage: [x,error]=kaczmarz(A,b,x0,Kmax,x_true) 
   %  A,b,  the coefficient matrix and right-hand vector of the linear system;
   %  Kmax, the maximum iterative number;
   %  x_ture, the true solution of Ax=b.
   %
   %  by Chuan-gang Kang,september 29,2018.
  
   if (nargin<4),x0(1:size(A,2),1)=0; end
   if (nargin<5),Kmax=10000;end   % in the use of function, the x_true version is necessary.

   % xt is true solution
   [nr,nc]=size(A);
   nb=length(b);
   delta=0.1;
   %delta=0;
   tau=1.1;   
   % generate the disturbed right-hand
   xi=randn(nb,1);  % xi is randomized vector and the its dimension is in consistence with right-hand b 
   b_delta=b+delta*max(abs(b));%*xi;

   xk=x0;
   for k=1:Kmax
       if mod(k,nr)==0
           i=nr;
       else
           i=mod(k,nr);
       end
       % compute the norm \|a_i\|^2
       noi=sum(A(i,:).^2);
       % update x_k 
       xk=xk+(b_delta(i)-A(i,:)*xk)*A(i,:)'/noi;
       % compute the error between xk and x_true
       error(k)=norm(xk-x_true)^2;
   end
   x=xk;
end
  


   
   

% the prgrame of randomized kaczmarz1 method
% by Chuan-gang Kang,september 29,2018.

function [x,error]=kaczmarz1(A,b,x0,Kmax,x_true)
   % A,b are coefficient matrix and right-hand vector of the linear system;
   % Kmax is the maximum iterative number;
   % x_ture is the true solution of Ax=b.
   
   if (nargin<4),x0(1:size(A,2),1)=0; end
   if (nargin<5),Kmax=10000;end   % in the use of function, the x_true version is necessary.

   % xt is true solution
   [nr,nc]=size(A);
   nb=length(b);
%    delta=0.1;
   delta=0;
   tau=1.1;   
   % generate the disturbed right-hand
   xi=randn(nb,1);  % xi is randomized vector and the its dimension is in consistence with right-hand b 
   b_delta=b+delta*max(abs(b));%*xi;
   
   % define the row range of random choice
   rorr=[1:nr]'; % rorr is the abbr. range of choice of random row. 

   % compute norm(a_i)^2 and sum(a_ij^2);
   nof=0;  % nof returns the norm of frobenius of coefficient matrix;
   nos(1:nr,1)=0;
   for i=1:nr
       nos(i)=sum(A(i,1:nc).^2); % nos(i) is the square sum of the norm vector of the ith equation;
       nof=nof+nos(i);
   end

   % returns the ratio of ||a_i||^2/||A||_F
   prob=nos/nof;   % prob is, in fact, the probabilty to choose the equation
 
   xk=x0;
   for k=1:Kmax
       ni=random_number(rorr,prob); 
       % compute the norm \|a_i\|^2
       noi=sum(A(ni,:).^2);
       % update x_k 
       xk=xk+(b_delta(ni)-A(ni,:)*xk)*A(ni,:)'/noi;
       % compute the error between xk and x_true
       error(k)=norm(xk-x_true)^2;
   end
   x=xk;
end
  


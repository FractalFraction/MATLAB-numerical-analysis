function [hess,U]=Hessenberg(A)
'Algorithm that returns a hessenberg form of a matrix with its unitary matrix' 
[n,m]=size(A);
U = eye(m);

for j=1:n-2
    H=ModifiedHouseholder(A,j,n);
    A=(H*A)*H;
    U = H*U;
  
hess = A;

end 
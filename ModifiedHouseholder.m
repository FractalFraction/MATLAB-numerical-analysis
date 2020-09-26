function h = ModifiedHouseholder(a,j,n)
'A modified version of the Householder algorithm, used to bring matrices into Hessenberg form'


v = a(:,j);
v(1:j) =[];
v(1) = v(1) - norm(v);
h = eye(n); 
h(j+1:n,j+1:n) = eye(n-j) - (2*eye(n-j)*(v*v')/(v'*v));   
    
end
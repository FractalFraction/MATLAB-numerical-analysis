function [Q,R] = Householder(a)
n = size(a,1);
m = size(a,2);
%Check the matrix is full column rank.
if n >= m
    q = eye(n);

    for j = 1:m
        h = HouseholderStep(a,j,n);
        a = h*a;
        q = h*q;
    end
    Q = q';
    R = a;
    
end

function h = HouseholderStep(a,j,n)
v = a(:,j);
    if j==1
        v(1) = v(1) - norm(v);
        h = eye(n) - (2*(v*v')/(v'*v));

    else
        v(1:j-1) =[];
        v(1) = v(1) - norm(v);
        h = eye(n); 
        h(j:n,j:n) = eye(n+1-j) - (2*eye(n+1-j)*(v*v')/(v'*v));   
    end
    
    end

end
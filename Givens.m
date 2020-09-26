function [Q,R] = Givens(x)
%The main function implementing Given's Rotations
n = size(x,1);
m = size(x,2);
%Check if matrix is of full column rank
if n>=m
    G = eye(n);
    %Gives number of rows and columns of matrix x
    for j = 1:m
    %Go along columns left to right
        for i = n:-1:n-(m - j)
        %Go up the rows, from the bottom and stop at the diagonal position
        g = GivensDeterminer(x,i,j,n,m)
        G = g*G        
        %Cumulative Givens rotation matrices        
        x = g*x
        end
    end
    %Finding Q,R from the Givens Matrices. 

end
    R = x;
    Q = transpose(G);
end

function g = GivensDeterminer(x,i,j,n,m)
%For the algorithm, we take val2 to be the value of the current position in
%the matrix. val1 is the value of thr diagonal entry of the same column.
%Determines if a Givens Matrix is needed. If it is, then this determines the axis of rotation.
val2 = x(i,j);
if val2 == 0    
    g = eye(n);

else 
   val1 = x(n-(m - j +1),j);
   [p,q] = GivensAngle(val1,val2); 
   g = eye(n);
   %Inserting the rotation component of a Given's matrix. Always shifts
   %weight to the diagonal element - avoid calculating unnecessary Given's
   %Matrices
   g(n-(m - j +1),n-(m-j + 1)) = p;
   g(i,n-(m-j + 1)) = q;
   g(n-(m-j+1),i) = -q;
   g(i,i) = p;
end

function [c,s] = GivensAngle(val1,val2)
%Calculates the cos and sin of the angle needed for a given step of the Given's Rotations   
    a = [val1, -val2; val2, val1];
    b = [sqrt((val1)^2 + (val2)^2); 0];
    y = linsolve(a,b);
    c = y(1); 
    s = y(2);
end

end
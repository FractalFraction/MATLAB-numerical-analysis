function [Q,R] = GramSchmidt(x)
%A function calculating reduced QR decompostion by Clasical Gram
%Schmidt Approach
%Finding size of matrix
n = size(x,1);
m = size(x,2);

Q = zeros(n,m);
R = zeros(m);

%Working out the 1st column of the algorithm
v = x(:,1);
no = norm(v);

%Find q1
R(1,1) = no;
Q(:,1) = v/(no);

    for j = 2:m
        a = x(:,j);
        v = a;
        r_list = zeros(m,1);
        
            for i = 1:j-1
                %Finding v and column of r for each step.
                r_list(i) = (Q(:,i)'*a); 
                v = v - (r_list(i))*Q(:,i)               
                
        no = norm(v);
        %calculating diagonal entry of r column
        r_list(j) = no; 
        %calculating 
        Q(:,j) = v/(no); 
        R(:,j) = r_list;
        end 
    end
    Q = Q
    R = R
end
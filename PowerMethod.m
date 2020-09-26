function [q,evals]= PowerMethod(A,N)
q=A(:,1);
q=q/norm(q);
evals=zeros(1,N);

for i=1:N
    q=A*q;
    lambda=(q'*(A*q))/(q'*q);
    evals(i)=lambda;
end
end
            
       
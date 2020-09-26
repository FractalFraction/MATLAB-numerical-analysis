function [q,evals]= PM_Normal(A,N)
q=A(:,1);
q=q/norm(q);
evals=zeros(1,N);

for i=1:N
    p=A*q;
    q=p/norm(p);
    lambda=(q'*(A*q))/(q'*q);
    evals(i)=lambda;
end
end
function [Q,D]=myQRbasic(A,N)
for i=1:N
    [Q,R]=qr(A);
    A=R*Q;

D=A;
Q=Q;
end

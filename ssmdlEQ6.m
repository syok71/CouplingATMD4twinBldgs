% function program for constitution of state equation to make A,B,C,D matrix
% for two buildings
% for earthquake loading
% for n-ATMDs between 2 Buildings


function [A0,Bw,C0,D0]=ssmdlEQ6(ndof1,ndof2,n,iM,C,K)

N=ndof1+ndof2+n;
O1=zeros(N);
O3=zeros(N,1);
I=eye(N);
A0=[O1 I;-iM*K -iM*C];
Bw=[O3;-ones(N,1)];
D0=zeros(ndof1+ndof2,1);
C0=[eye(ndof1) zeros(ndof1,1+ndof2);zeros(ndof2,ndof1) eye(ndof2) zeros(ndof2,1)];
for j=2:ndof1
   C0(j,j-1)=-1;
end
for j=ndof1+2:ndof1+ndof2
   C0(j,j-1)=-1;
end
C0=[C0 zeros(ndof1+ndof2,N)]; % Relative Story Drift
% This program calculates natural frequencies accorind to given mass and stiffness matrix.
% Mass and stiffness matrices don't need to be symmetric!
% Eigen-vectors are normalized in the way which all the diagonal terms of the madal mass matrix are equal to unit value.

function [w,f,Eivt]=eigkm(K,M)

[evt,ww]=eig(K,M);
ww=diag(ww);
[ww2,index]=sort(ww);
n=max(size(ww));
Eivt=zeros(n,n);
for jj=1:n
   Eivt(:,jj)=evt(:,index(jj));
   sc=Eivt(:,jj)'*M*Eivt(:,jj);
   Eivt(:,jj)=Eivt(:,jj)/sqrt(sc);
   if Eivt(1,jj)<0
       Eivt(:,jj)=-Eivt(:,jj);
   end
end
w=sqrt(ww2); % rad/sec
f=w/2/pi; % Hertz
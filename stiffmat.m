% Construct Stiffness Matrix
function K=stiffmat(k)

n=length(k);
K=zeros(n,n);
K(n,n)=k(n);
for j=1:n-1
   K(j,j+1)=-k(j+1);
   K(j+1,j)=-k(j+1);
   K(j,j)=k(j)+k(j+1);
end
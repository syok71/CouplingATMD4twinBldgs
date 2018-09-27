function  [a0,a1,xi,C]=rayDamp2(wn,ii,jj,xii,xjj,M,K)

A=[1/wn(ii) wn(ii)
   1/wn(jj) wn(jj)]*0.5;
sol_1=inv(A)*[xii;xjj];
a0=sol_1(1,1);
a1=sol_1(2,1);
n=length(wn);
for j=1:n
   xi(j,1)=0.5*a0/wn(j)+a1*wn(j)*0.5;
end
C=a0*M+a1*K;

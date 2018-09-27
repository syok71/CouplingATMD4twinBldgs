function [Bu_btwA,LocAtmd]=bubwtatmd(ATMDFLoor,ndof1,ndof2,im1ATMDm2)

N=length(ATMDFLoor);
for j=1:N
    if ATMDFLoor(j) > ndof1 || ATMDFLoor(j) > ndof2
        disp('Link Floor Error')
        return
    end
end

LocAtmd=zeros(ndof1+ndof2+N,2*N);

for j=1:N
    LocAtmd(ATMDFLoor(j),2*j-1)=1;
    LocAtmd(ATMDFLoor(j)+ndof1,2*j)=1;
    LocAtmd(ndof1+ndof2+j,2*j-1)=-1;
    LocAtmd(ndof1+ndof2+j,2*j)=-1;
end    

Bu_btwA=[zeros(size(LocAtmd)); im1ATMDm2*LocAtmd];
function xdecoded=decoding001(G1,Nvar,BitArray,xmax,xmin)

id0=0;
xdecoded=zeros(Nvar,1);

for ii=1:Nvar
    bb0=id0+1;
    id0=BitArray(ii)+id0;
    Bj=bin2dec(G1(bb0:id0));
    Nbit2=BitArray(ii);
    xdecoded(ii,1)=xmin(ii)+Bj*(xmax(ii)-xmin(ii))/(2^Nbit2-1);
end
    
function Ka1=addatmaKC(N,K0,LFr,n1,n2,k1L,k1R)

% Construct K, C matrix for ATMD-Linked-two-Building

if LFr > n1 || LFr > n2
    disp('Link Floor Error')
    return
end
    
Ka1=addatmd0(N,K0,1,k1L+k1R);
Ka1(LFr,LFr)=Ka1(LFr,LFr)+k1L;
Ka1(n1+LFr,n1+LFr)=Ka1(n1+LFr,n1+LFr)+k1R;
Ka1(N+1,LFr)=Ka1(N+1,LFr)-k1L;
Ka1(LFr,N+1)=Ka1(LFr,N+1)-k1L;
Ka1(N+1,n1+LFr)=Ka1(N+1,n1+LFr)-k1R;
Ka1(n1+LFr,N+1)=Ka1(n1+LFr,N+1)-k1R;
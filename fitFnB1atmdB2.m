function fxj=fitFnB1atmdB2(OBJModels,Nch,G1,Nvar,BitArray,xmax,xmin)

% Load Earthquake Model
t1=OBJModels.t1;t2=OBJModels.t2;t3=OBJModels.t3;
Xg1=OBJModels.Xg1;Xg2=OBJModels.Xg2;Xg3=OBJModels.Xg3;

% Structure Model
ndof1=OBJModels.ndof1;
ndof2=OBJModels.ndof2;
%K1=OBJModels.K1;
%C1=OBJModels.C1;
%K2=OBJModels.K2;
%C2=OBJModels.C2;
Ndof=OBJModels.Ndof;
Q=OBJModels.Q;
iM1atmdM2=OBJModels.iM1atmdM2;
C0_12=OBJModels.C0_12;
K0_12=OBJModels.K0_12;
n=OBJModels.n;

% Constraints
Ulimit=OBJModels.Ulimit;
D1limit=OBJModels.D1limit;% meter
D2limit=OBJModels.D2limit;%meter
ZB1limit=OBJModels.ZB1limit;% meter
ZB2limit=OBJModels.ZB2limit;%meter

xj=zeros(Nch,1);
fxj=xj;

for jn3=1:Nch
    
    xdecoded=decoding001(G1(jn3,:),Nvar,BitArray,xmax,xmin);% Decoding process
    
    ATMDFLoor=xdecoded(1);% Decoding process
    k1L=xdecoded(2);% Decoding process
    c1L=xdecoded(3);% Decoding process
    k1R=xdecoded(4);% Decoding process
    c1R=xdecoded(5);% Decoding process
    Ra1=xdecoded(6);% Decoding process
    Ra2=xdecoded(7);% Decoding process
    Re=3+xdecoded(8);% Decoding process
    r1=Ra1*10^-Re;
    r2=Ra2*10^-Re;
    R=[r1 0;0 r2];
    
    Ka1=addatmaKC(Ndof,K0_12,ATMDFLoor,ndof1,ndof2,k1L,k1R);
    Ca1=addatmaKC(Ndof,C0_12,ATMDFLoor,ndof1,ndof2,c1L,c1R);

    [Aa1,Bwa1,C0a1,D0a1]=ssmdlEQ6(ndof1,ndof2,n,iM1atmdM2,Ca1,Ka1);
    [Bu_a1,LocAtmd_a1]=bubwtatmd(ATMDFLoor,ndof1,ndof2,iM1atmdM2);
    Gatmd=lqr(Aa1,Bu_a1,Q,R);
    LQRsys=ss(Aa1-Bu_a1*Gatmd,Bwa1,C0a1,D0a1);
    
    % Simulation
    [Yatmd1,t1,Xstate1]=lsim(LQRsys,Xg1,t1);
    [Yatmd2,t2,Xstate2]=lsim(LQRsys,Xg2,t2);
    [Yatmd3,t3,Xstate3]=lsim(LQRsys,Xg3,t3);
    ZB1maxEQ1=max(abs(Xstate1(:,1:ndof1)))';
    ZB1maxEQ2=max(abs(Xstate2(:,1:ndof1)))';
    ZB1maxEQ3=max(abs(Xstate3(:,1:ndof1)))';
    ZB2maxEQ1=max(abs(Xstate1(:,ndof1+1:Ndof)))';
    ZB2maxEQ2=max(abs(Xstate2(:,ndof1+1:Ndof)))';
    ZB2maxEQ3=max(abs(Xstate3(:,ndof1+1:Ndof)))';
    ZB1TOPmaxEQ1=ZB1maxEQ1(ndof1);
    ZB1TOPmaxEQ2=ZB1maxEQ2(ndof1);
    ZB1TOPmaxEQ3=ZB1maxEQ3(ndof1);
    ZB2TOPmaxEQ1=ZB2maxEQ1(ndof2);
    ZB2TOPmaxEQ2=ZB2maxEQ2(ndof2);
    ZB2TOPmaxEQ3=ZB2maxEQ3(ndof2);
    ZB1max=max([ZB1TOPmaxEQ1 ZB1TOPmaxEQ2 ZB1TOPmaxEQ3]);
    ZB2max=max([ZB2TOPmaxEQ1 ZB2TOPmaxEQ2 ZB2TOPmaxEQ3]);
    
    DB1max1=max((abs(Yatmd1(:,1:ndof1))))';
    DB2max1=max((abs(Yatmd1(:,ndof1+1:Ndof))))';
    DB1max2=max((abs(Yatmd2(:,1:ndof1))))';
    DB2max2=max((abs(Yatmd2(:,ndof1+1:Ndof))))';
    DB1max3=max((abs(Yatmd3(:,1:ndof1))))';
    DB2max3=max((abs(Yatmd3(:,ndof1+1:Ndof))))';
    DB1max=max([max(DB1max1) max(DB1max2) max(DB1max3)]);
    DB2max=max([max(DB2max1) max(DB2max2) max(DB2max3)]);

    Umax1=max(abs(Gatmd*Xstate1')');
    Umax2=max(abs(Gatmd*Xstate2')');
    Umax3=max(abs(Gatmd*Xstate3')');
    Umax=max([Umax1 Umax2 Umax3]);% maximum control force
    
    
    %fxj(jn3,1)=DB1max/D1limit+DB2max/D2limit+Umax/(5*Ulimit);
    %fxj(jn3,1)=DB1max/D1limit+DB2max/D2limit;
    fxj(jn3,1)=ZB1max/ZB1limit+ZB2max/ZB2limit;
    
    if Umax>Ulimit
       fxj(jn3,1)=100.0;
    end
    if fxj(jn3,1)>100.0
       fxj(jn3,1)=100.0;
    end
            
end
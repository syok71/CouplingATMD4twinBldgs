% 2 Building Model

clear
close all

%% Load Earthquakes
load elcen50EQ005S04short
t1=t;Xg1=Xg;
load kobeEQ005S04short
t2=t;Xg2=Xg;
load northREQ005S04short
t3=t;Xg3=Xg;

%% Structure Properties : Building 1
ndof1=10;
Height1=3.8;% height of each floor (meters)
fHeight1=Height1*[1:1:ndof1];
m1=123.2*10^3;% ton
M1=m1*eye(ndof1);iM1=M1^-1;
c01=ones(ndof1,1);
k01=1845.5*10^5*ones(ndof1,1);%N/m
K1=stiffmat(k01);
% Eigen Analysis and Construct Rayleigh-Damping Matrix
[w1,f1,Eivt1]=eigkm(K1,M1);
M1modal=(Eivt1(:,1)'*M1*Eivt1(:,1))/(Eivt1(ndof1,1)^2);
[ca01,ca11,xi01,C1]=rayDamp2(w1,1,2,0.01,0.01,M1,K1);

%% Structure Properties : Building 2
ndof2=10;
Height2=3.8;% height of each floor (meters)
fHeight2=Height2*[1:1:ndof2];
m2=110.7*10^3;% ton
M2=m2*eye(ndof2);iM2=M2^-1;
c02=ones(ndof2,1);
k02=2545.5*10^5*ones(ndof2,1);%N/m
K2=stiffmat(k02);
% Eigen Analysis and Construct Rayleigh-Damping Matrix
[w2,f2,Eivt2]=eigkm(K2,M2);
M2modal=(Eivt2(:,1)'*M2*Eivt2(:,1))/(Eivt2(ndof2,1)^2);
[ca02,ca12,xi02,C2]=rayDamp2(w2,1,2,0.01,0.01,M2,K2);

Ndof=ndof1+ndof2;
M0_12=addatmd0(ndof1,M1,ndof2,M2);
K0_12=addatmd0(ndof1,K1,ndof2,K2);
C0_12=addatmd0(ndof1,C1,ndof2,C2);

%% Packing and Svaing Object Models
OBJModels.t1=t1;
OBJModels.t2=t2;
OBJModels.t3=t3;
OBJModels.Xg1=Xg1;
OBJModels.Xg2=Xg2;
OBJModels.Xg3=Xg3;
OBJModels.ndof1=ndof1;
OBJModels.ndof2=ndof2;
OBJModels.K1=K1;
OBJModels.C1=C1;
OBJModels.K2=K2;
OBJModels.C2=C2;
OBJModels.Ndof=Ndof;
OBJModels.M0_12=M0_12;
OBJModels.C0_12=C0_12;
OBJModels.K0_12=K0_12;

save Models OBJModels

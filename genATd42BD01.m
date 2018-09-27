clear
close all;

warning off;

% Load Model Objectives
load Models;
Ndof=OBJModels.Ndof;
N=Ndof;
M0_12=OBJModels.M0_12;

% ATMD Mass matrix
n=1;% No. of ATMD
OBJModels.n=n;
Ntdof=Ndof+n;
Nst=2*Ntdof;
m_atmd=15000*2;
M1atmdM2=addatmd0(Ndof,M0_12,n,m_atmd);iM1atmdM2=M1atmdM2^-1;
OBJModels.iM1atmdM2=iM1atmdM2;

% ATMD Q matrix
Q=eye(Nst);
for j=1:n
    Jdx1=N+j;Jdx2=Ntdof+N+j;
    Q(Jdx1,Jdx1)=0;Q(Jdx2,Jdx2)=0;
end
OBJModels.Q=Q;

% Constraints
OBJModels.D1limit=0.05;% meter
OBJModels.D2limit=0.05;%meter
OBJModels.ZB1limit=0.28;% meter
OBJModels.ZB2limit=0.30;%meter
OBJModels.Ulimit=1350*1000;% Newton

save Models OBJModels

% Damping & Stiffness of ATMD
k_max=(2.5572e+06)*10;
k_min=5000.0;
c_max=(2.3267e+05)*10;
c_min=500.0;

% ATMD R
Ra_min=1;Ra_max=10;
Re_min=1;Re_max=16;

% Design Variables
Nvar=8;% Number of Design Variables
nALoc_bit=3;
nk_bit=40;
nc_bit=40;
nRa_bit=16;
nRe_bit=4;
BitArray=[nALoc_bit nk_bit nc_bit nk_bit nc_bit nRa_bit nRa_bit nRe_bit];
xmax=[10 k_max c_max k_max c_max Ra_max Ra_max Re_max];
xmin=[3 k_min c_min k_min c_min Ra_min Ra_min Re_min];

% Define Parameters of GA
Nch=50;% Number of Choromosomes for 1 generation
Pcr=0.85;% Crossover Probability (80%~95%)
MaxCpt=2;% Number of Maximum Crossover Point
Pmu=0.20;% Mutation Probability (generally 0.5%~5.0%)
Nsame=3000;% Total Number of Repeat Best Solution - Converge Requirement
Max_iter=100000;

G001=0;
%G001='10011110100100001100101110101000111101110011010101011111110101010101001100010110101100010';
[BestGene,BestG0,Ghistory,JiOpt,JiAve,Niter,Tlap]=GeneB1atmdB2(OBJModels,Nvar,xmax,xmin,BitArray,Nch,Pcr,Pmu,G001,Nsame,Max_iter);

% Plot Iteration Process
close all
figure
subplot(2,1,1), plot(JiOpt)
ylabel('Maximum fitness')
subplot(2,1,2), plot(JiAve)
xlabel('No. of iteration')      
ylabel('Average fitness')

save optTESTres3
save BestGeneCode BestGene Nvar BitArray xmax xmin

clear

SimBest



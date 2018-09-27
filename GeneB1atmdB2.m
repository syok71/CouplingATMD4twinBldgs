

function [BestGene,BestG0,Ghistory,JiOpt,JiAve,Niter,Tlap]=GeneB1atmdB2(OBJModels,Nvar,xmax,xmin,BitArray,Nch,Pcr,Pmu,G001,Nsame,Max_iter)

sameindex=0;
BestSol2=zeros(Max_iter,1);
fxj_mean=zeros(Max_iter,1);

%1.Generate Random Population of 'Nch' Choromosomes
Nbit=sum(BitArray);% number of total bit
G1(Nch,Nbit)='0';
Ghistory(Max_iter,Nbit)='0';
for ii=1:Nch
   for jj=1:Nbit
       if ii==1
           G1(ii,jj)='0';
       else
           G1(ii,jj)=num2str(rem(round(rand(1)*100),2));
       end
   end
end
if G001==0
    G001=G1(2,:);
    G1(2,:)=G001;% Initial guess is not used.
else
    G1(2,:)=G001;% Initial G0 is used.
end
New=G1;

tt0=clock;

for isol=1:Max_iter
      
   if sameindex > Nsame
      break
   end
  
   % 2. FITNESS Evaluation

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % User Define Fitness Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   fxj=fitFnB1atmdB2(OBJModels,Nch,G1,Nvar,BitArray,xmax,xmin);
   % User Define Fitness Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
           
   [BestSol2(isol,1),Bindex] = min(fxj);

   if isol==1
       BestG0=G1(Bindex,:);
       BestGene=BestG0;
   end
   
   if BestGene==G1(Bindex,:)
      sameindex=sameindex+1;
   else
      sameindex=0;
   end
   
   BestGene=G1(Bindex,:);
       
   Ghistory(isol,:)=BestGene;
   fxj_mean(isol,1)=mean(fxj);
   fxj_max=max(fxj);
   fit=fxj_max-fxj;
   S3=sum(fit);
      
   %3. Creat a New Population
	
	% 3.1 Selection
	% 3.1.1 Elitism
	[~,ii]=max(fit);
	New(1,:)=G1(ii,:);
	tem_fit=fit;
    tem_fit(ii)=0.0;
	[~,ii]=max(tem_fit);
	New(2,:)=G1(ii,:);
	
   	for iiii=2:Nch/2
                       
        % 3.1.2 Roulette Wheel Selection
		s3=0;rr3=rand(1)*S3;
		for ii=1:Nch
		   s3=s3+fit(ii);
		   if s3 > rr3
		      jj=ii;
		      break
		   end
		end
		P1=G1(jj,:);
  		s3=0;rr3=rand(1)*S3;
		for ii=1:Nch
		   s3=s3+fit(ii);
		   if s3 > rr3
		      jj=ii;      
		      break
		   end
        end
		P2=G1(jj,:);
		
		% 3.2 Crossover 
		rr3=rand(1);
        
        if rr3<Pcr
            r2=rem(round(rr3*100),2);
            Cpos=randperm(Nbit-1);
            temp1=P1;temp2=P2;
            if r2==0 % in case of 1 point Crossover
                Cposition=Cpos(1,1);
                P1=[temp1(1:Cposition) temp2(Cposition+1:Nbit)];
                P2=[temp2(1:Cposition) temp1(Cposition+1:Nbit)];
            else % in case of 2 points Crossover
                Cpos1=Cpos(1,1);Cpos2=Cpos(1,2);
                Cp=sort([Cpos1 Cpos2]);
                P1=[temp2(1:Cp(1)) temp1(Cp(1)+1:Cp(2)) temp2(Cp(2)+1:Nbit)];
                P2=[temp1(1:Cp(1)) temp2(Cp(1)+1:Cp(2)) temp1(Cp(2)+1:Nbit)];
            end
        end
	
		% 3.3 Mutation
        for ii=1:Nbit
            
		   r1=rand(1);
		   r2=rand(1);
           if r1<Pmu
               if P1(ii)=='0';
                   P1(ii)='1';
               else
                   P1(ii)='0';
               end 
           end
           if r2<Pmu
               if P2(ii)=='0';
                   P2(ii)='1';
               else
                   P2(ii)='0';
               end
           end
           
        end
        
        New(2*iiii-1,:)=P1;
      	New(2*iiii,:)=P2;
      
    end
    
    % 3.4 Update Newgeneration
    G1=New;
    if rem(isol,10)==0
        disp([num2str(isol) '-iteraton  ' BestGene])
    end
       
end

Tlap=etime(clock,tt0);

Niter=isol-1;
Ghistory=Ghistory(1:Niter,:);

JiOpt=BestSol2(1:Niter)';
JiAve=fxj_mean(1:Niter)';


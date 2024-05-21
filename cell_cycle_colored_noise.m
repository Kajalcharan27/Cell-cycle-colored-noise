% Paper: Subtle alteration in transcriptional memory governs the lineage-level cell cycle duration heterogeneities of mammalian cells
% Author: Kajal Charan and Sandip Kar
% e-mail about the code: kajalcharan97@gmail.com,sandipkar@iitb.ac.in
clc;
clear all;
close all;
aa=tic;% calculate execution time

nk=100; % number of lineages to simulate
 tou_G=10.*ones(1,4); % correlation time for lineages
 tou_S=10.*ones(1,4);
 seed_arr=[10 1001 2002]; % seed array 
 s_array=[0.001 0.005 0.01 0.015 0.02 0.025 0.03 0.035]; % noise strengths array 
 for see=1:length(seed_arr)
 for ns=1:length(s_array)
 ns1=s_array(ns); % noise strength of first TR
ku=4;
mycol={'r','b','m','g','m','k','b','r','g'}; % for plotting purpose
seed=seed_arr(see);
file_name=sprintf('t_a_u=%d s=%.3f seed=%d.dat',tou_G(1),ns1,seed); %generate file to store output
fid2=fopen(file_name,'w');
%tname=sprintf('%s,n=%d',baseFileName(1:end-4),N);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 d=2.8; % cell cycle period scalling parameter
 sf=1; % Scalling factor
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk=1:nk
    try
   zd=zeros(14,1);
    
%   tstart and tfinal
   % tstart = 0.0;
    tfinal = 150.0;
     sf=1;
%	initial values of mrna and proteins
	
tcv=0.05;
sig=1.0;
  ark(1)=0.0037*d; % mean transcription rates 
    ark(2)=0.005*d;
    ark(3)=0.2*d;
    ark(4)=0.5*d;  
    ark;
    for i=1:ku; % generating log normal distribution for TR 
        rsigma(i)=sqrt(log(1.0+(tcv*tcv)));
        rmu(i)=log(ark(i))-(rsigma(i)*rsigma(i))/2;
    end % mu and sigma for initial transcription rate
rng(seed);
for j=1:ku
    
        rk1=rand();
        rk2=rand();
        api=4.0*atan(1.0);
        if (rk1<=0.0);
            rk1=0.00001;
        end
        ba=-2.0*sig*sig*log(rk1);
        ba1=cos(2.0*api*rk2);
        sff=sqrt(ba)*ba1;
        c_arr(1,j)=(exp(rmu(j)+rsigma(j)*sff));
end % initial transcription rate from lognormal distribution
%	
dt=0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[dk1m,param(1)]=deal((c_arr(kk,1))*sf);%%% 
[dk1m,param(1)]=deal(c_arr(1,1)*sf);
param(1);
    [dk1dm,param(2)]=deal(0.058*d);
    [dk1,param(3)]=deal(0.4*d);
    [dk2a,param(4)]=deal(0.04*d);
    [dk2b,param(5)]=deal(2.0*d/sf);

%	cdh1	
	[dk3a,param(6)]=deal(1.0*d);
    [dk3b,param(7)]=deal(8.0*d);%----
    [dk4,param(8)]=deal(40.0*d);
    [dj3,param(9)]=deal(0.04*sf);
    [dj4,param(10)]=deal(0.04*sf);
    [dk3m,param(11)]=deal((c_arr(1,4))*sf);%%%
    [dk3dm,param(12)]=deal(0.5*d);
    [dk3dt,param(13)]=deal(1.0*d);%----
%	#cdc20m
    [dk5am,param(14)]=deal((c_arr(1,2))*sf);%%%%
    [dk5bm,param(15)]=deal((c_arr(1,3))*sf);%%%%%
    [dk5dm,param(16)]=deal(1.386*d);
    [dj5,param(17)]=deal(0.3*sf);
    [dn,param(18)]=deal(4.0);
    
%	#cdc20p
    [dk5a,param(19)]=deal(1.0*d);
    [dk6,param(20)]=deal(0.05*d);
    
%	#cdc20a
    [dk7,param(21)]=deal(1.4*d);%----
    [dk8,param(22)]=deal(0.5*d*sf);
    [dj7,param(23)]=deal(0.001*sf);
    [dj8,param(24)]=deal(0.001*sf);
    [dmad,param(25)]=deal(1);
%	#IE
	[dk9,param(26)]=deal(0.1*d/sf);
    [dk10,param(27)]=deal(0.02*d);
    [dk11,param(28)]=deal(0.045*d*sf);
    [dk12,param(29)]=deal(2.27*d/sf);
    [dk13,param(30)]=deal(0.004*d);
    [dk3,param(31)]=deal(1.28*d*sf);%--
    [dgf,param(32)]=deal(2.0); 
    [dkmm,param(33)]=deal(0.2);%-- 
    [dkeff,param(34)]=deal(1.0);
    param(35)=sf;
    [dk5cm,param(36)]=deal(1.0);
    [dj5c,param(37)]=deal(0.02);%----
    %[dtou,param(38)]=deal(10.0);
    [dtou,param(38)]=deal(tou_G(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial condtions
%set to some very large number

idiv=1;
    t_arr(idiv)=0.0;
	[inbm,zd(1,1)]=deal(0.05*sf);
	[incycb,zd(2,1)]=deal(0.3*sf); 
	[inc20m,zd(3,1)]=deal(0.07*sf);
	[incdc20t,zd(4,1)]=deal(0.59*sf);
	[incdhm,zd(5,1)]=deal(0.99*sf);
    [incdht,zd(6,1)]=deal(0.99*sf);
    [incdh1,zd(7,1)]=deal(0.006*sf);
	[incdc20a,zd(8,1)]=deal(0.036*sf);
	[iniep,zd(9,1)]=deal(0.41*sf);
    [incdt1,zd(10,1)]=deal(0.06*sf);     
    [inc_noise1(idiv),zd(11,1)]=deal(0.0);
    [inc_noise2(idiv),zd(12,1)]=deal(0.0);
    [inc_noise3(idiv),zd(13,1)]=deal(0.0);
    [inc_noise4(idiv),zd(14,1)]=deal(0.0);
z_arr(:,idiv)=zd(:,1);
seed;
rng(seed) ; %start from the same seed each time
seed=seed+10;
%create global variables for use later
    
    tevent   = [];
    zevent   = [];
    ievent   = [];
% using events -- very helpful in computing precisely when cell division occurs.
%options = odeset('RelTol',1e-5,'AbsTol',1e-5);

options = odeset('Events',@colored_new_events_cellcycle,'RelTol',1e-5,'AbsTol',1e-5);

 
 
gen=2; % generations 
    l1=0;
    cycle_time=[];
   ij=0;
% while tstart<tfinal
     for ii=0:gen
      
         for jj=1:2^(ii)
             param(38)=tou_G(ii+1);
            
             time     = [];
             zsol     = [];
             
              if (mod(jj,2))~=0
                    ij=ij+1;
              end
              ij;
              z0=z_arr(:,ij);
              tstart=t_arr(ij);
              
              paramk(1)=param(1);
              param(1)=paramk(ij);
              param2(1)=param(11);
              param(11)=param2(ij);
              param3(1)=param(14);
              param(14)=param3(ij);
              param4(1)=param(15);
              param(15)=param4(ij);
              
              tspan  = [tstart:dt:tfinal];
              step=size(tspan);
              % white noise generation from random numbers
              noise_1=ns1*randn(step);
              noise_2=ns1*randn(step);%k3m*
              noise_3=ns1*randn(step);% k5am*
              noise_4=ns1*randn(step);%k5bm*
              wt.tspan=tspan;
              wt.noise1=noise_1;
              wt.noise2=noise_2;
              wt.noise3=noise_3;
              wt.noise4=noise_4;

    ie=0;
    
    while (ie~=1) & (tstart < tfinal)
        param(1);
             [t,z,te,ze,ie]=ode15s(@(t,z)colored_new_cellcycle_rhs(t,z,wt,param),[tstart:dt:tfinal],z0,options);


     %[t,z]=ode15s(@cellcycle_rhs,[tstart tfinal],z0,options,param);
             %set tstart to the last time calculated
             
            %fill up the global variables
              time     = [time;t];
              zsol     = [zsol;z]; 
              tevent   = [tevent;te];
              zevent   = [zevent;ze];
              ievent   = [ievent;ie];
              
             z0=zsol(end,:);
              tstart=t(end); 
             
              
              %end
              ie;t(end);
              if (isempty(ie))
                  ie=4;
              elseif (ie(end)==1)
                 clear max;
                    % to calculate G1 duration from the maxima of variable
                    % cdt1
                     div_time=t(end)-t_arr(ij);
                     max_cdt1=max(zsol(:,10));
                     g1=time(find(zsol(:,10)==max_cdt1));
                     g1_time=g1-t_arr(ij);
                     sg2m_time=div_time-g1_time;
                     [div_time g1_time sg2m_time];
                  
                  %param(38)=tou_G(ii+1);% G1 correlation time
              elseif (ie(end)==2)
                  cv=2;
                  % TR variation before Mitosis
                  param(1)
                  param(1)=trvariation(param(1),cv);
                  param(1)
                  param(11);
                  param(11)=trvariation(param(11),cv);
                  param(11);
                  param(14)=trvariation(param(14),cv);
                  param(15)=trvariation(param(15),cv);
                  
              elseif (ie(end)==3)
                  param(38)=tou_S(ii);%SG2M correlation time
              end
    end
    
    if (ie(end)==1)
        idiv=idiv+1;
        paramk(idiv)=param(1);
        param2(idiv)=param(11);
        param3(idiv)=param(14);
        param4(idiv)=param(15);
                  f=0.5;
              % equal partitioning of cell cycle components
                  z_arr(1:end-5,idiv)=(zsol(end,1:end-5))*f;
                  z_arr(end-4,idiv)=zsol(end,end-4)*f;
                 
                  z_arr(end-3:end,idiv)=zsol(end,end-3:end);
                 
                  t_arr(idiv)=t(end);
                 
                  %plot(time,(zsol(:,11)+param(1)),'color',mycol{ij},'LineStyle','-', 'Marker', 'none','LineWidth',0.5);
                   %hold on;
    else
        
        idiv=idiv+1;
        t_arr(idiv)=t(end);
        z_arr=zsol(end,:);
        
       %plot(time,(zsol(:,11)+param(1)),'color',mycol{ij},'LineStyle','-', 'Marker', 'none','LineWidth',0.5); 
    end
    %{
    if div_time>40
        flag=1;
        break
        
    end
   %}
        
             
    %% Plot and store
    % 
    
    %yyaxis left;
           %plot(time,(zsol(:,11)),'color',mycol{ij},'LineStyle','-', 'Marker', 'none','LineWidth',0.5);
              
             %hold on;
             
             %xlswrite(filename,[time zsol(:,11)])
            % yyaxis right; 
            %plot(time,zsol(:,8),'color',mycol{2},'LineStyle','-', 'Marker', 'none')
             %ylim([-0.0001 0.0001])
             %set(gca,'FontSize',20,'FontName','Times New Roman','FontWeight','b');
             
             
             if t(end)<72
                % plot(time,(zsol(:,2)),'color',mycol{ij},'LineStyle','-', 'Marker', 'none','LineWidth',0.5);
              
                 if ii>0
                     
                     % store the output in file
                    fprintf(fid2,'%.2f %.2f %.2f %5f %3f %4f\r',div_time, g1_time, sg2m_time,kk,ii,jj);
                 end
             end
             %{
                 if ii==1 % store data for first generation 
               if jj==1 
                   mat_mrna=[time zsol(:,1)];
                   gen1_mat_mrna(1:size(mat_mrna,1),4*kk-3:4*kk-2)=mat_mrna;
                   mat=[time zsol(:,11)];
                   gen1_mat(1:size(mat,1),4*kk-3:4*kk-2)=mat;
               elseif jj==2
                   jj
                   mat_mrna=[time zsol(:,1)];
                   gen1_mat_mrna(1:size(mat_mrna,1),4*kk-1:4*kk)=mat_mrna;
                   mat=[time zsol(:,11)];
                   gen1_mat(1:size(mat,1),4*kk-1:4*kk)=mat;
               end
               % store data for second generation (time,c_noise1)
               % 8 column --> 1 Lineage 
            elseif ii==2
                if jj==1 
                   mat_mrna=[time zsol(:,1)];
                   gen2_mat_mrna(1:size(mat_mrna,1),8*kk-7:8*kk-6)=mat_mrna;
                   mat=[time zsol(:,11)];
                   gen2_mat(1:size(mat,1),8*kk-7:8*kk-6)=mat;
               elseif jj==2
                   jj
                   mat_mrna=[time zsol(:,1)];
                   gen2_mat_mrna(1:size(mat_mrna,1),8*kk-5:8*kk-4)=mat_mrna;
                   mat=[time zsol(:,11)];
                   gen2_mat(1:size(mat,1),8*kk-5:8*kk-4)=mat;
               elseif jj==3
                   jj
                   mat_mrna=[time zsol(:,1)];
                   gen2_mat_mrna(1:size(mat_mrna,1),8*kk-3:8*kk-2)=mat_mrna;
                   mat=[time zsol(:,11)];
                   gen2_mat(1:size(mat,1),8*kk-3:8*kk-2)=mat;
               elseif jj==4
                   jj
                   mat_mrna=[time zsol(:,1)];
                   gen2_mat_mrna(1:size(mat_mrna,1),8*kk-1:8*kk)=mat_mrna;
                   mat=[time zsol(:,11)];
                   gen2_mat(1:size(mat,1),8*kk-1:8*kk)=mat;
               end
            end
            %}
           [div_time g1_time sg2m_time kk ii jj]
           
         end
         %{
         if flag==1
             kk
             break
         end
          %}  
        
     end
     hold off
    catch
    end
     end
    %{
    writematrix(gen1_mat,'seed=10_10h_0.001_gen1_file.csv')
    writematrix(gen2_mat,'seed=10_10h_0.001_gen2_file.csv')
    writematrix(gen1_mat_mrna,'seed=10_10h_0.001_gen1_file_mRNA.csv')
    writematrix(gen2_mat_mrna,'seed=10_10h_0.001_gen2_file_mRNA.csv')
    %}

     
       fclose(fid2);
 %end    
    box on;
    
 end 
 end
 whole_time=toc(aa)   
   

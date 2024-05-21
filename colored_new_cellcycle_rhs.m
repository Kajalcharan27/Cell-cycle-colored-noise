% Paper: Subtle alteration in transcriptional memory governs the lineage-level cell cycle duration heterogeneities of mammalian cells
% Author: Kajal Charan and Sandip Kar
% e-mail about the code: kajalcharan97@gmail.com,sandipkar@iitb.ac.in
function zp=colored_new_cellcycle_rhs(t,z,wt,param)
%
% all the variables in the cell cycle problem
%	initial values of mrna and proteins
	
%clear c_noise;
	
	bm=z(1);
	cycb=z(2); 
    c20m=z(3);
	cdc20t=z(4);
    cdhm=z(5);
    cdht=z(6);
	cdh1=z(7);
	 
	cdc20a=z(8);
	iep=z(9);
    cdt1=z(10);
	c_noise1=z(11);
	c_noise2=z(12);
	c_noise3=z(13);
	c_noise4=z(14);
	

% all the parameters in the cell cycle problem
	%parameter values

	%cycb
	k1m=param(1);
    k1dm=param(2);
    k1=param(3);
	k2a=param(4);
	k2b=param(5);
	%cdh1	
	k3a=param(6);
	k3b=param(7);	
	k4=param(8);
	j3=param(9);
	j4=param(10);
    k3m=param(11);
    k3dm=param(12);
    k3dt=param(13);
	%cdc20m
	k5am=param(14);
	k5bm=param(15);
	k5dm=param(16);
	j5=param(17);
	n=param(18);
    %cdh1
	k5a=param(19);
	k6=param(20);
    %cdc20a
	k7=param(21);
	k8=param(22);
	j7=param(23);
	j8=param(24);
	mad=param(25);
	%ie
    k9=param(26);
	k10=param(27);
    k11=param(28);
    k12=param(29);
    k13=param(30);
    k3=param(31);
    gf=param(32);
    kmm=param(33);
    keff=param(34);
    sf=param(35);
    k5cm=param(36);
    j5c=param(37);
    tou=param(38);
   
    tou_1=tou;
    tou_2=tou;
    tou_3=tou;
    tou_4=tou;
	
    %if Cdt1==2000
     %  param(1)=param(1)+0.1; 
    %end

 
% all the differential equations in the population space
%Volume
%dV/dt=kve*V*GF/(Vmv+GF)
%dV/dt=kv*GF/(Vmv+GF)

%s=10000
%Dp1p=((Dpt)-DE-iDEP);
%E2Fp=(Y-DE-iDE-iDEP);
%Rb=RbT-iRb1-iRb2-iDE-iDEP;
%CycBi=CycBT-CycBa;
%C25i=C25T-C25a;
%Wee1i=Wee1T-Wee1a;w_value1=wt.noise1(find(wt.tspan>=t,1,'first'));
w_value1=wt.noise1(find(wt.tspan>=t,1,'first'));
w_value2=wt.noise2(find(wt.tspan>=t,1,'first'));
w_value3=wt.noise3(find(wt.tspan>=t,1,'first'));
w_value4=wt.noise4(find(wt.tspan>=t,1,'first'));
%c_noise=(y*w_value)-(y*c_noise);


zp(1)=(k1m+c_noise1)*gf/(kmm+(keff*gf))-k1dm*bm;
zp(2)=k1*bm-k2a*cycb-k2b*cycb*cdh1;
zp(3)=(k5am+c_noise2)+(k5bm+c_noise3)*((cycb/j5)^(n))/((k5cm+(gf*j5c))*(1+((cycb/j5)^(n))))-k5dm*c20m;
zp(4)=k5a*c20m-k6*cdc20t;
zp(5)=(k3m+c_noise4)-k3dm*cdhm;
zp(6)=k3a*cdhm-k3dt*cdht;
zp(7)=(k3+k3b*cdc20a)*(cdht-cdh1)/(j3+cdht-cdh1)-k4*cycb*cdh1/(j4+cdh1)-k3dt*cdh1;
zp(8)=k7*iep*(cdc20t-cdc20a)/(j7+cdc20t-cdc20a)-k8*mad*cdc20a/(j8+cdc20a)-k6*cdc20a;

zp(9)=k9*cycb*(sf-iep)-k10*iep;

zp(10)=k11-k12*cycb*cdt1-k13*cdt1;
zp(11)=(w_value1/tou_1)-(c_noise1/tou_1);
zp(12)=(w_value2/tou_2)-(c_noise2/tou_2);
zp(13)=(w_value3/tou_3)-(c_noise3/tou_3);
zp(14)=(w_value4/tou_4)-(c_noise4/tou_4);
%r=0+(1-0)*rand(1);

zp=zp(:);


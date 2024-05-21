% Paper: Subtle alteration in transcriptional memory governs the lineage-level cell cycle duration heterogeneities of mammalian cells
% Author: Kajal Charan and Sandip Kar
% e-mail about the code: kajalcharan97@gmail.com,sandipkar@iitb.ac.in

function [akn] = trvariation(akn,vf)
%TRVARIATION Summary of this function goes here
%   Detailed explanation goes here
%rng(seed);
	%akn=akm;
    ku=1;
    
j=1;

 		rk3 = rand();

        if (rk3 <= 0.5)
            if (rk3>=0.4)
			akn(j)=akn(j)-(0.01*vf*akn(j));
			
            elseif (rk3>=0.3)
			akn(j)=akn(j)-(0.02*vf*akn(j));
			

            elseif (rk3>=0.2)
			akn(j)=akn(j)-(0.03*vf*akn(j));

            elseif (rk3>=0.1)
			akn(j)=akn(j)-(0.04*vf*akn(j));
			
            elseif (rk3>=0.0)
			akn(j)=akn(j)-(0.05*vf*akn(j));
			
            end
            
        elseif (rk3 >= 0.5)

            if (rk3<=0.6)
			akn(j)=akn(j)+(0.01*vf*akn(j));
			
            elseif (rk3<=0.7)
			akn(j)=akn(j)+(0.02*vf*akn(j));
		

            elseif (rk3<=0.8)
			akn(j)=akn(j)+(0.03*vf*akn(j));
		
            elseif (rk3<=0.9)
			akn(j)=akn(j)+(0.04*vf*akn(j));
		
            elseif (rk3<=1.0)
			akn(j)=akn(j)+(0.05*vf*akn(j));
		
            end 

        
end




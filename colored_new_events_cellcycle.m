% Paper: Subtle alteration in transcriptional memory governs the lineage-level cell cycle duration heterogeneities of mammalian cells
% Author: Kajal Charan and Sandip Kar
% e-mail about the code: kajalcharan97@gmail.com,sandipkar@iitb.ac.in
function [value,isterminal,direction] = colored_new_events_cellcycle(t,z,param)
    
    
    %	initialization for period calculation
	maxc=0;
	maxg=0;
	tc=0;
	ncell=1;
	eg(ncell)=maxc;
	teg(ncell)=tc;
	mcount=1;
	ncount=0;
	  isc=10;
      
      
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
    %{
    if cdh1<=0.2
        
        value=[1];
        z(7);
        isterminal=[0];
        direction=[-1]; 
    else
        cdh1
        value=[cycb-0.1];
        z(7);
        isterminal=[0];
        direction=[-1]; 
    end
   
   
    celldiv=0;
    
  mcount=1;
  ncount=1;
    if cycb<=0.1 & cdh1>=0.5 & mcount==1
       
        celldiv=1;
        mcount=0;
        ncount=1;
    end         
    if cycb>=0.1 &cdh1>=0.5 & ncount==1
       mcount=1;
       ncount=0;
       celldiv=0;
    end
                        
    value=[celldiv-1];
    isterminal=[1];
    direction=[1];
            
%} 
value=[cycb-0.1;cycb-0.25;cdt1-0.6];
    isterminal=[1;1;1];
    direction=[-1;1;-1];




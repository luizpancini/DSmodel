function [adaptiveInertia,adaptiveSwarmSize,alwaysUpdateInf,c0,c1,c2,D,EM_swarms,initializePos,initializeVel,improveTol,maxFE,maxStalledIt,max_c0,min_c0,moveMethod,N,N_c,N_i,N_m,N_e,num_workers_parfor,reboundVel,run_parfor,set_GP,specified_filename,stopOnMaxStalledIt,setNMLOnMaxStalledIt,topology] = unpack_BLPSO_options(options)

adaptiveInertia = options.adaptiveInertia;           
adaptiveSwarmSize = options.adaptiveSwarmSize;          
alwaysUpdateInf = options.alwaysUpdateInf;            
c0 = options.c0;                       
c1 = options.c1;                        
c2 = options.c2;  
D = options.D;
EM_swarms = options.EM_swarms;                  
initializePos = options.initializePos;   
initializeVel = options.initializeVel;
improveTol = options.improveTol;
maxFE = options.maxFE;  
maxStalledIt = options.maxStalledIt;
max_c0 = options.max_c0;                  
min_c0 = options.min_c0;                 
moveMethod = options.moveMethod; 
N = options.N;  
N_c = options.N_c;       
N_i = options.N_i; 
N_m = options.N_m;       
N_e = options.N_e; 
num_workers_parfor = options.num_workers_parfor;  
reboundVel = options.reboundVel;   
run_parfor = options.run_parfor;
set_GP = options.set_GP;
specified_filename = options.specified_filename;
stopOnMaxStalledIt = options.stopOnMaxStalledIt;         
setNMLOnMaxStalledIt = options.setNMLOnMaxStalledIt;       
topology = options.topology;
        
end
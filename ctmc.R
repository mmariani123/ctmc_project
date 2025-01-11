#!/usr/bin/env R

#CTMC project R code

#Below we run the MATLAB/OCTAVE scripts from OLIVIA:

############## MATLAB OCTAVE CODE: ######################################

#1.) Create_Parameter_Set.m
#2.) nParasiteGroups_ContRuptFunc.m

#Remember may need to load the staistics package in R
#pkg install -forge statistics ?
#% pkg load statistics ?

#work_dir='G:/My Drive/mariani_systems/malaria_ctmc_project/MS1ReviewCode/MS1ReviewCode'
#addpath(work_dir);
#par = Create_Parameter_Set;
#%out = CTMC_ContRuptFunc(par);
#nParasiteGroups_ContRuptFunc;

######## HEre we create the parameters:

#% Function creates structure of parameter values.  Number of genotypes,
#% simulations, and max fitness bias can also be modified here.
#function param = Create_Parameter_Set

#The adaptivetau package in R implements both an exact solution 
#and an approximate solution known as the “adaptive tau-leaping algorithm” 
#[Cao et al., 2007]. Similar functionality exists in the publicly-available 
#GillespieSSA R package [Pineda-Krch, 2008]; however, our new implementation 
#is much faster, due in part to its underlying implementation in C++.

library('adaptivetau')  #Philip Johnson, 2024
library('futureverse')  #Henrik Bengtsson, 2024, https://www.futureverse.org/
library('pracma')

params = list(
N        = 2,           #% Number of different starting genotypes
a        = 24/(20/60),  #% failure rate of male gametocytes (per day)
b        = 24/(25/60),  #% failure rate of female gametocytes (per day)
r        = 0.08,        #% fertilization of male and female gametes (per day)
mu_z     = 1,           #% death rate of zygotes (per day)
sigma_z  = 24/19,       #% transformation rate of zygotes (per day)
sigma_e  = 0.6,         #% transformation rate of ookinetes (per day)
mu_e     = 1.4,         #% death rate of ookinetes (per day)
mu_o     = 0,           #% death rate of oocysts (per day) **changed from paper**
n        = 3e3,         #% number of sporozoites per oocyst
p        = 0.2,         #% proportion of sporozoites that make it to salivary gland
alpha    = 0.39,        #% fraction of male gametes that are viable
beta     = 1,           #%0.96;  #% fraction of female gametes that are viable
eta      = 1,           #% number of female gametes per female gametocyte
rho      = 8,           #% number of male gametes per male gametocytes
k        = 1/7,
t0       = 10,
max_bias = 0.1,         #%0.5; #%0.1;
NumSim   = 2            #%10000;
)

##############################################################################

#Here we recode the nParasiteGroups_ContRuptFunc;

#% -- Simulate CTMC model of within-vector parasite dynamics using modified
#% Gillespie Algorithm (outlined in Supplementary Materials), with
#% continuous rupture function.

#close all
#clear

#SEED = rng('shuffle'); # set seed
set.seed(42)

#% Set parameters equal to values defined in Create_Parameter_Set.m:
#par = Create_Parameter_Set;

N = params[['N']]
a = params[['a']]
b = params[['b']]
r = repmat(params[['r']],params[['N']],params[['N']])
#r = matrix(params[['r']],params[['N']],params[['N']])
mu_z = params[['mu_z']]
sigma_z = params[['sigma_z']]
sigma_e = params[['sigma_e']]
mu_e = params[['mu_e']]
mu_o = params[['mu_o']]
n = params[['n']]
p = params[['p']]
alpha = params[['alpha']]
beta = params[['beta']]
eta = params[['eta']]
rho = params[['rho']]
k = params[['k']]
max_bias = params[['max_bias']]
t0 = params[['t0']]

#% proportion of gametocytes that are male of each genotype
percent_male = repmat(0.25,N,1)

#% Define vector of fitness bias for each genotype: 
#(this can be modified to consider >3 genotypes)
#% For N>3, must define vector of fitness biases.
if(N == 1){
  bias = 0
}else if(N == 2){
  bias = c(0,max_bias)
}else if(N == 3){
  bias = c(0, 0.1, 0.5)
}

#bias = bias';
a_vec     = repmat(a,N,1)*(1 + bias);
b_vec     = repmat(b,N,1)*(1 + bias);
alpha_vec = repmat(alpha,N,1)*(1 - bias);
beta_vec  = repmat(beta,N,1)*(1 - bias);
mu_z_vec  = repmat(mu_z,N,1)*(1 + bias); # zygote mortality rate for zygotes whose parents are of the same genotype
mu_e_vec  = repmat(mu_e,N,1)*(1 + bias); # ookinete mortality rate for ookinetes " ...
mu_o_vec  = repmat(mu_o,N,1)*(1 + bias); # oocyst mortality rate oocysts " ...

# Create symmetric matrices of biased parameters.  The (i,j) entry is the
# parameter value for a parasite that has a genotype (i) female, and
# genotype (j) male parent.

#mu_z_mat = NaN(N,N)
#mu_e_mat = NaN(N,N)
#mu_o_mat = NaN(N,N)

mu_z_mat = repmat(NaN,2)
mu_e_mat = repmat(NaN,2)
mu_o_mat = repmat(NaN,2)

for(i in 1:N)
    for(j in 1:N)
      mu_z_mat[i,j] = mean(mu_z_vec[i],mu_z_vec[j]);
      mu_e_mat[i,j] = mean(mu_e_vec[i],mu_e_vec[j]);
      mu_o_mat[i,j] = mean(mu_o_vec[i],mu_o_vec[j]);
    end
end

# ---- Run NumSim number of simulations for each value of G0 in G0_vec as
# defined in CreateParameterSet.m:

# Number of simulations %
NumSim = params[['NumSim']]

# Max number of iterations for each simulation %
maxiter = 1e3

# Vector of Initial gametocyte densities %
G0_vec = seq(from=150,by=50,to=450)

# Max time (in days)
Tfinal = 21;

# Time-step to increment by if the time to next time-step is greater than
#% dt.
dt = 0.1

# Fraction of G0 of each subtype (Nx1 vector):
strainprop = repmat(1/N,N,1)

#TransitionIDs = NaN(NumSim,maxiter)
TransitionIDs = repmat(NaN,NumSim,maxiter)

for(ii in 1:length(G0_vec))

    G0 = round(strainprop*G0_vec[ii]) 
    # Vector of initial gametocyte densities.  
    G0[i] #is the # of genotype i gametocytes.

    Gf0 = ceil((1 - percent_male)*G0*eta*beta) 
    
    # Male gamete densities of each genotype
    Gm0 = ceil(percent_male*G0*rho*alpha)

    # Initialize Time_data vector.
    #Time_data = NaN(maxiter,NumSim) 
    Time_data = repmat(NaN,maxiter,NumSim)
    
    #Each row is a different genotype; 
    #each column is a different time step; 
    #each layer is a different simulation
    
    #Male_data       = zeros(N,maxiter,NumSim) 
    #Female_data     = zeros(N,maxiter,NumSim)   #% Same as for Male_data
    #is the # of Zygotes with female parent i, male parent j, at timestep k, in simulation l
    #Zygote_data     = zeros(N,N,maxiter,NumSim) #% 4D array: Zygote(i,j,k,l)
    #Ookinete_data   = zeros(N,N,maxiter,NumSim) #% "
    #Oocyst_data     = zeros(N,N,maxiter,NumSim) #% "
    #Sporozoite_data = zeros(N,N,maxiter,NumSim) #% "
    #Burst_time      = zeros(N,N,maxiter,NumSim) #% "

    Male_data       = array(0, c(N, N, maxiter, NumSim))  
    Female_data     = array(0, c(N, N, maxiter, NumSim))
    #is the # of ...
    Zygote_data     = array(0, c(N, N, maxiter, NumSim))
    Ookinete_data   = array(0, c(N, N, maxiter, NumSim))
    Oocyst_data     = array(0, c(N, N, maxiter, NumSim))
    Sporozoite_data = array(0, c(N, N, maxiter, NumSim)) 
    Burst_time      = array(0, c(N, N, maxiter, NumSim))
    
    #tic
    start.time <- Sys.time()
    
    for (j in 1:NumSim)

      #Male_data(:,1,j) = Gm0; % Set initial conditions for male and female gamete data.
      #Female_data(:,1,j) = Gf0; % "
    
      # Set initial conditions for male and female gamete data.
      Male_data[1:N,1,j,1:NumSim]   = Gm0 
      Female_data[1:N,1,j,1:NumSim] = Gf0 #% "

      # Initialize burst_count 3-D array.
      #burst_count(:,:,j) = zeros(N,N) 
      burst_count = array(0, c(N,N,NumSim))
      
      t = 0
      t[1] = 0
      time = t[1]

      for(i in 2:maxiter)

        #y1 = rand;
        #y2 = rand;

        y1 = runif(1)
        y2 = runif(2)
        
        #State variables
        #Column vector where each row is a different genotype, 
        #and the elements are the # of Male gametes at the previous time-step (t_{i-1}) 
        #in simulation j.  So, m(k) is the number of genotype-k males at time t_{i-1} 
        #in sim j.
        m = Male_data[1:N,i-1,j,1:NumSim]
        f = Female_data[1:N,i-1,j,1:NumSim]
        # z(k,l) = # zygotes with female-k parent and male-l parent at time t_{i-1} 
        # in sim j.
        z = Zygote_data(:,:,i-1,j)
        e = Ookinete_data(:,:,i-1,j)
        o = Oocyst_data(:,:,i-1,j);
        s = Sporozoite_data(:,:,i-1,j)

        # Break if states m through o are extinct:
        state_variables = c(m,
                            f,
                            reshape(z,N*N,1),
                            reshape(e,N*N,1),
                            reshape(o,N*N,1))
        
        if(sum(state_variables(:)) == 0){
          LastIter(j) = i-1;
          break
        }

        # Assign data values at current time-step, i, 
        # to values at previous time-step, i-1.
        Male_data(:,i,j)         = m
        Female_data(:,i,j)       = f
        Zygote_data(:,:,i,j)     = z
        Ookinete_data(:,:,i,j)   = e
        Oocyst_data(:,:,i,j)     = o
        Sporozoite_data(:,:,i,j) = s

        #Define bursting function.
        # probability of bursting at time 'time'
        bursting = 1/(1 + (exp(t0 - t(i-1)))) 
        #Create vectors/matrices of all possible transitions, 
        #arranged by type: death, mating, maturation, bursting.
        
        Male_death = a_vec.*m
        Female_death = b_vec.*f
        Zygote_death = mu_z_mat.*z
        Ookinete_death = mu_e_mat.*e
        Oocyst_death = mu_o_mat.*o
        Mating = r*(f*m)
        Zyg_maturation = sigma_z*z
        Ook_maturation = sigma_e*e
        Ooc_Bursting = k.*bursting.*o

        #Calculate all possible transitions:
        transitions = c(Male_death,
                        Female_death,
                        reshape(Zygote_death,N*N,1),
                        reshape(Ookinete_death,N*N,1),
                        reshape(Oocyst_death,N*N,1),
                        reshape(Mating,N*N,1),
                        reshape(Zyg_maturation,N*N,1),
                        reshape(Ook_maturation,N*N,1),
                        reshape(Ooc_Bursting,N*N,1)
                       )
  
  
        ## Create vector of indices marking different transition types in
        ## the transition vector:
        transition_type_indices = c(0,
                                    cumsum(
                                      c(length(Male_death),
                                        length(Female_death),
                                        length(reshape(Zygote_death,N*N,1)),
                                        length(reshape(Ookinete_death,N*N,1)),
                                        length(reshape(Oocyst_death,N*N,1)),
                                        length(reshape(Mating,N*N,1)),
                                        length(reshape(Zyg_maturation,N*N,1)),
                                        length(reshape(Ook_maturation,N*N,1)),
                                        length(reshape(Ooc_Bursting,N*N,1))
                                       )
                                     )
                                    )
    
        # Calculate vector of transition probabilities %
        total_transitions = sum(transitions)
        trans_probs = transitions/total_transitions
      
        # Calculate time of next event.  If time > Tfinal, break!
        stoch_time_step = -log(y1)/(total_transitions)
        time1 = stoch_time_step +t(i-1) # time until next event;
        time2 = t(i-1) + dt
        
        if(min(time1,time2) > Tfinal)
        {
          LastIter(j) = i-1
          break
        }elseif(stoch_time_step > dt){
          t(i) = time2
        }elseif(stoch_time_step <= dt){
          t(i) = time1
        }
        
        # ---------------------- %
        # Determine next event and update population matrices %
        
        choose_transition = find(cumsum(trans_probs)>y2,1,'first')
        
        
        TransitionIDs(j,i) = choose_transition;
        
        if(choose_transition <= transition_type_indices(2)){
        
          #% Male gamete death
          CT = choose_transition - transition_type_indices(1)
          Male_data(CT,i,j) = m(CT) - 1
        
        }else if((choose_transition > transition_type_indices(2)) && (choose_transition <= transition_type_indices(3))){
        
          #% Female gamete death 2
          CT = choose_transition - transition_type_indices(2)
          Female_data(CT,i,j) = f(CT) - 1
        
        }else if((choose_transition > transition_type_indices(3)) && (choose_transition <= transition_type_indices(4))){
        
          #% Zygote Death 3
          #CT_mat = reshape(transition_type_indices(3)+1:transition_type_indices(4),N,N)
          #C(ix1,ix2) = find(CT_mat == choose_transition)
        
          Zygote_data(ix1,ix2,i,j) = z(ix1,ix2) - 1
        
        }else if((choose_transition > transition_type_indices(4)) && (choose_transition <= transition_type_indices(5))){
        
          #% Ookinete Death 4
          #CT_mat = reshape(transition_type_indices(4)+1:transition_type_indices(5),N,N)
          #[C(x1,ix2) = find(CT_mat == choose_transition)
          Ookinete_data(ix1,ix2,i,j) = e(ix1,ix2) - 1
        
        }else if((choose_transition > transition_type_indices(5)) && (choose_transition <= transition_type_indices(6))){
        
          #% Oocyst Death 5
          #CT_mat = reshape(transition_type_indices(5)+1:transition_type_indices(6),N,N);
          C(ix1,ix2) = find(CT_mat == choose_transition)
        
          Oocyst_data(ix1,ix2,i,j) = o(ix1,ix2) - 1
        
        }else if((choose_transition > transition_type_indices(6)) && (choose_transition <= transition_type_indices(7))){
  
          #% Mating 6
          #CT_mat = reshape(transition_type_indices(6)+1:transition_type_indices(7),N,N);
          C(ix1,ix2) = find(CT_mat == choose_transition)
        
          Male_data(ix2,i,j) = m(ix2) - 1
          Female_data(ix1,i,j) = f(ix1) - 1
          Zygote_data(ix1,ix2,i,j) = z(ix1,ix2) + 1
        
        }else if((choose_transition > transition_type_indices(7)) && (choose_transition <= transition_type_indices(8))){
          
          #% Zygote maturation 7
          CT_mat = reshape(transition_type_indices(7)+1:transition_type_indices(8),N,N)
          C(ix1,ix2) = find(CT_mat == choose_transition)
        
          Zygote_data(ix1,ix2,i,j) = z(ix1,ix2) - 1
          Ookinete_data(ix1,ix2,i,j) = e(ix1,ix2) + 1
        
        }else if((choose_transition > transition_type_indices(8)) && (choose_transition <= transition_type_indices(9))){
            
          #% Ookinete maturation 8
          CT_mat = reshape(transition_type_indices(8)+1:transition_type_indices(9),N,N)
          C(ix1,ix2) = find(CT_mat == choose_transition)
        
          Ookinete_data(ix1,ix2,i,j) = e(ix1,ix2) - 1
          Oocyst_data(ix1,ix2,i,j) = o(ix1,ix2) + 1
        
        }else if((choose_transition > transition_type_indices(9)) && (choose_transition <= transition_type_indices(10))){
          
          #% Oocyst bursting 9
          #CT_mat = reshape(transition_type_indices(9)+1:transition_type_indices(10),N,N);
          C(ix1,ix2) = find(CT_mat == choose_transition)
        
          Oocyst_data(ix1,ix2,i,j) = o(ix1,ix2) - 1;
          Sporozoite_data(ix1,ix2,i,j) = s(ix1,ix2) + my_binornd(poissrnd(n),p)
          Burst_time(ix1,ix2,i,j) = time1
          burst_count(ix1,ix2,j) = burst_count(ix1,ix2,j) + 1
        
        }

Time_data(1:length(t),j) = t
        
rm(t)

#TOC

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
        
#%% Uncomment to Overwrite Data in Directory RawData/ContRuptFunc ---------- %
        
filename = paste0('RawData/ContRuptFunc/nParasiteGroupsData_G0',
           num2str(G0_vec(ii)),
           '_max_bias',
           num2str(max_bias),
           '_NumStrains',
           num2str(N),
           '.mat')

save(filename,
     'Time_data',
     'Male_data',
     'Female_data',
     'Zygote_data',
     'Ookinete_data',
     'Oocyst_data',
     'Sporozoite_data',
     'Burst_time',
     'burst_count',
     'TransitionIDs',
     'SEED')
      
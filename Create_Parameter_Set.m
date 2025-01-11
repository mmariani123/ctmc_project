% Function creates structure of parameter values.  Number of genotypes,
% simulations, and max fitness bias can also be modified here.
function param = Create_Parameter_Set


param.N = 2; % Number of different starting genotypes

param.a = 24/(20/60); % failure rate of male gametocytes (per day)
param.b = 24/(25/60); % failure rate of female gametocytes (per day)
param.r = 0.08; % fertilization of male and female gametes (per day)
param.mu_z = 1; % death rate of zygotes (per day)
param.sigma_z = 24/19; % transformation rate of zygotes (per day)
param.sigma_e = 0.6; % transformation rate of ookinetes (per day)
param.mu_e = 1.4; % death rate of ookinetes (per day)
param.mu_o = 0; % death rate of oocysts (per day) **changed from paper**
param.n = 3e3; % number of sporozoites per oocyst
param.p = 0.2; % proportion of sporozoites that make it to salivary gland
param.alpha = 0.39; % fraction of male gametes that are viable
param.beta = 1;%0.96; % fraction of female gametes that are viable
param.eta = 1; % number of female gametes per female gametocyte
param.rho = 8; % number of male gametes per male gametocytes

param.k = 1/7;
param.t0 = 10;

param.max_bias = 0.1;%0.5;%0.1;

param.NumSim = 1;%10000;

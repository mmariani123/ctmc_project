% -- Simulate CTMC model of within-vector parasite dynamics using modified
% Gillespie Algorithm (outlined in Supplementary Materials), with
% continuous rupture function.
function SpDat=CTMC_ContRuptFunc(par)
% close all
% clear
%
% SEED = rng('shuffle');% set seed
%
% % Set parameters equal to values defined in Create_Parameter_Set.m:
% par = Create_Parameter_Set;

N = par.N;
a = par.a;
b = par.b;
r = repmat(par.r,N,N);
mu_z = par.mu_z;
sigma_z = par.sigma_z;
sigma_e = par.sigma_e;
mu_e = par.mu_e;
mu_o = par.mu_o;
n = par.n;
p = par.p;
alpha = par.alpha;
beta = par.beta;
eta = par.eta;
rho = par.rho;
k = par.k;
max_bias = par.max_bias;
bias = par.bias;
t0 = par.t0;
male_prop = par.male_prop;

percent_male = repmat(male_prop,N,1); % proportion of gametocytes that are male of each genotype

% Define vector of fitness bias for each genotype:
% if N==1;
%     bias = 0;
% elseif N>1
%     bias_tmp = zeros(1,N);%0:1/(N-1):1;
%     bias = bias_tmp*max_bias;
% send


bias = bias';
a_vec = repmat(a,N,1).*(1 + bias);
b_vec = repmat(b,N,1).*(1 + bias);
alpha_vec = repmat(alpha,N,1).*(1 - bias);
beta_vec = repmat(beta,N,1).*(1 - bias);
mu_z_vec = repmat(mu_z,N,1).*(1 + bias); % zygote mortality rate for zygotes whose parents are of the same genotype
mu_e_vec = repmat(mu_e,N,1).*(1 + bias); % ookinete mortality rate for ookinetes " ...
mu_o_vec = repmat(mu_o,N,1).*(1 + bias); % oocyst mortality rate oocysts " ...


% Create symmetric matrices of biased parameters.  The (i,j) entry is the
% parameter value for a parasite that has a genotype (i) female, and
% genotype (j) male parent.

mu_z_mat = NaN(N,N);
mu_e_mat = NaN(N,N);
mu_o_mat = NaN(N,N);

for i = 1:N;
    for j = 1:N;

        mu_z_mat(i,j) = mean([mu_z_vec(i),mu_z_vec(j)]);
        mu_e_mat(i,j) = mean([mu_e_vec(i),mu_e_vec(j)]);
        mu_o_mat(i,j) = mean([mu_o_vec(i),mu_o_vec(j)]);
    end
end



% --------------------------------------------------------------------%
% ---- Run simulation for values defined in CreateParameterSet.m:

% Max fitness bias %
max_bias = par.max_bias;
% Max simulation time %
Tfinal = par.Tfinal;
% Max number of iterations for each simulation %
maxiter = par.maxiter;

% Time-step to increment by if the time to next time-step is greater than
% dt.
dt = 0.1;

% Fraction of G0 of each subtype (Nx1 vector):
strainprop = repmat(1/N,N,1);

%TransitionIDs = NaN(maxiter);



G0 = round(strainprop*par.G0); % Vector of initial gametocyte densities.  G0(i) is the # of genotype i gametocytes.

Gf0 = ceil((1 - percent_male).*G0.*eta.*beta); % Female gamete densities of each genotype
Gm0 = ceil(percent_male.*G0.*rho.*alpha); % Male gamete densities of each genotype

Time_data = NaN(maxiter,1); % Initialize Time_data vector.

Male_data = zeros(N,maxiter); % Each row is a different genotype; each column is a different time step; each layer is a different simulation
Female_data = zeros(N,maxiter); % Same as for Male_data
Zygote_data = zeros(N,N,maxiter); % 4D array: Zygote(i,j,k,l) is the # of Zygotes with female parent i, male parent j, at timestep k, in simulation l
Ookinete_data = zeros(N,N,maxiter); % "
Oocyst_data = zeros(N,N,maxiter); % "
Sporozoite_data = zeros(N,N,maxiter); % "
Burst_time = zeros(N,N,maxiter); % "

SpDat = zeros(1,4); % Initialize Matrix whose columns are (1) female parent ID, (2) male parent ID, (3) Time of oocyst rupture, and (4) number of new sporozoites

%tic

Male_data(:,1) = Gm0; % Set initial conditions for male and female gamete data.
Female_data(:,1) = Gf0; % "


burst_count(:,:) = zeros(N,N); % Initialize burst_count 3-D array.

t(1) = 0;
time = t(1);
    transition_type_indices = [0;cumsum([N;N;N*N;...
        N*N;N*N;N*N;...
        N*N;N*N;N*N])];

    j = 0;
for i = 2:maxiter

    y1 = rand;
    y2 = rand;

    % State variables
    m = Male_data(:,i-1); % Column vector where each row is a different genotype, and the elements are the # of Male gametes at the previous time-step (t_{i-1}) in simulation j.  So, m(k) is the number of genotype-k males at time t_{i-1} in sim j.
    f = Female_data(:,i-1);
    z = Zygote_data(:,:,i-1); % z(k,l) = # zygotes with female-k parent and male-l parent at time t_{i-1} in sim j.
    e = Ookinete_data(:,:,i-1);
    o = Oocyst_data(:,:,i-1);
    s = Sporozoite_data(:,:,i-1);

    % Break if states m through o are extinct:
    state_variables = [m;f;reshape(z,N*N,1);...
        reshape(e,N*N,1);reshape(o,N*N,1)];
    if max(state_variables(:)) == 0
        LastIter = i-1;
        break
    end

    % Assign data values at current time-step, i, to values at previous time-step, i-1.
    Male_data(:,i) = m;
    Female_data(:,i) = f;
    Zygote_data(:,:,i) = z;
    Ookinete_data(:,:,i) = e;
    Oocyst_data(:,:,i) = o;
    Sporozoite_data(:,:,i) = s;

    % Define bursting function.
    bursting = 1/(1 + (exp(t0 - t(i-1)))); % probability of bursting at time 'time'
    % Create vectors/matrices of all possible transitions, arranged by type: death, mating, maturation, bursting.
    Male_death = a_vec.*m;
    Female_death = b_vec.*f;
    Zygote_death = mu_z_mat.*z;
    Ookinete_death = mu_e_mat.*e;
    Oocyst_death = mu_o_mat.*o;
    Mating = r.*(f*m');
    Zyg_maturation = sigma_z*z;
    Ook_maturation = sigma_e*e;
    Ooc_Bursting = k.*bursting.*o;


    % Calculate all possible transitions:
    transitions = [Male_death;Female_death;reshape(Zygote_death,N*N,1);...
        reshape(Ookinete_death,N*N,1);reshape(Oocyst_death,N*N,1);reshape(Mating,N*N,1);...
        reshape(Zyg_maturation,N*N,1);reshape(Ook_maturation,N*N,1);reshape(Ooc_Bursting,N*N,1)];


    % Create vector of indices marking different transition types in
    % the transition vector:
    %transition_type_indices = [0;cumsum([length(Male_death);length(Female_death);length(reshape(Zygote_death,N*N,1));...
    %    length(reshape(Ookinete_death,N*N,1));length(reshape(Oocyst_death,N*N,1));length(reshape(Mating,N*N,1));...
    %    length(reshape(Zyg_maturation,N*N,1));length(reshape(Ook_maturation,N*N,1));length(reshape(Ooc_Bursting,N*N,1))])];

%     transition_type_indices = [0;cumsum([N;N;N*N;...
%         N*N;N*N;N*N;...
%         N*N;N*N;N*N])];


    % Calculate vector of transition probabilities %
    total_transitions = sum(transitions);
    trans_probs = transitions/total_transitions;

    % Calculate time of next event.  If time > Tfinal, break!
    stoch_time_step = -log(y1)/(total_transitions);
    time1 = stoch_time_step +t(i-1); % time until next event;
    time2 = t(i-1) + dt;

    if min(time1,time2) > Tfinal

        LastIter = i-1;

        break

    elseif stoch_time_step > dt
        t(i) = time2;
    elseif stoch_time_step <= dt
        t(i) = time1;


        % ---------------------- %
        % Determine next event and update population matrices %

        choose_transition = find(cumsum(trans_probs)>y2,1,'first');


        %TransitionIDs(i) = choose_transition;

        if choose_transition <= transition_type_indices(2);
            % Male gamete death 1
            CT = choose_transition - transition_type_indices(1);
            Male_data(CT,i) = m(CT) - 1;

        elseif choose_transition > transition_type_indices(2) && choose_transition <= transition_type_indices(3)
            % Female gamete death 2
            CT = choose_transition - transition_type_indices(2);
            Female_data(CT,i) = f(CT) - 1;

        elseif choose_transition > transition_type_indices(3) && choose_transition <= transition_type_indices(4)
            % Zygote Death 3
            CT_mat = reshape(transition_type_indices(3)+1:transition_type_indices(4),N,N);
            [ix1,ix2] = find(CT_mat == choose_transition);

            Zygote_data(ix1,ix2,i) = z(ix1,ix2) - 1;

        elseif choose_transition > transition_type_indices(4) && choose_transition <= transition_type_indices(5)
            % Ookinete Death 4
            CT_mat = reshape(transition_type_indices(4)+1:transition_type_indices(5),N,N);
            [ix1,ix2] = find(CT_mat == choose_transition);

            Ookinete_data(ix1,ix2,i) = e(ix1,ix2) - 1;

        elseif choose_transition > transition_type_indices(5) && choose_transition <= transition_type_indices(6)
            % Oocyst Death 5
            CT_mat = reshape(transition_type_indices(5)+1:transition_type_indices(6),N,N);
            [ix1,ix2] = find(CT_mat == choose_transition);

            Oocyst_data(ix1,ix2,i) = o(ix1,ix2) - 1;

        elseif choose_transition > transition_type_indices(6) && choose_transition <= transition_type_indices(7)
            % Mating 6
            CT_mat = reshape(transition_type_indices(6)+1:transition_type_indices(7),N,N);
            [ix1,ix2] = find(CT_mat == choose_transition);

            Male_data(ix2,i) = m(ix2) - 1;
            Female_data(ix1,i) = f(ix1) - 1;
            Zygote_data(ix1,ix2,i) = z(ix1,ix2) + 1;

        elseif choose_transition > transition_type_indices(7) && choose_transition <= transition_type_indices(8)
            % Zygote maturation 7
            CT_mat = reshape(transition_type_indices(7)+1:transition_type_indices(8),N,N);
            [ix1,ix2] = find(CT_mat == choose_transition);

            Zygote_data(ix1,ix2,i) = z(ix1,ix2) - 1;
            Ookinete_data(ix1,ix2,i) = e(ix1,ix2) + 1;

        elseif choose_transition > transition_type_indices(8) && choose_transition <= transition_type_indices(9)
            % Ookinete maturation 8
            CT_mat = reshape(transition_type_indices(8)+1:transition_type_indices(9),N,N);
            [ix1,ix2] = find(CT_mat == choose_transition);

            Ookinete_data(ix1,ix2,i) = e(ix1,ix2) - 1;
            Oocyst_data(ix1,ix2,i) = o(ix1,ix2) + 1;

        elseif choose_transition > transition_type_indices(9) && choose_transition <= transition_type_indices(10)
            % Oocyst bursting 9
            CT_mat = reshape(transition_type_indices(9)+1:transition_type_indices(10),N,N);
            [ix1,ix2] = find(CT_mat == choose_transition);

            Oocyst_data(ix1,ix2,i) = o(ix1,ix2) - 1;
            new_sporo = my_binornd(poissrnd(n),p);
            Sporozoite_data(ix1,ix2,i) = s(ix1,ix2) + new_sporo;
            Burst_time(ix1,ix2,i) = time1;
            burst_count(ix1,ix2) = burst_count(ix1,ix2) + 1;

            j = j+1;
            SpDat(j,:) = [ix1,ix2,time1,new_sporo]; % ix1 = female ID, ix2 = male ID

        end

    end


    Time_data(1:length(t)) = t;

end

%toc

% --- Compute input for Diversity model:
%   (1) Initial male and female gamete densities.
%   (2) Total ruptured oocysts by day Tfinal.
%   (3) Total number of sporozoites surviving to salivary glands by day
%   Tfinal.

Tfinal_index = find(Time_data(:) <= Tfinal,1,'last');

T_index = NaN(1,Tfinal);
BT_tmp1 = Burst_time>0;
BurstCount = NaN(N,N,Tfinal);


for i = 1:Tfinal
     T_index(i) = find(Time_data(:) <= i,1,'last');
     [t1,t2,t3]=size(sum(BT_tmp1(:,:,1:T_index(i)),3));
     [s1,s2,s3]=size(BurstCount(:,:,i));
     %if t1~=s1 & t2~=s2 & t3~=s3
         %disp([t1,t2,t3])
         %disp([s1,s2,s3])
     %end
     BurstCount(:,:,i) = sum(BT_tmp1(:,:,1:T_index(i)),3);
 end


 FGam = Female_data(:,1);
 MGam = Male_data(:,1);
 Sporo = Sporozoite_data(:,:,T_index);


 S1 = reshape(permute(Sporo,[2,1,3]),1,N^2,Tfinal);
 [x y z] = size(S1);
 for j=1:z
     S2(j,1:y) = S1(1,:,j);
 end

 BC1 = reshape(permute(BurstCount,[2,1,3]),1,N^2,Tfinal);
 [x y z] = size(BC1);
 for j=1:z
     BC2(j,1:y) = BC1(1,:,j);
 end

 % --- Uncomment to create plots of cumulative sporos and ruptured oocysts:
  BC11 = reshape(BurstCount(1,1,:),Tfinal,1);
  BC12 = reshape(BurstCount(1,2,:),Tfinal,1);
  BC21 = reshape(BurstCount(2,1,:),Tfinal,1);
  BC22 = reshape(BurstCount(2,2,:),Tfinal,1);

  S11 = reshape(Sporo(1,1,:),Tfinal,1);
  S12 = reshape(Sporo(1,2,:),Tfinal,1);
  S21 = reshape(Sporo(2,1,:),Tfinal,1);
  S22 = reshape(Sporo(2,2,:),Tfinal,1);

  figure
  hold all
  stairs(1:Tfinal,BC11)
  stairs(1:Tfinal,BC12)
  stairs(1:Tfinal,BC21)
  stairs(1:Tfinal,BC22)

  figure
  hold all
  stairs(1:Tfinal,S11)
  stairs(1:Tfinal,S12)
  stairs(1:Tfinal,S21)
  stairs(1:Tfinal,S22)



 %     % Uncomment to Overwrite Data in Directory RawData/ContRuptFunc ------------ %
 %
 %         filename = ['RawData/ContRuptFunc/nParasiteGroupsData_G0' ...
 %             num2str(G0_vec(ii)) '_max_bias' num2str(max_bias) '_NumStrains' num2str(N) '.mat'];
 %         save(filename,'Time_data','Male_data','Female_data','Zygote_data',...
 %             'Ookinete_data','Oocyst_data','Sporozoite_data','Burst_time','burst_count','TransitionIDs','SEED')
 %





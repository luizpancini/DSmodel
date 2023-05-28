function [bestGlobalAim,bestGlobalPos,P,V] = BLPSO_speedcore(BL_INPUTS,PSO_options)

% Unpack PSO algorithm options
[adaptiveInertia,adaptiveSwarmSize,alwaysUpdateInf,c0,c1,c2,D,EM_swarms,initializePos,initializeVel,improveTol,maxFE,maxStalledIt,max_c0,min_c0,moveMethod,N,N_c,N_i,N_m,N_e,num_workers_parfor,reboundVel,run_parfor,set_GP,specified_filename,stopOnMaxStalledIt,setNMLOnMaxStalledIt,topology] = unpack_BLPSO_options(PSO_options);

% Bounds according to airfoil
[LB,UB,LBmat,UBmat,Vmax,Vmaxmat] = set_BLPSO_bounds(N,BL_INPUTS.airfoil);

% Initialize arrays
[Aim,bestAimIter,bestGlobalAim,bestOwnAim,bestOwnPos,bestInformantsAim,bestInformantsPos,bestMemInformantsAim,bestMemInformantsPos,bestExpInformantsAim,bestExpInformantsPos,informants,cumStallIt,FE,moveMethodFlag,N_c0,N_i0,N_max,N_min,N_e_range,N_m_range,noiseFlag,P,V,topologyFlag,axes1,axes_size,lw,ms] = ...
 initialize_BLPSO(D,EM_swarms,initializePos,initializeVel,LB,maxFE,moveMethod,N,N_c,N_i,N_e,N_m,num_workers_parfor,topology,UB,Vmax,run_parfor,set_GP,specified_filename);

% Initialize iteration and stalled iterations counters
it = 0; stalledIt = 0;   

% Loop until maximum number of function evaluations is reached
while FE < maxFE 
    
    % Update iteration counter
    it = it+1;
    
    %% Find aims (function values)
    % Loop over particles
    if run_parfor
        parfor p=1:N
            Aim(p) = optimize_BL(P(p,:),BL_INPUTS);
        end
    else
        for p=1:N
            Aim(p) = optimize_BL(P(p,:),BL_INPUTS);
        end
    end
    % Update number of function evaluations
    FE = FE + N;
    
    %% Update best positions and aims
    % TFs for each particle having improved its position with respect to the last iteration
    newOwnBest = Aim < bestOwnAim;
    oldOwnBest = true(N,1) - newOwnBest;
    particlesToUpdateInf = oldOwnBest + alwaysUpdateInf*newOwnBest;
    % Update particle's best personal position and aim
    bestOwnPos(newOwnBest,:) = P(newOwnBest,:);
    bestOwnAim(newOwnBest) = Aim(newOwnBest);
    % Update global best position and aim
    prevBestGlobalAim = bestGlobalAim;
    [bestGlobalAim,bestParticleNow] = min(bestOwnAim);
    bestGlobalPos = bestOwnPos(bestParticleNow,:);
    % TF for best global aim being improved with respect to previous
    % iteration (given the tolerance for improvement)
    improvedBestGlobalAim = bestGlobalAim < (1-improveTol)*prevBestGlobalAim;
    % Update stalled iterations counter
    if ~improvedBestGlobalAim
        stalledIt = stalledIt + 1;
    else
        stalledIt = 0;
    end
    % Update array for current iteration information on best aim
    bestAimIter(it,:) = [bestGlobalAim, bestParticleNow];
        
    %% Stalled number of iterations check / No Man's Land reset
    if stalledIt > maxStalledIt 
        if stopOnMaxStalledIt           % Stop
            return;
        elseif setNMLOnMaxStalledIt     % Set a random particle into the "No Man's Land"
            % Select a random particle
            p = randi(N);
            % For each dimension
            for d=1:D
                % Sort positions of all particles in that dimension, including the borders
                sorted_Pd = [LB(d); sort(P(:,d),'ascend'); UB(d)];
                % Find index where there is the greatest gap
                [~,ind_max_diff] = max(diff(sorted_Pd));
                % Set particle's position in the center of the gap
                P(p,d) = mean([sorted_Pd(ind_max_diff) sorted_Pd(ind_max_diff+1)]);
            end
            % Calculate aim and update FE
            Aim(p) = optimize_BL(P(p,:),BL_INPUTS);
            FE = FE + 1;
        end
    end
    
    %% Adaptation
    N_changed = 0;  % TF for the swarm size having changed this iteration
    % Adapt swarm size - add a particle
    if adaptiveSwarmSize && N < N_max && stalledIt >= maxStalledIt
        % Update TF
        N_changed = 1;
        % Find the position in the center of the current "No Man's Land" (center of greatest gap in each dimension)
        P_nml = nan(1,D);
        % For each dimension
        for d=1:D
            % Sort positions of all particles in that dimension, including the borders
            sorted_Pd = [LB(d); sort(P(:,d),'ascend'); UB(d)];
            % Find index where there is the greatest gap
            [~,ind_max_diff] = max(diff(sorted_Pd));
            % Set new particle's position in the center of the gap
            P_nml(d) = mean([sorted_Pd(ind_max_diff) sorted_Pd(ind_max_diff+1)]);
        end
        % Calculate aim and update FE
        Aim_Pnml = optimize_BL(P_nml,BL_INPUTS);
        FE = FE + 1;
        % Update / initialize swarm size dependent variables
        N = N+1;
        if ~EM_swarms % Single swarm
            % Add particle at the end of arrays
            indNewPart = N;
            % Initialize
            bestInformantsAim(N,1) = 1e9;
            bestInformantsPos(N,:) = zeros(1,D);
        else 
            % Add particle to the smaller swarm
            if N_m <= N_e   % Add a memory
                % Add particle at the begin of arrays
                indNewPart = 1;
                % Update
                N_m = N_m+1;
                N_m_range = 1:N_m; 
                N_e_range = N_m+1:N;
                % Shift arrays
                P(2:N,:) = P;
                V(2:N,:) = V;
                Aim(2:N,1) = Aim;
                bestOwnAim(2:N,1) = bestOwnAim;
                bestOwnPos(2:N,:) = bestOwnPos;
                newOwnBest(2:N,1) = newOwnBest;
                oldOwnBest(2:N,1) = oldOwnBest;
                particlesToUpdateInf(2:N,1) = particlesToUpdateInf;
                bestMemInformantsAim(2:N_m,1) = bestMemInformantsAim;
                bestMemInformantsPos(2:N_m,:) = bestMemInformantsPos;
                % Initialize
                bestMemInformantsAim(1,1) = 1e9;
                bestMemInformantsPos(1,:) = zeros(1,D);
            else            % Add an explorer
                % Add particle at the end of arrays
                indNewPart = N;
                % Update
                N_e = N_e+1;
                N_m_range = 1:N_m; 
                N_e_range = N_m+1:N;
                % Initialize
                bestExpInformantsAim(N_e,1) = 1e9;
                bestExpInformantsPos(N_e,:) = zeros(1,D);
            end
        end
        % Update variables for new particle
        Aim(indNewPart,1) = Aim_Pnml;
        P(indNewPart,:) = P_nml;
        V(indNewPart,:) = zeros(1,D);
        bestOwnAim(indNewPart) = Aim_Pnml;
        bestOwnPos(indNewPart,:) = P_nml;
        newOwnBest(indNewPart,1) = 1;
        oldOwnBest(indNewPart,1) = 0;
        particlesToUpdateInf(indNewPart,1) = 1;
        % Update other swarm size dependent variables
        LBmat(N,:) = LB;
        UBmat(N,:) = UB;
        Vmaxmat(N,:) = Vmaxmat(1,:);
    end
    % Adapt swarm size - remove a particle
    if adaptiveSwarmSize && N > N_min && nnz(newOwnBest) >= N/2 && it > 1 % If the number of particles is greater than the number of dimensions, and if at least half have improved their positions      
        % Update TF
        N_changed = 1;
        % Find current worst particle
        [~,worstParticleNow] = max(bestOwnAim);
        % Update swarm size dependent variables
        N = N-1;
        Aim(worstParticleNow) = [];
        P(worstParticleNow,:) = [];
        V(worstParticleNow,:) = [];
        bestOwnAim(worstParticleNow) = [];
        bestOwnPos(worstParticleNow,:) = [];
        newOwnBest(worstParticleNow) = [];
        oldOwnBest(worstParticleNow) = [];
        particlesToUpdateInf(worstParticleNow) = [];
        if ~EM_swarms % Single swarm
            bestInformantsAim(worstParticleNow) = [];
            bestInformantsPos(worstParticleNow,:) = [];
        else
            % Determine if worst is an explorer or memory and update accordingly
            if worstParticleNow <= N_m  % Is a memory
                N_m = N_m-1;
                N_m_range = 1:N_m; 
                N_e_range = N_m+1:N; 
                ind = worstParticleNow; % Index for the particle in the memory swarm
                bestMemInformantsAim(ind) = [];
                bestMemInformantsPos(ind,:) = [];
            else                        % Is an explorer
                N_e = N_e-1;
                N_m_range = 1:N_m; 
                N_e_range = N_m+1:N;
                ind = worstParticleNow-N_m; % Index for the particle in the explorer swarm
                bestExpInformantsAim(ind) = [];
                bestExpInformantsPos(ind,:) = [];
            end
        end
        % Update other swarm size dependent variables
        LBmat(1,:) = [];
        UBmat(1,:) = [];
        Vmaxmat(1,:) = [];
    end
    % Adapt number of informants
    if topologyFlag == 2 % topology = "adaptive_N_randi"
        if ~improvedBestGlobalAim  % No improvement
            % Enhance communication
            if ~EM_swarms
                N_i = min(N,N_i+2);     % Increase number of informants
            else
                N_c = max(1,N_c-2);     % Decrease number of communication channels
            end
        else            % Improvement
            % Decrease communication
            if ~EM_swarms
                N_i = max(N_i0,N_i-1);  % Decrease number of informants
            else
                N_c = min(N_c0,N_c+1);  % Increase number of communication channels
            end
        end
    end
    % Adapt inertia weight
    if adaptiveInertia
        % Update cummulative counter for no improvement iterations
        if improvedBestGlobalAim
            cumStallIt = max(0,cumStallIt-1); % Decrease counter
        else
            cumStallIt = cumStallIt+1;        % Increase counter
        end
        % Update inertia weight according to counter, enforcing limits
        if cumStallIt < 2
            c0 = max(min_c0,min(max_c0,1.1*c0));
        elseif cumStallIt > 5
            c0 = max(min_c0,min(max_c0,0.75*c0));
        end
    end
    
    %% Topology (method of communication with informants)
    % Topology for single swarm
    if ~EM_swarms
        switch topologyFlag
            case {1,2} % topology = "fixed_N_randi" or topology = "adaptive"
                % For each particle, select N_i random informants, including the particle itself
                informants = [(1:N)' randi(N,N,N_i-1)]; % Does not guarantee N_i unique informants, but is much faster than the method below
                % For each particle, select exactly N_i unique random informants, including the particle itself
%                 informants = zeros(N,N_i);  % Initialize
%                 informants(:,1) = 1:N;      % First informant is self
%                 for p=1:N
%                     % Determine unique random remaining informants
%                     rem_inf = randperm(N-1,N_i-1);
%                     % Add 1 to indices that are >= current particle index (so that the last can be an informant and the particle itself is not chosen again)
%                     rem_inf(rem_inf >= p) = rem_inf(rem_inf >= p) + 1;
%                     % Set remaining informants
%                     informants(p,2:end) = rem_inf;
%                 end
            case {3,4}      % topology = "ring" or topology = "SL_ring"
                if it == 1 || N_changed  % Fixed informants, so only need to be set once in case the swarm size is also fixed
                    % Relative indices of informants
                    if topologyFlag == 3        % [topology = "ring"] For each particle, informants are its neighbors by index in a ring, besides itself
                        ind = [-1 1];
                    elseif topologyFlag == 4    % [topology = "SL_ring"] Particle p informs itself and particles p+1 and p-2, and is informed by particles p-1 and p+2
                        ind = [-1 2];
                    end
                    % Informants
                    informants = zeros(N,3);  % Initialize
                    informants(:,1) = 1:N;    % First informant is self
                    % Second column has informants in the counter-clockwise direction,
                    % third column has informants in the clockwise direction
                    for k=1:length(ind)
                        informants(:,k+1) = circshift(informants(:,1),-ind(k));
                    end
                end
            case 5  % topology == "VR_ring"
                % Randomly permuted informants in a ring
                informants = zeros(N,3);  % Initialize
                informants(:,1) = 1:N;    % First informant is self
                seq = randperm(N,N);      % New sequence, running in the clockwise direction
                % Second column has informants in the counter-clockwise direction, 
                % third column has informants in the clockwise direction
                informants(seq,2:3) = [circshift(seq,1)', circshift(seq,-1)']; 
        end
        % For each particle, sort informants by their best personal aim and get the best informant's position
        [bestNewInformantsAim,bestNewInformantsInd] = min(bestOwnAim(informants),[],2);
        bestNewInformantsPos = bestOwnPos(informants(sub2ind([N,N_i],1:N,bestNewInformantsInd')),:);
        % Update informants according to specifications (only if position has not improved, or always)
        bestInformantsPos = bestNewInformantsPos .* particlesToUpdateInf + bestInformantsPos .* (1-particlesToUpdateInf);
        bestInformantsAim = bestNewInformantsAim .* particlesToUpdateInf + bestInformantsAim .* (1-particlesToUpdateInf);
    % Topology for memory/explorer swarms
    else
        % For each particle in the swarms, determine a random channel of communication
        m_channel = randi(N_c,N_m,1);
        e_channel = randi(N_c,N_e,1);
        % Initialize arrays for best informants in each channel
        bestMemPosOnCh = nan(N_c,D);    % Position of the best memories in each channel
        bestMemAimOnCh = nan(N_c,1);    % Aims of the best memories in each channel
        bestExpPosOnCh = nan(N_c,D);    % Position of the best explorers in each channel
        bestExpAimOnCh = nan(N_c,1);    % Aims of the best explorers in each channel
        % Loop over channels
        for c = 1:N_c
            % Memories
            memOnCh = N_m_range(m_channel==c);                          % List of memories in the current channel
            if isempty(memOnCh)                                         % In case no memory tunned to the current channel
                mostTunnedCh = mode(m_channel);                         % Find the most tunned channel by the memories
                ind = find(m_channel==mostTunnedCh,1,'first');          % And first corresponding memory that tunned to that channel
                memOnCh = N_m_range(ind);                               % Set that explorer to the current channel
                m_channel(ind) = c;                                     % And realize the change on the channel list
            end
            [bestMemAimOnCh(c),bestMemInd] = min(bestOwnAim(memOnCh));  % Aim and index of the best memory in current channel
            bestMemPosOnCh(c,:) = bestOwnPos(memOnCh(bestMemInd),:);    % Position of the best memory in current channel
            % Explorers
            expOnCh = N_e_range(e_channel==c);                          % List of explorers in the current channel
            if isempty(expOnCh)                                         % In case no explorer tunned to the current channel
                mostTunnedCh = mode(e_channel);                         % Find the most tunned channel by the explorers
                ind = find(e_channel==mostTunnedCh,1,'first');          % And first corresponding explorer that tunned to that channel
                expOnCh = N_e_range(ind);                               % Set that explorer to the current channel
                e_channel(ind) = c;                                     % And realize the change on the channel list
            end
            [bestExpAimOnCh(c),bestExpInd] = min(bestOwnAim(expOnCh));  % Aim and index of the best explorer in current channel
            bestExpPosOnCh(c,:) = bestOwnPos(expOnCh(bestExpInd),:);    % Position of the best explorer in current channel
        end
        % For each explorer and memory, pick the best informant from the other swarm available on its current channel
        bestNewExpInformantsPos = bestMemPosOnCh(e_channel,:);
        bestNewMemInformantsPos = bestExpPosOnCh(m_channel,:);
        bestNewExpInformantsAim = bestMemAimOnCh(e_channel,:);
        bestNewMemInformantsAim = bestExpAimOnCh(m_channel,:);
        % Particles in each swarm to update informants
        expToUpdateInf = newOwnBest(N_e_range) + alwaysUpdateInf*oldOwnBest(N_e_range);
        memToUpdateInf = newOwnBest(N_m_range) + alwaysUpdateInf*oldOwnBest(N_m_range);
        % Update informants position and aims according to specifications
        bestExpInformantsPos = bestNewExpInformantsPos .* expToUpdateInf + bestExpInformantsPos .* (1-expToUpdateInf);
        bestMemInformantsPos = bestNewMemInformantsPos .* memToUpdateInf + bestMemInformantsPos .* (1-memToUpdateInf);
        bestExpInformantsAim = bestNewExpInformantsAim .* expToUpdateInf + bestExpInformantsAim .* (1-expToUpdateInf);
        bestMemInformantsAim = bestNewMemInformantsAim .* memToUpdateInf + bestMemInformantsAim .* (1-memToUpdateInf);
    end
    
    %% Movement (update velocities and set new positions)
    % Explorer/Memory swarms' variables
    if EM_swarms
        % Memory swarm
        Pm = P(N_m_range,:);                    % Memories' positions
        bestMemPos = bestOwnPos(N_m_range,:);   % Memories' best personal positions
        bestMemAim = bestOwnAim(N_m_range,:);   % Memories' best personal aims
        % Explorer swarm
        Pe = P(N_e_range,:);                    % Explorers' positions
        Ve = V(N_e_range,:);                    % Explorers' velocities
    end    
    % Update velocities and positions
    switch moveMethodFlag
        % Uniform random distribution (standard PSO)
        case 1 % move_method = "uniform"          
            if ~EM_swarms
                V = c0*V + c1*rand(N,D).*(bestOwnPos-P) + c2*rand(N,D).*(bestInformantsPos-P);
                P = P + V;
            else
                % Memory swarm
                Vm = c1*rand(N_m,D).*(bestMemPos-Pm) + c2*rand(N_m,D).*(bestMemInformantsPos-Pm);
                Pm = Pm + Vm;
                % Explorer swarm
                Ve = c0*Ve + c2*rand(N_e,D).*(bestExpInformantsPos-Pe);
                Pe = Pe + Ve;
            end
        % Clerc's Gaussian pivots
        case 2 % move_method = "pivots" or move_method = "noisy_pivots"
            if ~EM_swarms  % Single swarm
                % Position beyond current one, in the direction of the best personal
                Pp = P + c1*(bestOwnPos-P);
                % Position beyond current one, in the direction of the best informant
                Pg = P + c2*(bestInformantsPos-P);
                % Barycenter of current, personal best and best informant's positions (resumes to that of only current and personal best if personal and informant's best are coincident)
                G = 1/3*(Pp+Pg+P) + 1/6*(Pp+P-2*Pg).*(bestOwnAim == bestInformantsAim);
                % Radii of hyperspheres
                rho = vecnorm(G-P,2,2);
                % Noise
                if noiseFlag
                    % Weigth factors on particles' personal best position
                    cp = bestOwnAim./(bestOwnAim+bestInformantsAim);
                    % Weigth factors on best informant's position
                    cg = bestInformantsAim./(bestOwnAim+bestInformantsAim);
                    % Standard deviation for noise
                    sigma = cp-cg;
                    % If the particle's aim improved, allow Gaussian noise
                    noise = sigma.*randn(N,1).*(newOwnBest);
                else
                    noise = 0;
                end
                % Update velocity
                V = c0*V + G + rand_sphere(N,D,rho) - P;
                % Update position
                P = P + (1+noise).*V;
            else
                % Memory swarm
                Pp = Pm + c1*(bestMemPos-Pm);
                Pg = Pm + c2*(bestMemInformantsPos-Pm);
                G = 1/3*(Pp+Pg+Pm) + 1/6*(Pp+Pm-2*Pg).*(bestMemAim == bestMemInformantsAim);
                rho = vecnorm(G-Pm,2,2);
                if noiseFlag
                    cp = bestMemAim./(bestMemAim+bestMemInformantsAim);
                    cg = bestMemInformantsAim./(bestMemAim+bestMemInformantsAim);
                    sigma = cp-cg;
                    noise = sigma.*randn(N_m,1).*(newOwnBest(N_m_range));
                else
                    noise = 0;
                end
                Vm = G + rand_sphere(N_m,D,rho) - Pm;
                Pm = Pm + (1+noise).*Vm;
                % Explorer swarm
                Pg = Pe + c2*(bestExpInformantsPos-Pe);
                G = 1/2*(Pg+Pe);
                rho = vecnorm(G-Pe,2,2);
                if noiseFlag
                    noise = randn(N_e,1).*(newOwnBest(N_e_range));
                else
                    noise = 0;
                end
                Ve = c0*Ve + G + rand_sphere(N_e,D,rho) - Pe;
                Pe = Pe + (1+noise).*Ve;
            end        
    end
    % Copy to total particle position and velocity arrays
    if EM_swarms
        P = [Pm; Pe];
        V = [Vm; Ve];
    end
    
    %% Check bounds on position and velocity
    % Find out of bound variables
    below_LB = P < LB;
    above_UB = P > UB;
    below_Vmin = V < -Vmax;
    above_Vmax = V > Vmax;
    out_of_bounds = logical(below_LB+above_UB);
    % Bring particle's variables to within bounds
    P(below_LB) = LBmat(below_LB);
    P(above_UB) = UBmat(above_UB);
    % Set rebound velocity in those variables and limit to maximum
    V(out_of_bounds) = -reboundVel*V(out_of_bounds);
    V(below_Vmin) = -Vmaxmat(below_Vmin);
    V(above_Vmax) = Vmaxmat(above_Vmax);
    
    %% Information display
    % Setup display header
    if it == 1 
        fprintf('\n Iteration     FEs     Best aim      Best current aim   Mean aim    Stalled it. \n');
    end
    % Print display
    fprintf('   %3d      %5d      %.2e          %.2e       %.2e         %d\n', it, FE, bestGlobalAim, bestAimIter(it,1), mean(bestOwnAim), stalledIt);
    % Plot
    cmap = jet(N);
    for p=1:N
        plot(it,Aim(p),'o','Color',cmap(p,:),'MarkerSize',ms,'Parent',axes1);
        title(["Best global aim: " + num2str(bestGlobalAim,'%.4f') + ", best aim on iteration " + num2str(it) + ": " + num2str(bestAimIter(it,1),'%.4f')],'FontWeight','bold','FontSize',axes_size*0.6);
    end
    plot(1:it,bestAimIter(1:it,1),'k-','LineWidth',lw,'Parent',axes1);
    drawnow;
end

end
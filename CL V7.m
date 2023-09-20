%% Code Canals

%% Cleaning

close all
clear all
clc

%% Configuration

format long

%% Options

predefined_route_selection = 2;

plots = 1; % 1 if we want to do plots or 0 if we do not want plots
individual_plots = 0; % 1 if we want individual plots to generate the ones that will be on the paper 0 if we want joined plots
real_time_plot = 1;
video = 1; % real_time_plot must be also active to generate video

cannal_model = 2; % 1 = reduced model % 2 = ampied model

soft_hard_constraints = 1; % 1 = soft constraints % 2 = hard constraints 
%When it goes with hard constraints--> better to tighten both SR and PD
%otherwise PD does not show good results 

segments_used = 1; % 1 = only water reaches % 2 = using all segments

measuring_noise = 0; % If 1 we will consider y(k)=eye(nx)*x(k) + noise with sigma_xi 
measuring_noise_MPC = 0; % If 1 we will consider y(k)=eye(nx)*x(k) + noise with sigma_xi even in MPC (don't make sense because MPC is not doing anything to decrease sigma_xi) 

J_postcomputed = 1; % if 0 use already computed J and if 1 used post computed J

do_SR = 1;
do_PD = 1;
do_MPC = 1;

unlimited_boundaries_PD = 0; % Do predefined route version, with no limit boundaries tightening (if there is no measuring noise, this will be the same than classical MPC where 
                             % only one segment is measured at a time and in the others there is noise)  

%% Load

if cannal_model == 1
    
    load('reduced_cannal','A','C','D','Bu');

elseif cannal_model == 2
    
    load('amplied_cannal','A','C','D','Bu');
    
end

nx = size(A,2);  % number of states
nu = size(Bu,2); % number of inputs
nw = size(D,2);  % number of disturbances
ny = size(C,1);  % number of outputs

%% Parameters (must be set in the simulation)

% Simulation time

Tsim = 60;  %300 

% Sample time

sampletime = 0.5; % hours

% Prediction horizon

Np   = 10;

% Routes Prediction horizon reduction

nred = 7;

% States

if cannal_model == 1

    x   = [0.4 0.2 0.25 0.5 0.13 0.34]';  % Initial x. just a random size 
    initial_sigma_x_segment = [0 0 0 0 0 0]'; % Initial sigma of the x

elseif cannal_model == 2

    x = 0.5*[0.54 0.46 0.37 0.16 0.47 0.65 0.48 0.8 0.14 0.42 ...
             0.43 0.38 0.76 0.79 0.18 0.48 0.44 0.64 0.71 0.75 ...
             0.27 0.68 0.65 0.16 0.12 0.50 0.96 0.34 0.58 0.22]';
%     x = [0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00]';
    initial_sigma_x_segment = [0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05...
                              0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05...
                              0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05]';
    
end

% Disturbances

if cannal_model == 1

    mu_disturbances = [0 0.1]'; % Represent the mu of the disturbances. Array of dimension O x 1
    sigma_disturbances = [0.01 0.01]'; % Represent the sigma of the disturbances. Array of dimension O x 1
    
elseif cannal_model == 2
    
%     sigma_disturbances = [0.3 0.3 0.3 0.3 0.3 0.3 0.3]';
%     mu_disturbances = [zeros(1,20) -0.9*ones(1,30) zeros(1,30);
%                        zeros(1,30) -0.8*ones(1,20) zeros(1,30);
%                        zeros(1,40) -0.6*ones(1,10) zeros(1,30);
%                        zeros(1,35) -0.5*ones(1,25) zeros(1,20);
%                        zeros(1,25) +0.3*ones(1,35) zeros(1,20);
%                        zeros(1,35) +0.5*ones(1,15) zeros(1,30);
%                        zeros(1,35) +0.4*ones(1,25) zeros(1,20)];
    sigma_disturbances = [0.3 0.3 0.3 0.3 0.3 0.3 0.3]';
    mu_disturbances = [zeros(1,20) -0.6*ones(1,30) zeros(1,30);
                       zeros(1,30) -0.5*ones(1,20) zeros(1,30);
                       zeros(1,40) -0.3*ones(1,10) zeros(1,30);
                       zeros(1,35) -0.2*ones(1,25) zeros(1,20);
                       zeros(1,25) +0.4*ones(1,35) zeros(1,20);
                       zeros(1,35) +0.3*ones(1,15) zeros(1,30);
                       zeros(1,35) +0.4*ones(1,25) zeros(1,20)];

end

% Constraints (lets consider for now that all the x are fixed by the same minimum and maximum level) 

maximum_slack = 0.6;
% 
% xmin = maximum_slack*[-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
%                       -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
%                       -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]; %-0.1
% xmax = maximum_slack*[+1 +1 +1 +1 +1 +1 +1 +1 +1 +1 ...
%                       +1 +1 +1 +1 +1 +1 +1 +1 +1 +1 ...
%                       +1 +1 +1 +1 +1 +1 +1 +1 +1 +1]; %+0.1
                  
xmin = maximum_slack*[-1 -inf -1   -inf -inf -inf -1   -inf -1   -inf...
                      -1 -inf -inf -inf -inf -inf -inf -inf -inf -inf...
                      -1 -inf -inf -inf -1   -inf -inf -inf -inf -inf]; %-0.1
xmax = maximum_slack*[+1 +inf +1   +inf +inf +inf +1   +inf +1   +inf...
                      +1 +inf +inf +inf +inf +inf +inf +inf +inf +inf...
                      +1 +inf +inf +inf +1   +inf +inf +inf +inf +inf];

sat_max_constraint = +0.1;
sat_min_constraint = -0.1;

% Cost

Q   =  C'*eye(ny)*C; % 200*C'*C; % 300*eye(nx);  %100 is assigned to the weight of error--- X'[Q]X
R   =  0.2*eye(nu);   % 1 is assigned to the weight of U_gate--- U'[R]U

% Soft restriction weight

H_soft = 1e10;

% Cannal

allsegments = [2 1 6 5 4 3 8 7 10 9 20 19 18 17 16 15 14 13 12 11 24 23 22 21 30 29 28 27 26 25]; % Order of the x's in the cannal

% Cannal distances

distances_reaches = [1300 400 800 2000 1700 1600 1700]; % Distance between the watter reaches
distances_allsegments = [650 650 88.88 88.88 88.88 133.33 400 400 1000 1000 167.9 167.9 167.9 167.9 167.9 167.9 167.9 167.9 167.9 188.88 355.55 355.55 355.55 533.33 340 340 340 340 340 0];

% Robot parameters

initial_segment = 1; % Initial position of the robot. NOTE THAT IF WE ARE NOT GOING TO USE ALLSEGMENTS THE INITIAL SEGMENT SHOULD BE A WATER REACH 

robot_velocity = 12; % in m/s

energy_consumption = 0.003; % units of SOC (battery) for each unit of distance

% initial_battery = 28; % Initial battery of the robot
initial_battery = 40;
battery_rechage = 2.5; % units of SOC (battery) for each sample time without moving 
                       % first used 0.3 but it was way too low
%% Check

if size(x,1) ~= nx || size(initial_sigma_x_segment,1) ~= nx
    
    printf('Error in x size\n')
    
end

if size(mu_disturbances,1) ~= nw || size(sigma_disturbances,1) ~= nw
    
    printf('Error in disturbances size\n')
    
end

%% Automatic configuration

% Video configuration

execution_time = datestr(now);

if video

    video_name_SR = strcat('Video_SR_',execution_time,'.avi'); 
    video_name_SR = strrep(video_name_SR,':','-');
    % create the video writer with 1 fps
    writerObj_SR = VideoWriter(video_name_SR);
    % set the seconds per image
    writerObj_SR.FrameRate = 10;  
    
    video_name_PD = strcat('Video_PD_',execution_time,'.avi'); 
    video_name_PD = strrep(video_name_PD,':','-');
    % create the video writer with 1 fps
    writerObj_PD = VideoWriter(video_name_PD);
    % set the seconds per image
    writerObj_PD.FrameRate = 10; datestr(now) 

end

%% Automatic Variables

% robot battery

robot_battery_SR(1) = initial_battery;
robot_battery_PD(1) = initial_battery;

% Current segment

current_segment = initial_segment;

% Distance limit

distance_limit = robot_velocity*(3600)*sampletime; % distance limit in meters

% Create distances matrix

if segments_used == 1

    distances_used = distances_reaches;
    
elseif segments_used == 2
    
    distances_used = distances_allsegments;
    
end

distances = zeros(length(distances_used));

for i = 1:size(distances,1)
    
    for j = 1:size(distances,2)
        
        if distances(i,j) == 0
        
            distances(i,j) = 0;

                for k = i:j-1

                    distances(i,j) = distances(i,j) + distances_used(k);

                end

            distances(j,i) = distances(i,j);
        
        end
        
    end
    
end

% water reaches

water_reaches = [];
 
for i = 1:size(A,1)
     
    if A(i,i) == 1

        water_reaches = [water_reaches, i];

    end

end

if segments_used == 1

    segments = water_reaches;
    
elseif segments_used == 2
    
    segments = allsegments;
    
end

ns = length(segments); % Number of segments

% Compute a and b

[a,b] = compute_a_and_b(A,D,sigma_disturbances); % Not neccesary any longer

% Predefined Route

Predefined_Route = [];

if predefined_route_selection == 1

    while size(Predefined_Route,2) < (Tsim + Np + 1)

          Predefined_Route = [Predefined_Route, segments];
          
    end
    
elseif predefined_route_selection == 2
    
    cont = 1;
        
    while size(Predefined_Route,2) < (Tsim + Np + 1)
        
        if cont == 1

              Predefined_Route = [Predefined_Route, segments];
              
        else
            
            if mod(cont,2) == 0
                
                  Predefined_Route = [Predefined_Route, segments(ns-1:-1:1)];
                  
            elseif mod(cont,2) == 1
                
                  Predefined_Route = [Predefined_Route, segments(2:1:ns)];
                  
            end

        end
        
        cont = cont + 1;
        
    end
    
end

% Estimated noise (we are going to create the estimated noise as a vector containing the corresponding \mu of the noise)

if size(mu_disturbances,2) == 1 % Only when mu_disturbance is constant in time (if not, is computed inside the temporal loop)

    for n = 1:Np

        for d = 1:nw

            W_aux(d,n)   = mu_disturbances(d); 

        end

    end

    W = W_aux(:);

end

%% Separate things for the comparison

initial_sigma_x_segment_SR = initial_sigma_x_segment;
initial_sigma_x_segment_PD = initial_sigma_x_segment;
initial_sigma_x_segment_MPC = initial_sigma_x_segment;

sigma_xi_SR = initial_sigma_x_segment_SR;
sigma_xi_PD = initial_sigma_x_segment_PD;
sigma_xi_MPC = initial_sigma_x_segment_MPC;

[Routes,distance_step,distance_traveled,wasted_energy,energy_step] = limited_route_generator2(current_segment,segments,Np-nred,distances,initial_battery,energy_consumption,distance_limit,battery_rechage);
[Route] = complete_Routes(Routes,Np,segments);

nr = size(Route,1);

x_real_SR  = x; % SR stands for Select Route
x_real_PD  = x; % PD stands for Predefined Route
x_real_MPC = x; % Normal MPC approach

x_controller_SR  = x; % SR stands for Select Route
x_controller_PD  = x; % PD stands for Predefined Route
x_controller_MPC = x; % Normal MPC approach

%% Build matrixes (is better to do this outside the t loop because they only need to be computed once)

% Au = [eye(nu); -eye(nu)];

Au = [0.1 0	  0	  0	  0	  0	  0	  0;
      0	  1	  0	  0	  0	  0	  0	  0;
      0	  0	  1	  0	  0	  0	  0	  0;
      0	  0	  0	  1	  0	  0	  0	  0;
      0	  0	  0	  0	  1	  0	  0	  0;
      0	  0	  0	  0	  0	  1	  0	  0;
      0	  0	  0	  0	  0	  0	  1	  0;
      0	  0	  0	  0	  0	  0	  0	  1;
     -0.1 0	  0	  0	  0	  0	  0	  0;
      0	 -1	  0	  0	  0	  0	  0	  0;
      0	  0  -1	  0	  0	  0	  0	  0;
      0	  0	  0  -1	  0	  0	  0	  0;
      0	  0	  0	  0  -1	  0	  0	  0;
      0	  0	  0	  0	  0  -1	  0	  0;
      0	  0	  0	  0	  0	  0  -1	  0;
      0	  0	  0	  0	  0	  0	  0  -1];

Ax = [eye(nx); -eye(nx)];
bu = ones(2*nu,1);

bx = [];

for i = 1:nx
    
    bx = [bx; xmax(i)];
    
end

for i = 1:nx
    
    bx = [bx; -xmin(i)];
    
end
 
% Build Gx

for i=1:Np
    
    Gx((i-1)*nx+1:i*nx,:)=A^i;
    
end

% Build Gu

Gu = zeros(Np*nx,Np*nu);

for i = 1:Np
    
    for j=1:i
        
        Gu((i-1)*nx+1:(i)*nx,(j-1)*nu+1:j*nu) = (A^(i-j))*Bu;
        
    end
    
end

% Build Gw

Gw = zeros(Np*nx,Np*nw);

for i = 1:Np
    
    for j = 1:i
        
        Gw((i-1)*nx+1:(i)*nx,(j-1)*nw+1:j*nw) = A^(i-j)*D;
        
    end
    
end

% Build Q_hat

Q_hat = kron(eye(Np),Q); % The same than in the commented part

% Build R_hat

R_hat = kron(eye(Np),R);

% Build cost function H

H = Gu'*Q_hat*Gu + R_hat;

% Constraints

% u constraints

Au_hat = kron(eye(Np),Au);
bu_hat = kron(ones(Np,1),bu);

% x constraints

Ax_hat = kron(eye(Np),Ax);
bx_hat = kron(ones(Np,1),bx); 

% Transform to U constraints

AU = [Ax_hat*Gu; Au_hat];

if size(mu_disturbances,2) == 1 % Only when mu_disturbance is constant in time (if not, is computed inside the temporal loop, since W must be updated)

    bU_MPC = [bx_hat-Ax_hat*Gx*x-Ax_hat*Gw*W; bu_hat];

end

%%  Real time plot variables

upper_bound_SR = [];
lower_bound_SR = [];

for index_route = 1:nr

    predicted_upper_bound_SR{index_route} = [];
    predicted_lower_bound_SR{index_route} = [];

end

upper_bound_PD = [];
lower_bound_PD = [];

predicted_upper_bound_PD = [];
predicted_lower_bound_PD = [];

Route_travelled_SR = [];
Route_travelled_PD = [];

%% Temporal loop

for t = 1:Tsim

    fprintf('Simulation time step %d out of %d \n', t, Tsim);
    
    ST = clock;
    
    %% Update Estimated disturbance
    
    for n = 1:Np

        for d = 1:nw

            W_aux(d,n)   = mu_disturbances(d,t+n-1); 

        end

    end

    W = W_aux(:);
    
    %% F must be updated inside the loop since it depends on x. Note that W is the predicted w and in all cases we use the \mu of the disturbances to predict it
    
    if do_SR
    
%         F_SR  =  x_real_SR(:,t)'*Gx'*Q_hat*Gu + W'*Gw'*Q_hat'*Gu;
        F_SR  =  x_controller_SR(:,t)'*Gx'*Q_hat*Gu + W'*Gw'*Q_hat'*Gu;

    end

    if do_PD
        
%         F_PD  =  x_real_PD(:,t)'*Gx'*Q_hat*Gu + W'*Gw'*Q_hat'*Gu;
        F_PD  =  x_controller_PD(:,t)'*Gx'*Q_hat*Gu + W'*Gw'*Q_hat'*Gu;

    end

    if do_MPC

%         F_MPC = x_real_MPC(:,t)'*Gx'*Q_hat*Gu + W'*Gw'*Q_hat'*Gu;
        F_MPC = x_controller_MPC(:,t)'*Gx'*Q_hat*Gu + W'*Gw'*Q_hat'*Gu;

    end 
    
    %% Disturbance in t
    
    for d = 1:nw
    
        w(d,t) = normrnd(mu_disturbances(d,t),sigma_disturbances(d));
        
    end
    
    %% For the Select Best Route case

    if do_SR
        
        if t > 1 % Compute the possible routes from the new position and state
            
            [Routes,distance_step,distance_traveled,wasted_energy,energy_step] = limited_route_generator2(Visited_segment_SR(t-1),segments,Np-nred,distances,robot_battery_SR(t),energy_consumption,distance_limit,battery_rechage);
            [Route] = complete_Routes(Routes,Np,segments);
            
            Routes_in_iteration{t} = Route;

            nr = size(Route,1);
        
        end
    
        % Prepare variables for parfor
        
        for index_route = 1:nr
            
            aux1{index_route} = zeros(nx,1);
            aux2{index_route} = zeros(nx,1);
            
            bx_hat_SR{t,index_route} = zeros(size(bx_hat));
            bU_SR{t, index_route} = zeros(size(bx_hat,1)+size(bu_hat,1),1);
            
            U_SR{index_route,t} = zeros(sum(size(AU)),1); % why its 105?
            
            J_SR(index_route) = 0; 
            
            predicted_upper_bound_SR{index_route} = [];
            predicted_lower_bound_SR{index_route} = [];
            
            predicted_sigma_SR{index_route} = [];
            
        end
        
        index_route = 1;
        
        parfor index_route = 1:nr

            Route_predicted_SR = Route(index_route,:);

            % Predict sigma

            predicted_sigma_SR{index_route}(:,:,t) = predict_sigma_V2(sigma_xi_SR(:,t),Route_predicted_SR,Np,A,D,sigma_disturbances); % Function that predicts how sigma is evolving in the route            
                
            % Build the bx for all the routes
            
            bx_hat_SR{t,index_route} = []; % We are building bx_hat directly
    
            for k = 1:Np
    
                for i = 1:nx
                  
                  aux1{index_route}(i,1) = min(sat_min_constraint, xmin(i) + 3.*predicted_sigma_SR{index_route}(i,k,t));
                  aux2{index_route}(i,1) = max(sat_max_constraint, xmax(i) - 3.*predicted_sigma_SR{index_route}(i,k,t));
    
                end
                
                predicted_upper_bound_SR{index_route} = [predicted_upper_bound_SR{index_route} aux2{index_route}(:,1)];
                predicted_lower_bound_SR{index_route} = [predicted_lower_bound_SR{index_route} aux1{index_route}(:,1)];
                
                bx_hat_SR{t, index_route} = [bx_hat_SR{t, index_route}; +aux2{index_route}; -aux1{index_route}];
    
            end
    
            % Updated Agregated constraint
    
            bU_SR{t, index_route} = [bx_hat_SR{t, index_route}-Ax_hat*Gx*x-Ax_hat*Gw*W; bu_hat];
            
            % Compute MPC in every route
            
            if soft_hard_constraints == 1
                
                [U_SR{index_route,t},J_SR(index_route)] = soft_MPC(H, F_SR, AU, bU_SR{t, index_route}, H_soft);                
                
            elseif soft_hard_constraints == 2
                
                [U_SR{index_route,t},J_SR(index_route)] = hard_MPC(H, F_SR, AU, bU_SR{t, index_route});                
                
            end      
    
        end
        
         % Store J_SR
        
        J_SR_precomputed{t} = J_SR';
        
        % Predicted x for SR per route
        
        for index_route = 1:nr

            X_SR{index_route,t} = Gx*x_controller_SR(:,t) + Gu*U_SR{index_route,t}(1:Np*nu) + Gw*W;    
            
            J_SR_postcomputed{t}(index_route,1) = X_SR{index_route,t}'*Q_hat*X_SR{index_route,t} + (U_SR{index_route,t}(1:Np*nu))'*R_hat*U_SR{index_route,t}(1:Np*nu);
        
        end
          
        % Select best route
        
        if J_postcomputed == 1
            
            J_SR_used{t} = J_SR_postcomputed{t}; % IT WORKS BETTER
            
        elseif J_postcomputed == 0
            
            J_SR_used{t} = J_SR_precomputed{t}; 
            
        end
        
        [selected_J(t),selected_route(t)] = min(J_SR_used{t});
        
        nr_min_J(t) = sum(J_SR_used{t} == selected_J(t));
        
        Candidate_route{t} = (J_SR_used{t} == selected_J(t)).*Route;
        Candidate_route{t} = Candidate_route{t}(any(Candidate_route{t},2),:);
        
        route_id = (J_SR_used{t} == selected_J(t)).*linspace(1,nr,nr)';
        route_id = nonzeros(route_id);
        
        if nr_min_J(t) > 1
            
            selected_route(t) = route_id(randi(nr_min_J(t)));
            
        end
    
        % Segment visited by robot
        
        Visited_segment_SR(t) = Route(selected_route(t),1);
        
        % Robot battery 
        
        robot_battery_SR(t+1) = energy_step(selected_route(t),1);
        
        % Update bounds
        
        upper_bound_SR = [upper_bound_SR predicted_upper_bound_SR{selected_route(t)}(:,1)];
        lower_bound_SR = [lower_bound_SR predicted_lower_bound_SR{selected_route(t)}(:,1)];
    
        % Update sigma_xi
        
        for i = 1:nx
    
            sigma_xi_SR(i,t+1) = a(i)*sigma_xi_SR(i,t) + b(i);
    
        end
    
        if Visited_segment_SR(t) >0
    
            % Set sigma of visited segment to 0
            
            sigma_xi_SR(Visited_segment_SR(t),t+1) = 0;
    
        end
        
        % Predict x in the selected route
   
        X_SR_selected(:,t) = Gx*x_controller_SR(:,t) + Gu*U_SR{selected_route(t),t}(1:Np*nu) + Gw*W;

        predicted_x_SR_selected(:,:,t) = reshape(X_SR_selected(:,t),nx,Np);
       
    end
      
    %% For the Predefined Route Case

    if do_PD
    
        % Build the bx for the route
    
        bx_hat_PD{t} = []; % We are building bx_hat directly
        
        % start predicted boundaries
        
        predicted_upper_bound_PD = [];
        predicted_lower_bound_PD = [];
        
        % Visited segment
        
        Visited_segment_PD(t) = Predefined_Route(t);    
        
        % Update robot battery
        
        if t == 1
            
            robot_battery_PD(t+1) = update_battery_robot_with_route(segments,initial_segment,Visited_segment_PD(t),distances,robot_battery_PD(t),battery_rechage,energy_consumption);
            
            
        else
            
            robot_battery_PD(t+1) = update_battery_robot_with_route(segments,Visited_segment_PD(t-1),Visited_segment_PD(t),distances,robot_battery_PD(t),battery_rechage,energy_consumption);
            
        end
        
        % If not enough battery
        
        if robot_battery_PD(t+1) < 0
            
            robot_battery_PD(t+1) = 0;
            
            Predefined_Route(t+1:end) = Predefined_Route(t:end-1);
            
        end

        % Next route
        
        Route_predicted_PD = Predefined_Route(t:t+Np-1);
        
        % Predict sigma
        
        predicted_sigma_PD(:,:,t) = predict_sigma_V2(sigma_xi_PD(:,t),Route_predicted_PD,Np,A,D,sigma_disturbances); % Function that predicts how sigma is evolving in the route
        
        % start mpc problem
        
        if ~unlimited_boundaries_PD
    
            for k = 1:Np

                for i = 1:nx

                    aux1_PD(i,1) = min(sat_min_constraint, xmin(i) + 3.*predicted_sigma_PD(i,k,t));
                    aux2_PD(i,1) = max(sat_max_constraint, xmax(i) - 3.*predicted_sigma_PD(i,k,t));


                end

                % expected boundaries

                predicted_upper_bound_PD = [predicted_upper_bound_PD aux2_PD(:,1)];
                predicted_lower_bound_PD = [predicted_lower_bound_PD aux1_PD(:,1)];

                bx_hat_PD{t} = [bx_hat_PD{t}; aux2_PD; -aux1_PD];

            end

            % Update boundaries

            upper_bound_PD = [upper_bound_PD predicted_upper_bound_PD(:,1)];
            lower_bound_PD = [lower_bound_PD predicted_lower_bound_PD(:,1)];
            
            % Updated Agregated constraint
    
            bU_PD{t} = [bx_hat_PD{t}-Ax_hat*Gx*x-Ax_hat*Gw*W; bu_hat];
        
        elseif unlimited_boundaries_PD
            
            predicted_upper_bound_PD = xmax'*ones(1,Np);
            predicted_lower_bound_PD = xmin'*ones(1,Np);
            upper_bound_PD = [upper_bound_PD xmax']; 
            lower_bound_PD = [lower_bound_PD xmin'];
            
            bU_PD{t} = [bx_hat-Ax_hat*Gx*x-Ax_hat*Gw*W; bu_hat];
            
        end
    
        % Compute MPC in every 
        
        if soft_hard_constraints == 1
        
            [U_PD{t},J_PD(t)] = soft_MPC(H, F_PD, AU, bU_PD{t}, H_soft);
        
        elseif soft_hard_constraints == 2
            
            [U_PD{t},J_PD(t)] = hard_MPC(H, F_PD, AU, bU_PD{t});
            
        end
    
        % Update sigma_xi
        
        for i = 1:nx
    
            sigma_xi_PD(i,t+1) = a(i)*sigma_xi_PD(i,t) + b(i);
    
        end
    
        if Visited_segment_PD(t) > 0
    
            % Set sigma of visited segment to 0
            
            sigma_xi_PD(Visited_segment_PD(t),t+1) = 0;
    
        end
        
        % Predicted x for PD

        X_PD(:,t) = Gx*x_controller_PD(:,t) + Gu*U_PD{t}(1:Np*nu) + Gw*W;

        predicted_x_PD(:,:,t) = reshape(X_PD(:,t),nx,Np);
        
    end
    
    %% Solve MPC
    
    if size(mu_disturbances,2) > 1
        
        bU_MPC = [bx_hat-Ax_hat*Gx*x-Ax_hat*Gw*W; bu_hat];
        
    end
    
        for i = 1:nx
    
            sigma_xi_MPC(i,t+1) = a(i)*sigma_xi_MPC(i,t) + b(i);
    
        end

    if do_MPC
        
        if soft_hard_constraints == 1
    
            [U_MPC{t},J_MPC{t}] = soft_MPC(H, F_MPC, AU, bU_MPC, H_soft);
            
        elseif soft_hard_constraints == 2
            
            [U_MPC{t},J_MPC{t}] = hard_MPC(H, F_MPC, AU, bU_MPC);
            
        end
        
        % Predicted x for MPC

        X_MPC(:,t) = Gx*x_controller_MPC(:,t) + Gu*U_MPC{t}(1:Np*nu) + Gw*W;

        predicted_x_MPC(:,:,t) = reshape(X_MPC(:,t),nx,Np);

    end   
        
    %% Update x for the next iteration for SR

    if do_SR
        
        applied_u_SR(:,t) = U_SR{selected_route(t),t}(1:nu);
        
        x_real_SR(:,t+1) = A*x_real_SR(:,t) + Bu*applied_u_SR(:,t) + D*w(:,t);
        x_controller_SR(:,t+1) = A*x_controller_SR(:,t) + Bu*applied_u_SR(:,t) + D*mu_disturbances(:,t);
        
        x_controller_SR(Visited_segment_SR(t),t+1) = x_real_SR(Visited_segment_SR(t),t+1);
                
        if measuring_noise
            
            for i = 1:nx
            
                x_real_SR(i,t+1) = x_real_SR(i,t+1) + normrnd(0,sigma_xi_SR(i,t));
            
            end
                        
        end
        
        % real J
    
        realJ_SR(t) = x_real_SR(:,t)'*Q*x_real_SR(:,t) + applied_u_SR(:,t)'*R*applied_u_SR(:,t);

    end
    
    %% Update x for the next iteration for PD

    if do_PD
        
        applied_u_PD(:,t) = U_PD{t}(1:nu);
        
        x_real_PD(:,t+1) = A*x_real_PD(:,t) + Bu*applied_u_PD(:,t) + D*w(:,t);
        x_controller_PD(:,t+1) = A*x_controller_PD(:,t) + Bu*applied_u_SR(:,t) + D*mu_disturbances(:,t);
        
        x_controller_PD(Visited_segment_PD(t),t+1) = x_real_PD(Visited_segment_PD(t),t+1);
        
        if measuring_noise
            
            for i = 1:nx
            
                x_real_PD(i,t+1) = x_real_PD(i,t+1) + normrnd(0,sigma_xi_PD(i,t));
            
            end
                        
        end
        
        % real J
    
        realJ_PD(t) = x_real_PD(:,t)'*Q*x_real_PD(:,t) + applied_u_PD(:,t)'*R*applied_u_PD(:,t);

    end
    
    %% Update x for the next iteration for Classic MPC

    if do_MPC
    
        applied_u_MPC(:,t) = U_MPC{t}(1:nu);
        
        x_real_MPC(:,t+1) = A*x_real_MPC(:,t) + Bu*applied_u_MPC(:,t) + D*w(:,t);
        x_controller_MPC(:,t+1) = x_real_MPC(:,t+1);
        
        if measuring_noise_MPC
            
            for i = 1:nx
            
                x_real_MPC(i,t+1) = x_real_MPC(i,t+1) + normrnd(0,sigma_xi_MPC(i,t));
            
            end
                        
        end
        
        % real J
    
        realJ_MPC(t) = x_real_MPC(:,t)'*Q*x_real_MPC(:,t) + applied_u_MPC(:,t)'*R*applied_u_MPC(:,t);

    end
    
    %% Last time
    
    FT = clock;
    
    comp_time(t) = etime(FT,ST);
    
    last_t = t; % just in case simulation stops earlier than what it should
    
%     e(1,t)=max(abs(bU_PD{t}-bU_SR{t})); 
%     e(2,t)=max(abs(U_PD{t}-U_SR{t}));
%     e(3,t)=max(abs(x_real_PD(:,t)-x_real_SR_selected(:,t))); 
%     e(4,t)=max(abs(applied_u_SR(:,t)-applied_u_PD(:,t)));
%     e(5,t)=max(abs(x_real_SR(:,t)-x_real_PD(:,t)));
    
    %% Real time plot
    
    if real_time_plot
   
        Generate_real_time_plot
        
    end

end

if do_SR && do_PD && do_MPC

    comparison = [realJ_SR; realJ_PD; realJ_MPC];
    
    names = ["Proposed approach";"Predefined Route";"Classical MPC"];
    
    accumulated_stage_cost = sum(comparison,2);
    
    fprintf('The accumulated stage cost is:\n - Proposed approach -> %f \n - Predefined route -> %f \n - Classical MPC -> %f \n',accumulated_stage_cost(1),accumulated_stage_cost(2),accumulated_stage_cost(3))  
    
    [~,podium] = sort(accumulated_stage_cost);
    diffs = triu( accumulated_stage_cost(:).'- accumulated_stage_cost(:) );
    diffs = abs(diffs);
    biggest_difference = max( diffs(diffs>=0));
    
    fprintf('The best one is %s, the second one is %s and the worst one is %s. The biggest difference between the three methods is %f \n',names(podium(1)),names(podium(2)),names(podium(3)),biggest_difference)

end

if do_SR

    water_levels_SR = C*x_real_SR;
    
end

if do_PD
    
    water_levels_PD = C*x_real_PD;

end

if do_MPC
    
    water_levels_MPC = C*x_real_MPC;

end

water_level_mean_disturbances = D*mu_disturbances;

%% Plots

if plots

    Generate_plots
    
end

%% Close video

if video

    % close the writer object
    close(writerObj_SR);
    close(writerObj_PD);

end

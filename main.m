clc
clear all
close all
%% Purpose
% This script reproduces simulations from the article. Note that it may
% require significant time to accomplish all calculations. So the primary 
% goal of this script is to provide a source code of initial simulations.
% To see results obtained here, you could run Creating_Figures.m with 
% precalculated data. 


%% Global variables
global Results n previous_loadings thermal_memory

%% Choosing the case for simulation
% Select the number of cases (1-2) depending on cable constraints 
n = 2;

% Two possible cases:
% case 1 - Static Thermal Rating (STR) is used as a cable limit
% case 2 - Dynamic Thermal Rating (DTR) is used as a cable limit

%  In both cases, 5 market strategies are used:
% - "P50 strategy": Using P50 over the entire year (Reference for industry, “business as usual?)
% - "Variable quantile" : Changing the quantile each day. No obligation to follow the same quantile over the year
% - "Best fixed quantile": Calculating revenue for fixed quantile P01-P99. Choosing the quantile with maximum revenue.
% - "Optimal power profile":Using optimal power profile st wind installed capacity. Optimal power profile is calculated by fmincon
% - "Actual power profile": Using a measured power profile of wind farm (after application of contraints)


% The script by default is run for market prices 2018. if you would like to 
% run simulations for 2022 just uncomment lines 269-271. 

%% Launching the simulations
% Admissible current of cable
I_adm=871; % A

% Convert I_adm to Pnom (without overplanting)
P_nom=I2P(I_adm); 

% Vector of studied overplanting rates, pu
Overplanting_rate=[1 1.1 1.2 1.3 1.5 1.7 1.9 2]; % pu

% Overplanting capacities in MW
Overplanting_capacity=Overplanting_rate*P_nom;

switch n % switch the case 
    case 1 % if STR is used as cable constraints 
        tic
        for overplanting_idx=1:length(Overplanting_capacity) % for each overplanting rate 

            % Assign the installed capacity of offshore wind farm, MW
            Installed_capacity_offshore=Overplanting_capacity(overplanting_idx);

            % Do simulations
            [Results]=do_simulations(Installed_capacity_offshore);

            % Save results 
            Overplanting_results.(['Capacity_' num2str(round(Overplanting_capacity(overplanting_idx))) '_MW'])=Results;
            Results=[];

            % Save mat file
            filename=sprintf('main_simulations_STR_YEAR.mat');
            save(filename);

        end % end of for cycle 

        time=toc;
        
    case 2 % if DTR is used as cable constraints 
        tic
        for overplanting_idx=5:length(Overplanting_capacity)
            % Assign the installed capacity of offshore wind farm, MW
            Installed_capacity_offshore=Overplanting_capacity(overplanting_idx);

            % Do simulations
            [Results]=do_simulations(Installed_capacity_offshore);

            % Save results 
            Overplanting_results.(['Capacity_' num2str(round(Overplanting_capacity(overplanting_idx))) '_MW'])=Results;
            Results=[];
            previous_loadings=[];
            thermal_memory=[];

            % Save results as a mat file
            filename=sprintf('main_simulations_DTR_YEAR.mat');
            save(filename);

        end % end of for cycle 
        time=toc;

    otherwise % if n is not equal to 1 or 2
        error('Check n (the case number). It should be integer 1 or 2')

end % end of switch 



%% Embedded functioin simulations (used in script above)

function [Results]=do_simulations(Installed_capacity_offshore)
%% Description of the script
% This script solves a  operation problem when temperature or current limit is
% used as a constraint for loading of submarine cable.

% SYNTAX:
%   [Results]=do_simulations(Installed_capacity_offshore)

% INPUT:
% Installed_capacity_offshore

% OUTPUT:
% Results a structure variable 
%% Variables and constants
% Set global variable (used in any functions) but without specifying values
global previous_loadings optimization n
global thermal_memory days memoring I_adm theta_VECTOR_ISGT
global thermal_memory_preload

time_resolution=2;
% 1: 60 min
% 2: 15min

% Deciding  which constraints to use as a function of case n
if n==2% switch on temperature constraint and off the current constraint
    temperature_constraints_on=1;
    current_constraints_on=0;
else % switch off temperature constraint and on current_constraints
    temperature_constraints_on=0;
    current_constraints_on=1;
end

optimization_on=1; % permit the optimization with fmincon

load('Elia_Jan13_2016_Jan13_2017_powers.mat') % Elia 15-min data on powers
load('ENTSO_E_Jan13_2016_Jan13_2017_prices.mat') % Entso-e 30- min data on prices

% Find indexes corresponding to the studied horizon Jan13 2016-Jan13_2017
Start=t_year(1);
End=t_year(end);

days=1; % order number of studied day (needed in IEC6853_2.m!)

%% Constructing the preloadings
% Load Elia data o offshore power output and monitored capacity
load('Elia_Fev01_2012_Jan12_2017_preloads.mat')
load('Elia_Fev01_2012_Jul28_2021_monitored_capacity_offshore.mat')

% Note that the monitored_capacity is growing thus it is better convert it
% into Load factor in pu:
t_preload_start=t(102241); % Jan 1 2015 00:00

% Update preloads, t_preload, Monitored_capacity to start from Jan 1 2015 00:00
t_preload=t_preload(102241:end);
Preloads=Preloads(102241:end);
Monitored_capacity_preload=Monitored_capacity(102241:102240+length(Preloads));

% Find a preload load factor
Load_factor_preload=Preloads./Monitored_capacity_preload;

% Checking if there are NaN in the array of Load_factor_preload. This happens
% when no data is avaialble.
condition = isnan(Load_factor_preload);
if sum(condition)>0 % if there is at least one NaN in the array
    NaN_index=find(condition==1); % Find the index of NaN values
    Load_factor_preload(NaN_index)=0; % Set Load_factor to 0
end

% Datetime of very first day
Start_preload=t_preload_start; % not t_preload(1)

% Create a full horizon : preload + studied year
t_full=[t_preload; t_year];
Load_factor10_full=[Load_factor_preload; Load_factor10_year]; % P10
Load_factor50_full=[Load_factor_preload; Load_factor50_year]; % P50
Load_factor90_full=[Load_factor_preload; Load_factor90_year]; % P90
Load_factor_measur_full=[Load_factor_preload; Load_factor_measur_year];
Monitored_capacity_full=[Monitored_capacity_preload; Monitored_capacity_year];

% Find index of end date and first date
idx_Start=find(t_full==Start);
idx_End=find(t_full==End);
idx_Start_preload=find(t_full==Start_preload);

% Create an empty  power vector over the full horizon (Fev2 2012-Jan13 2017)
P_test=zeros(idx_End-idx_Start_preload+1,1);

% Set preload wind_farm_profile considering installed capacity
P_test(1:length(Load_factor_preload),1)=Load_factor_preload*339.4387;
disp('Preload at 339.4387 MW')

% Converting to current
I_test=P2I(P_test);
I_test(I_test>871)=871;
disp('I_test does not exceed 871 A')

if exist('thermal_memory_preload')
    thermal_memory_preload=[];
end

% Load preload temperature rises (after current 871 A was applied)
load('theta_VECTOR_preload_339MW_871A.mat')

% This theta_VECTOR_339MW_871A gives us the precalculated temperature rises 
% to save the time. However, if neccesary this calculation can be done 
% explicitly by uncommenting the code below. 

% if exist('theta_VECTOR_preload') && ~(length(theta_VECTOR_preload)==0)
%     % do nothing
% else % calculate theta_VECTOR (inside of cable_thermal_model_IEC_60853_2
%     [Tmax,Temperature,~]=cable_thermal_model_IEC_60853_2(I_test);
% end

% Create a thermal memory the full horizon (Jan15 2015-Jan13 2017)
thermal_memory=theta_VECTOR_preload;

% Delete a part of thermal memory during preload period
thermal_memory(1:length(Load_factor_preload),:)=[];

% Save a thermal memory_preload at the studied horizon (Jan13 2016-Jan11 2017)
thermal_memory_preload=thermal_memory;

% Setting the admissible current I_adm
I_adm=871;

% Convert to admissible power 
P_adm=I2P(I_adm);

close all

%% Load data on prices

% Day-ahaed and imbalance prices from ENTSO-E for 2015-2020
load('data_DA_IMB_years_2015 2020.mat')
load('DA_2022.mat')
load('Cb_2022.mat')

% create a figure
figure('DefaultAxesFontSize',14)
subplot(3,1,1)
hold on
plot(time_2018,DA_2018) % plot day-ahead prices for 2018
ylabel('Day-ahead price, Euro')
title('Day-ahead price')

subplot(3,1,2)
hold on
plot(time_2018,Cb_plus_2018) % plot imbalance+ prices for 2018
ylabel('Imbalance+ price, Euro')
title('Imbalance+ price')

subplot(3,1,3)
hold on
plot(time_2018,Cb_minus_2018) % plot imbalance- prices for 2018
ylabel('Imbalance- price, Euro')
title('Imbalance- price')

% Create a zero vector (considering the horizon of preloads)
Price_preload=zeros(length(Load_factor_preload),1);

% Create the array of prices at t_full (2018 year)
DA_2018=[Price_preload;DA_2018];
Cb_plus_2018=[Price_preload;Cb_plus_2018];
Cb_minus_2018=[Price_preload;Cb_minus_2018];


% Uncomment this section if simulation 2022 are needed 
% DA_2018=[Price_preload;DA_2022];
% Cb_plus_2018=[Price_preload;Cb_plus_2022];
% Cb_minus_2018=[Price_preload;Cb_minus_2022];

%% Perofrming day-by-day optimization
days=1; % set a first day (13 Jan 2016)

% While cycle: Solving the day-by-day problem
while ~(idx_Start+96*(days-1)+95>idx_Start+96*(365-1)+95) % till the last day (Jan13 2017) is reached
    %%  Preparing the data
    % Extract the indexes for beginning and the end of the day
    Day_start=idx_Start+96*(days-1);
    Day_end=Day_start+95;
    
    % Extrating the forecasted power profiles for J day
    P50_day=Load_factor50_full(Day_start:Day_end)*Installed_capacity_offshore;
    P10_day=Load_factor10_full(Day_start:Day_end)*Installed_capacity_offshore; % P10
    P90_day=Load_factor90_full(Day_start:Day_end)*Installed_capacity_offshore; % P90
    
    % Extrating the measured power profiles for J day
    P_measur_day=Load_factor_measur_full(Day_start:Day_end)*Installed_capacity_offshore;
    
    % Extrating the forecasted prices for J day in 2018
    DA_day_2018=DA_2018(Day_start:Day_end);
    Cb_plus_day_2018=Cb_plus_2018(Day_start:Day_end);
    Cb_minus_day_2018=Cb_minus_2018(Day_start:Day_end);
    
    % Checking if there are NaN in the forecasts. If so, we
    % assume that NaN == 0.
    condition1 = isnan(P10_day);
    condition2 = isnan(P50_day);
    condition3 = isnan(P90_day);

    if sum(condition1)>0 % if there is at least one NaN in P10
        NaN_index=find(condition1==1); % Find the index of NaN values
        P10_day(NaN_index)=0; % Set power production to 0 also known
    elseif sum(condition2)>0% if there is at least one NaN in P50
        NaN_index=find(condition2==1); % Find the index of NaN values
        P50_day(NaN_index)=0; % Set power production to 0 also known
    elseif sum(condition3)>0 % if there is at least one NaN in P90
        NaN_index=find(condition3==1); % Find the index of NaN values
        P90_day(NaN_index)=0; % Set power production to 0 also known
    else %there is no NaN in array
        % do nothing
    end
    
    % ------------------------Quantile start-------------------------------
    % load data of fitting for given day:
    load('Fitted_data.mat') % for P_cdf
    Fitting_result = Fitted_data.(sprintf('days_%d', days)).Fitting_result;

    % The fitting was additonally realized for Elia data to obtain 99
    % quantiles (instead of only 3 quantile initially given by Elia P90, P50
    % P10)    
    
    % Create an empty variable
    P_cdf=[];

    % Reconstruct power profiles for given quantiles:
    for i=1:96 % for each instant of time (96 values at 15-min resolution)

        % Extract cdf for given time step
        P_cdf_interm=Fitting_result.cdf_opt{i}; % for 712 MW (Elia data)
        
        % Note that fitting was done for 712.2 MW (Elia capacity) that is why we need to
        % convert data for given Installed_capacity_offshore
        P_cdf_interm=P_cdf_interm/712.2*Installed_capacity_offshore; % inMW for our overplanting rate
        
        % Extract quantiles from 1 to 99
        P_cdf_interm=P_cdf_interm(2:100);
        
        % Combine the P_cdf_interm with P_cdf
        P_cdf=cat(1,P_cdf,P_cdf_interm);
        
    end
    
    % Check if there are non-physical values (negative power or power
    % greater than installed capacity)
    P_cdf(P_cdf>Installed_capacity_offshore)=Installed_capacity_offshore;
    P_cdf(P_cdf<0)=0;
    


    % plot power profiles
    subplot (2,1,1)
    hold on
    plot(t_full(Day_start:Day_end),[P10_day,P50_day,P90_day],'linewidth',2);
    Installed_capacity_offshore_vector=linspace(Installed_capacity_offshore,Installed_capacity_offshore,96)';
    plot(t_full(Day_start:Day_end),Installed_capacity_offshore_vector,'linewidth',1);
    plot(t_full(Day_start:Day_end),P_measur_day,'linewidth',3);
    ylabel('Power, MW')
    
    % plot prices
    subplot (2,1,2)
    hold on
    plot(t_full(Day_start:Day_end),[DA_day_2018,Cb_plus_day_2018,Cb_minus_day_2018],'linewidth',0.5);
    ylabel('Prices, Euros')
    
    %% ---------------------Current constraints start----------------------
    if current_constraints_on==1
        % Set actual power profile as measured power profile
        Pactual=P_measur_day;
        
        % If actual power is greater than admissible power of cable
        if max(Pactual)>P_adm
            % Coonvert power to current
            I_actual=P2I(Pactual);
            % Find the index where actual current exceeeds the admissible
            % current
            idx=find(I_actual>I_adm);
            
            % Make actual current <= I adm
            I_actual(idx)=I_adm;
            
            % Convert back to actual power
            Pactual=I2P(I_actual);
        else % max(Pactual)=>P_adm
            % do nothing
        end

        % Revenue and Pplan for all quantiles
        for i=1:99
            if i<10
                optim_sol(i,1) = join(compose('P0%d', i));
            else
                optim_sol(i,1) = join(compose('P%d', i));
            end
        end
        for i=1:length(P_cdf)
            Pplan=P_cdf(:,i);
            Revenue_cdf_2018(i)=sum(Pplan.*DA_day_2018-Cb_minus_day_2018.*max(0,Pplan-Pactual)+Cb_plus_day_2018.*max(0,Pactual-Pplan));
        end        
        
        % Revenue and Pplan for P50
        Pplan_P50_cdf=P_cdf(:,50);
        Pplan_P50_Elia=P50_day;
        Revenue_P50_Elia_2018=sum(Pplan_P50_Elia.*DA_day_2018-Cb_minus_day_2018.*max(0,Pplan_P50_Elia-Pactual)+Cb_plus_day_2018.*max(0,Pactual-Pplan_P50_Elia));
        Revenue_P50_cdf_2018=sum(Pplan_P50_cdf.*DA_day_2018-Cb_minus_day_2018.*max(0,Pplan_P50_cdf-Pactual)+Cb_plus_day_2018.*max(0,Pactual-Pplan_P50_cdf));
        
        % Objective function (maximization of revenue) for optimal power profile
        f2018 = @(x)-(sum(x.*DA_day_2018-Cb_minus_day_2018.*max(0,x-Pactual)+Cb_plus_day_2018.*max(0,Pactual-x)));
        fun2018 = @(x) f2018(x(:)); % convert any x to a column (needed for ga otherwise we get the error)
        
        %  Solve the optimization problem with fmincon
        if optimization_on==1
            lb=zeros(96,1); % lower bound for x
            ub=Installed_capacity_offshore_vector; % upped bound for x
            %     options = optimoptions('ga','MaxGenerations',1000000,'PlotFcn',@gaplotbestf);
            %     tic
            %     [P_plan_opt,Revenue_opt,exitflag] = ga(fun,96,[],[],[],[],lb,ub,[],options);
            %     time=toc;
            options = optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'Algorithm','sqp','PlotFcns',{@optimplotfval,@optimplotfvalconstr,@optimplotfval,@optimplotx,@optimplotconstrviolation});
            tic
            [P_plan_opt,Revenue_opt,exitflag] = fmincon(fun2018,ub,[],[],[],[],lb,ub,[],options);
            time=toc;
        end
        close all
        %% Saving results for the given day for current constraints
        % Save Installed capacity
        Results.CurrentConstraints.Installed_Capacity(days,1)=Installed_capacity_offshore;
        % Save case number
        Results.CurrentConstraints.case(days,1)=n;
        
        % Revenue for the best guess power profile(Pplan==Pactual)
        Revenue_best_guess_2018=fun2018(Pactual);

        % Find the best quantile for the given day
        idx_max_2018=find(Revenue_cdf_2018==max(Revenue_cdf_2018));


        if optimization_on==1
            % Save the results for given day
            Results.CurrentConstraints.P_plan_opt{days,1}=P_plan_opt;
            Results.CurrentConstraints.Revenue_opt{days,1}=-Revenue_opt;
            Results.CurrentConstraints.exitflag{days,1}=exitflag;
        end

        % Save the revenue for each quantile 
        Results.CurrentConstraints.Revenue_cdf_2018{days,1}=Revenue_cdf_2018;

        % Save the maximal revenue for the given day
        Results.CurrentConstraints.max_Revenue_cdf_2018(days,1)=max(Revenue_cdf_2018);

        % Save the optimal variable quantile for the given day
        Results.CurrentConstraints.Quantile2018{days,1}=idx_max_2018;

        % Save revenue for Pactual
        Results.CurrentConstraints.Revenue_best_guess_2018(days,1)=-Revenue_best_guess_2018;


        % Save Revenue for P50
        Results.CurrentConstraints.Revenue_P50_Elia_2018(days,1)=Revenue_P50_Elia_2018;
        Results.CurrentConstraints.Revenue_P50_cdf_2018(days,1)=Revenue_P50_cdf_2018;


        % Save the power profiles for each quantile
        Results.CurrentConstraints.P_cdf{days,1}=P_cdf;


        % Save power profile for P50 strategy
        Results.CurrentConstraints.Pplan_P50_Elia{days,1}=Pplan_P50_Elia;
        Results.CurrentConstraints.Pplan_P50_cdf{days,1}=Pplan_P50_cdf;

        % Save power profile Pactual for the given day(after
        % application of CurrentConstraints)
        Results.CurrentConstraints.Pactual{days,1}=Pactual;
        
        % Save power profile P_measur_day for the given day (before
        % application of CurrentConstraints)
        Results.CurrentConstraints.P_measur_day{days,1}=P_measur_day;
        
        % Save the day ahead prices for given day
        Results.CurrentConstraints.DA_day_2018{days,1}=DA_day_2018;
        
        % Save the imbalance+ prices for given day
        Results.CurrentConstraints.Cb_plus_day_2018{days,1}=Cb_plus_day_2018;
        
        % Save the imbalance- prices for given day
        Results.CurrentConstraints.Cb_minus_day_2018{days,1}=Cb_minus_day_2018;

        
    else % if current_constraints_on==0
        % do nothing
    end % end of if current_constraints_on==1
    %----------------------Current constraints end-------------------------
    
    
    %% ---------------------Temperature constraints start------------------
    if temperature_constraints_on==1 % if temperature constraint is on
        % Set actual power profile as measured power profile
        Pactual=P_measur_day;
        
        % If actual power is greater than admissible power of cable
        if max(Pactual)>P_adm
            % Convert from power to current
            I_actual=P2I(Pactual);
            
            % Consider the previous loadings (needed for Temperature)
            if exist('previous_loadings')&& ~(length(previous_loadings)==0)
                I_previous_loadings=P2I(previous_loadings);
                
                I_input=[I_previous_loadings;I_actual];
            else % no previous loadings
                I_input=I_actual;
            end
            % Using the popwer shedding
            [I_output]=power_curtailement(I_input);
            
            % Convert back to power
            Pactual=I2P(I_output);
            
            % Save the Pactual for given day
            Pactual=Pactual(1+96*(days-1):1+96*(days-1)+95,1);
        else % max(Pactual)=<P_adm
            % do nothing
        end
       
        % Revenue and Pplan for all quantiles
        for i=1:99
            if i<10
                optim_sol(i,1) = join(compose('P0%d', i));
            else
                optim_sol(i,1) = join(compose('P%d', i));
            end
        end
            

        % Find the Revenue for all quantiles in cdf
        for i=1:length(P_cdf)
            Pplan=P_cdf(:,i);
            Revenue_cdf_2018(i)=sum(Pplan.*DA_day_2018-Cb_minus_day_2018.*max(0,Pplan-Pactual)+Cb_plus_day_2018.*max(0,Pactual-Pplan));
        end

        % Revenue and Pplan for P50
        Pplan_P50_cdf=P_cdf(:,50);
        Pplan_P50_Elia=P50_day;
        Revenue_P50_Elia_2018=sum(Pplan_P50_Elia.*DA_day_2018-Cb_minus_day_2018.*max(0,Pplan_P50_Elia-Pactual)+Cb_plus_day_2018.*max(0,Pactual-Pplan_P50_Elia));
        Revenue_P50_cdf_2018=sum(Pplan_P50_cdf.*DA_day_2018-Cb_minus_day_2018.*max(0,Pplan_P50_cdf-Pactual)+Cb_plus_day_2018.*max(0,Pactual-Pplan_P50_cdf));
        
        % Revenue for optimal power profile
        f2018 = @(x)-(sum(x.*DA_day_2018-Cb_minus_day_2018.*max(0,x-Pactual)+Cb_plus_day_2018.*max(0,Pactual-x)));
        fun2018 = @(x) f2018(x(:)); % convert any x to a column (needed for ga otherwise we get the error)
        
        %  Solve the optimization problem with fmincon
        if optimization_on==1
            lb=zeros(96,1); % lower bound for x
            ub=Installed_capacity_offshore_vector; % upped bound for x
            %     options = optimoptions('ga','MaxGenerations',1000000,'PlotFcn',@gaplotbestf);
            %     tic
            %     [P_plan_opt,Revenue_opt,exitflag] = ga(fun,96,[],[],[],[],lb,ub,[],options);
            %     time=toc;
            options = optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'Algorithm','sqp','PlotFcns',{@optimplotfval,@optimplotfvalconstr,@optimplotfval,@optimplotx,@optimplotconstrviolation});
            tic
            [P_plan_opt,Revenue_opt,exitflag] = fmincon(fun2018,ub,[],[],[],[],lb,ub,[],options);
            time=toc;
            close all
        end

        % Revenue for the best guess power profile(Pplan==Pactual)
        Revenue_best_guess_2018=fun2018(Pactual);

        % Find the best quantile
        idx_max_2018=find(Revenue_cdf_2018==max(Revenue_cdf_2018));

        %% Saving the results for temperature constraints
        Results.TempConstraints.Installed_Capacity(days,1)=Installed_capacity_offshore;
        % Save case number
        Results.CurrentConstraints.case(days,1)=n;
        
        % Save the results for given day
        if optimization_on==1
            Results.TempConstraints.P_plan_opt{days,1}=P_plan_opt;
            Results.TempConstraints.Revenue_opt{days,1}=-Revenue_opt;
            Results.TempConstraints.exitflag{days,1}=exitflag;
        end
        
        % Save the revenue for each quantile
        Results.TempConstraints.Revenue_cdf_2018{days,1}=Revenue_cdf_2018;

        % Save the maximal revenue at that day among fixed quantile
        Results.TempConstraints.max_Revenue_cdf_2018(days,1)=max(Revenue_cdf_2018);

        % Save optimal fixed quantile for the given day
        Results.TempConstraints.Quantile2018{days,1}=idx_max_2018;

        % Save revenue for Pactual
        Results.TempConstraints.Revenue_best_guess_2018(days,1)=-Revenue_best_guess_2018;

        % Save the revenue for P50 market strategy (based on cdf data)
        Results.TempConstraints.Revenue_P50_cdf_2018(days,1)=Revenue_P50_cdf_2018;
        
        % Save the revenue for P50 market strategy (based on Elia data)
        Results.TempConstraints.Revenue_P50_Elia_2018(days,1)=Revenue_P50_Elia_2018;
        
        % Save power profile for P50
        Results.TempConstraints.Pplan_P50_Elia{days,1}=Pplan_P50_Elia;
        Results.TempConstraints.Pplan_P50_cdf{days,1}=Pplan_P50_cdf;


        % Save the actual power profile for the given day (after
        % application of temperature constraints)
        Results.TempConstraints.Pactual{days,1}=Pactual;
        
        % Save the measured power profile for the given day (before
        % application of temperature constraints)
        Results.TempConstraints.P_measur_day{days,1}=P_measur_day;
        
        % Save the day ahead prices for given day
        Results.TempConstraints.DA_day_2018{days,1}=DA_day_2018;
        
        % Save the imbalance+ prices for given day
        Results.TempConstraints.Cb_plus_day_2018{days,1}=Cb_plus_day_2018;
        
        % Save the imbalance- prices for given day
        Results.TempConstraints.Cb_minus_day_2018{days,1}=Cb_minus_day_2018;
        
        % Save previous optimal loading of cable (needed for IEC 853-2)
        previous_loadings(1+96*(days-1):1+96*(days-1)+95,1)=Pactual;
        
        % Update the thermal memory over the whole horizon (20 years)
        memoring=1; % switch the memoring status in IEC60853_2.m
        % Create zero vector wih the length = studied horizon
        P_test=zeros(idx_End-idx_Start+1,1);
        
        % Set previous_loadings to given interval
        P_test(1:96*days,1)=previous_loadings;
        
        % Converting to current
        I_test=P2I(P_test);
        
        % Finding Tmax Temeprature
        [Tmax,Temperature,~]=cable_thermal_model_IEC_60853_2(I_test);
        
        difference=Tmax-90;
        if Tmax>90
            disp(['Previous loadings break 90 limit for : ', num2str(difference)])
        end
        memoring=[]; % switch the memoring status off in IEC60853_2.m to
        % ---------------------Temperature constraints end-----------------
    else % if temperature constraint is off
        % do nothing
    end
    
    % Go to the next day
    days=days+1
    
    close all
end % end of while cycle

% Show the installed capacity
Installed_capacity_offshore
end


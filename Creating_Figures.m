clc
clear all
close all

%% Purpose
%  This script generates figures used in the journal article [1]

% [1] I. Daminov et al, "Economic performance of overplanted offshore wind
% farm under several commitment strategies and dynamic thermal ratings of
% submarine export cable" in Applied Energy, 2023
%% Figure 1  Illustration of the electrical infrastructure of OWF.
% In our study, only the submarine section of export cable is considered

% This figure was ploted in powerpoint without using MATLAB 
%% Figure 2 Cross-section and 3D view of a cable similar to the one used in
%% our paper. Reproduced with the permission of Cableizer

% This figure was ploted without using MATLAB. See the cableizer site:
% https://www.cableizer.com/

%% Figure 3 OWF forecast performance for quantiles: P10, P50 and P90.
% Clear workspace
clear all

% Load data
load('main_simulations_DTR_2018.mat')

% Create a empty vector for each day of the year.
E=zeros(365,1);

%---------------------------- P50 forecast ------------------------------
Fs=[];% forecast vector
As=[];% actual vector

for day=1:365 % for each day

    % Extract the forecasted power profile of OWF (P50 quantile)
    F=Overplanting_results.Capacity_339_MW.TempConstraints.Pplan_P50_Elia{day};
    Fs=[Fs;F]; % add forecasted power profile into annual vector

    % Extract the daily actual power profile of OWF
    A=Overplanting_results.Capacity_339_MW.TempConstraints.P_measur_day{day};
    As=[As;A]; % add actual power profile into annual vector

    % Find the difference between daily actual and forected power profiles
    diff=F-A;

    % Calculate the mean squared error
    MSE=sum((diff).^2)/length(diff);

    % Save RMSE
    E(day,1)=sqrt(MSE);
    % E(day,1) = rmse( F , A ) ;
end

% Find the difference between annual forecated and actual power profiles
diff=Fs-As;

% Find the mean absolute error
MAE_P50=sum(abs(diff))/length(diff);

% Find the mean squared error
MSE_P50=sum((diff).^2)/length(diff);

% Find the root-mean-squared error
RMSE_P50=sqrt(MSE_P50);

% Find the time when a forecast is higher than the OWF actual output
Fs_above_As=find(Fs>As);
Fs_above_As=length(Fs_above_As)/length(Fs)*100;

% Find % of time when the forecast is equal to the actual output
Fs_equal_As=find(Fs==As);
Fs_equal_As=length(Fs_equal_As)/length(Fs)*100;

% Find the time when a forecast is lower than the OWF actual output
Fs_below_As=find(Fs<As);
Fs_below_As=length(Fs_below_As)/length(Fs)*100;

% Show the results (above, equal and below)
F_relative_A_P50=[Fs_above_As Fs_equal_As Fs_below_As];

% Plot the figure for P50 quantile
figure('DefaultAxesFontSize',14,'InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized')
subplot(1,3,2) % create a subplot
scatter(Fs,As) % plot the scatter plot
xlabel('Forecasted output,MW')
ylabel('Actual output ,MW')
xlim([0,350])
ylim([0,350])
hold on
plot(1:max(As),1:max(As),'linewidth',2) % plot the diaginal line
title('P50 quantile')

%---------------------------- P90 forecast ------------------------------
Fs=[];% forecast vector
As=[];% actual vector

for day=1:365 % for each day

    % Extract the forecasted power profile of OWF
    F=Overplanting_results.Capacity_339_MW.TempConstraints.P_cdf{day, 1};
    F=F(:,90); % P90 quantile
    Fs=[Fs;F]; % add the forecasted power profile into annual vector

    % Extract the daily actual power profile of OWF
    A=Overplanting_results.Capacity_339_MW.TempConstraints.P_measur_day{day};
    As=[As;A]; % add the actual power profile into annual vector

    % Find the difference between daily forected and actual power profiles
    diff=F-A;

    % Calculate the mean squared error
    MSE=sum((diff).^2)/length(diff);

    % Save RMSE
    E(day,1)=sqrt(MSE);
    % E(day,1) = rmse( F , A ) ;
end

% Find the difference between annual forecated and actual power profiles
diff=Fs-As;

% Find the mean absolute error
MAE_P90=sum(abs(diff))/length(diff);

% Find the mean squared error
MSE_P90=sum((diff).^2)/length(diff);

% Find the root-mean-squared error
RMSE_P90=sqrt(MSE_P90);

% Find the time when a forecast is higher than the OWF actual output
Fs_above_As=find(Fs>As);
Fs_above_As=length(Fs_above_As)/length(Fs)*100;

% Find % of time when the forecast is equal to the actual output
Fs_equal_As=find(Fs==As);
Fs_equal_As=length(Fs_equal_As)/length(Fs)*100;

% Find the time when a forecast is lower than the OWF actual output
Fs_below_As=find(Fs<As);
Fs_below_As=length(Fs_below_As)/length(Fs)*100;

% Show the results (above, equal and below)
F_relative_A_P90=[Fs_above_As Fs_equal_As Fs_below_As];

% plot the figure for P90 quantile
subplot(1,3,3) % create a subplot
scatter(Fs,As) % the scatter plot
xlabel('Forecasted output,MW')
ylabel('Actual output ,MW')
xlim([0,350])
ylim([0,350])
hold on
plot(1:max(As),1:max(As),'linewidth',2) % diagonal line
title('P90 quantile')

%---------------------------- P10 forecast ------------------------------
Fs=[];% forecast vector
As=[];% actual vector

for day=1:365 % for each day

    % Extract the daily forecasted power profile of OWF
    F=Overplanting_results.Capacity_339_MW.TempConstraints.P_cdf{day, 1};
    F=F(:,10); % P10 quantile
    Fs=[Fs;F]; % add the forecasted power profile into annual vector

    % Extract the daily actual power profile of OWF
    A=Overplanting_results.Capacity_339_MW.TempConstraints.P_measur_day{day};
    As=[As;A]; % add the actual power profile into annual vector

    % Find the difference between daily forected and actual power profiles
    diff=F-A;

    % Calculate the mean squared error
    MSE=sum((diff).^2)/length(diff);

    % Save RMSE
    E(day,1)=sqrt(MSE);

    %     E(day,1) = rmse( F , A ) ;
end

% Find the difference between annual forecated and actual power profiles
diff=Fs-As;

% Find the mean absolute error
MAE_P10=sum(abs(diff))/length(diff);

% Find the mean squared error
MSE_P10=sum((diff).^2)/length(diff);

% Find the root-mean-squared error
RMSE_P10=sqrt(MSE_P10);

% Find the time when a forecast is higher than the OWF actual output
Fs_above_As=find(Fs>As);
Fs_above_As=length(Fs_above_As)/length(Fs)*100;

% Find % of time when the forecast is equal to the actual output
Fs_equal_As=find(Fs==As);
Fs_equal_As=length(Fs_equal_As)/length(Fs)*100;

% Find the time when a forecast is lower than the OWF actual output
Fs_below_As=find(Fs<As);
Fs_below_As=length(Fs_below_As)/length(Fs)*100;

% Show the results (above, equal and below)
F_relative_A_P10=[Fs_above_As Fs_equal_As Fs_below_As];

% plot the figure for P10 quantile
subplot(1,3,1) % create a subplot
scatter(Fs,As) % the scatter plot
xlabel('Forecasted output,MW')
ylabel('Actual output ,MW')
xlim([0,350])
ylim([0,350])
title('P10 quantile')
hold on
plot(1:max(As),1:max(As),'linewidth',2) % diagonal line

%% Figure 4 Power generation forecasts and actual measurements of the OWFs
% connected to Elia’s network

% Clear workspace
clear all

% Load the power profiles from Elia
load('Elia_Jan13_2016_Jan13_2017_powers.mat')

% Create figure
figure('InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized');

% Plot load factors of offshore wind farm(s)
plot(t_year,[Load_factor90_year,Load_factor50_year,Load_factor10_year,Load_factor_measur_year])

% Plot the legend
legend('P90','P50','P10','Measured')

% Display the ylabel
ylabel('Load Factor, pu')

% The Figure reprsents annual load factors o offshore wind farms (both
% forecast and measurements. To see the Figure 4 as it is done in paper, just
% zoom in the power profile at the needed interval

%% Figure 5 Market prices in France and different situations where the DA
%% price may be greater, inbetween or less than the imbalance prices

% Clear workspace
clear all

% Load day-ahead and imbalance prices
load('data_DA_IMB_years_2015 2020.mat')

% Create figure
figure('InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized');

% Plot day-ahead price and imbalance prices+-
plot(time_2018,[DA_2018,Cb_plus_2018,Cb_minus_2018])

% Plot the legend
legend('Day-ahead','Imbalance+','Imbalance-')

% Display the ylabel
ylabel('Price, EUR/MWh')

%% Figure 6 Equivalent thermal circuit with two cells and corresponding values

% This figure was ploted without using MATLAB in powerpoint


%% Figure 7 Validation of the thermal model against reference temperature profiles
% Clear workspace
clear all

% Load data
load('Reference_temperature.mat')
load('Input_data_for_thermal_test.mat')
pretest_days=365;


% Converth I_test from 1-hour into 15-min resolution
time_init =  InputData.Time;
time_interp = (time_init(1):minutes(15):time_init(end)+minutes(45))';  % time resolution 15 min
I_test = interp1(time_init, InputData.CurrentA, time_interp,'previous','extrap');
I_pretest=linspace(609,609,24*pretest_days*4)';
I_test=[I_pretest;I_test]; % vector of current load

% Uncomment these lines if you would like to do a full calculation
% tic
% % [T_max,T,I_adm]=cable_thermal_model_IEC_60853_2(I_test);
% time=toc

% Otherwise, the load the precalculated values (to avoid long processing
% time)
load('precalculated_Tmax_T_Iadm.mat')

% Create a datetime vector for T
t_start=InputData.Time(1)-days(pretest_days);
t_end=InputData.Time(end);
t_15min=[t_start:minutes(15):t_end]';
t_15min(end+1)= t_15min(end)+minutes(15);
t_15min(end+1)= t_15min(end)+minutes(15);
t_15min(end+1)= t_15min(end)+minutes(15);
t_15min(end+1)= t_15min(end)+minutes(15);

% Convert data into 6 min resolution (as reference temperature is in this
% resolution)
dt = minutes(6);
T=timetable(t_15min,T);
T_6min = retime(T,'regular','linear','TimeStep',dt);
T_6min = timetable2table(T_6min);
T_6min=table2array(T_6min(:,2));

I_test(2:end+1)=I_test;
I=timetable(t_15min,I_test);
I_6min = retime(I,'regular','previous','TimeStep',dt);
I_6min = timetable2table(I_6min);
I_6min=table2array(I_6min(:,2));
I_6min(240*pretest_days+1:end-10)=I_6min(240*pretest_days+11:end);

% Ploting the figure
figure('DefaultAxesFontSize',14,'InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized')
subplot(3,1,1)
plot(t,I_6min(240*pretest_days+1:end),'linewidth',2)
ylabel('Current, A')

subplot(3,1,2)
hold on
plot(t,T_6min(240*pretest_days+1:end),'linewidth',2)
plot(t,Reference_temperature.Conductor,'linewidth',2)
legend('MATLAB','Reference')
ylabel('Temperature, degC')

subplot(3,1,3)
Difference=T_6min(240*pretest_days+1:end)-Reference_temperature.Conductor;
Mean_Difference=mean(abs(Difference));

plot(t,Difference,'linewidth',2)
ylabel('Difference, degC')

% About the difference between this figure and article
disp(['Note that figures, obtained here, and the figure in article are different. ' ...
    'This is because the thermal model of cable, used in the article, was ' ...
    'additionally calibrated with confidential data. This open-acces model ' ...
    'does not include these data but still demonstrate an appropriate performance'])


%% Figure 8 Block scheme for performing the simulations

% This figure was ploted without using MATLAB in powerpoint

%% Figure 9 CAPEX and OPEX of OWF as a function of the overplanting rate
% Clear workspace
clear all

%  Set of  OWF capacities
Overplanting_capacity=[339.4387  373.3825  407.3264  441.2703  509.1580...
    577.0457  644.9334  678.8773]; % MW

% Share of cable export costs in all cable installation costs(for all cables:inner and export):
coef_export=0.5;

% GPB/EUR 1.1405 in 2019
exchange_rate_GBP_EUR=1.1405;

% Calculate the CAPEX and OPEX (€) for each overplanting rate
[CAPEX,OPEX]=CATAPULT_costs_2019(coef_export,exchange_rate_GBP_EUR);

% Ploting the figure
figure('DefaultAxesFontSize',14,'InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized')

yyaxis left
plot(Overplanting_capacity/Overplanting_capacity(1),CAPEX/1e6,'linewidth',2)
ylabel('Capital expenditures, M€')
xlabel('Overplanting rate,pu')

yyaxis right
plot(Overplanting_capacity/Overplanting_capacity(1),OPEX/1e6,'linewidth',2)
ylabel('Operational expenditures, M€')

%% Figure 10 - 17
% Clear workspace
clear all

% Cases:
% Capital letters represent component which differs from Reference case (case 1)

% 1 - Reference case: P50 market strategy; Current constraint; no overplanting (338 MW);

% 2 - P50 market strategy; Current constraint; OVERPLANTING (338-712 MW);

% 3 - P50 market strategy; TEMPERATURE CONSTRAINTS; OVERPLANTING (338-712 MW);

% 4 - VARIABLE MARKET STRATEGY; Current constraint; no overplanting (338 MW);

% VARIABLE MARKET STRATEGY:
% - "P50 strategy": Using P50 over the entire year (Reference for industry, “business as usual?)
% - "Variable quantile" : Changing the quantile each day. No obligation to follow the same quantile over the year
% - "Best fixed quantile": Calculating revenue for fixed quantile P01-P99. Choosing the quantile with maximum revenue.
% - "Optimal power profile":Using optimal power profile st wind installed capacity. Optimal power profile is calculated by fmincon
% - "Actual power profile": Using a measured power profile of wind farm (after application of contraints)

% 5 - VARIABLE MARKET STRATEGY; Current constraint; OVERPLANTING (338-712 MW);

% 6 - VARIABLE MARKET STRATEGY; TEMPERATURE CONSTRAINTS; OVERPLANTING (338-712 MW);

% Extracting the relevant data
for cases=5:6
    if cases==5
        load('main_simulations_STR_2018.mat')
        for capacity_idx=1:length(Overplanting_capacity)

            eval(['P50_Revenues_case5=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Revenue_P50_cdf_2018/4']); % /4 for MWh
            eval(['Annual_P50_Revenue_case5_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(P50_Revenues_case5)']);

            eval(['Variable_quantile_Revenues_case5=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.max_Revenue_cdf_2018/4']); % /4 for MWh
            eval(['Annual_Variable_quantile_Revenues_case5_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(Variable_quantile_Revenues_case5)']);

            eval(['Optimal_fixed_quantile_Revenues_case5=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Revenue_cdf_2018  ']);
            Annual_Optimal_fixed_quantile_Revenues_case5=zeros(1,99);
            for i=1:365
                Interm_revenue=Optimal_fixed_quantile_Revenues_case5{i, 1}/4; % /4 for MWh
                Annual_Optimal_fixed_quantile_Revenues_case5=Annual_Optimal_fixed_quantile_Revenues_case5+Interm_revenue;
            end
            Interm_idx=find(Annual_Optimal_fixed_quantile_Revenues_case5==max(Annual_Optimal_fixed_quantile_Revenues_case5));
            eval(['Annual_best_fixed_quantile_Revenues_case5_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=Annual_Optimal_fixed_quantile_Revenues_case5;']);
            eval(['Annual_Optimal_fixed_quantile_Revenues_case5_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=Annual_Optimal_fixed_quantile_Revenues_case5(Interm_idx);']);
            Interm_idx_min=find(Annual_Optimal_fixed_quantile_Revenues_case5==min(Annual_Optimal_fixed_quantile_Revenues_case5));
            MaxFixedQuantileCase5(capacity_idx)=Interm_idx;
            MinFixedQuantileCase5(capacity_idx)=Interm_idx_min;            %             Annual_Optimal_fixed_quantile_Revenues_case5=Annual_Optimal_fixed_quantile_Revenues_case5(Interm_idx);


            eval(['Optimal_power_profile_Revenues_case5=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Revenue_opt']);
            Optimal_power_profile_Revenues_case5=cell2mat(Optimal_power_profile_Revenues_case5)/4; % /4 for MWh
            eval(['Annual_Optimal_power_profile_Revenues_case5_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(Optimal_power_profile_Revenues_case5);']);
            eval(['Actual_power_profile_Revenues_case5=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Revenue_best_guess_2018/4']); % /4 for MWh
            eval(['Annual_Actual_power_profile_Revenues_case5_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(Actual_power_profile_Revenues_case5)']);
            Pactual=NaN;
            Pmeasur=NaN;
            for days=1:365
                eval(['Pactual_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Pactual(days);']);
                eval(['Pmeasur_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.P_measur_day(days);']);
                Pactual_interm=cell2mat(Pactual_interm);
                Pmeasur_interm=cell2mat(Pmeasur_interm);
                Pactual=[Pactual;Pactual_interm];
                Pmeasur=[Pmeasur;Pmeasur_interm];
            end
            Pactual(1,:)=[];
            Pmeasur(1,:)=[];
            E_measur=trapz(Pmeasur)/4; % /4 for MWh
            E_actual=trapz(Pactual)/4; % /4 for MWh
            E_curt=E_measur-E_actual;
            Curtailements=Pmeasur-Pactual;

            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Pmeasur=Pmeasur;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Pactual=Pactual;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Curtailements=Curtailements;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_measur=E_measur;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_actual=E_actual;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_curt=E_curt;'])
        end
    elseif cases==6
        load('main_simulations_DTR_2018.mat')
        for capacity_idx=1:length(Overplanting_capacity)

            eval(['P50_Revenues_case6=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Revenue_P50_cdf_2018/4']); % /4 for MWh
            eval(['Annual_P50_Revenue_case6_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(P50_Revenues_case6)']);

            eval(['Variable_quantile_Revenues_case6_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.max_Revenue_cdf_2018/4']); % /4 for MWh
            eval(['Annual_Variable_quantile_Revenues_case6_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(Variable_quantile_Revenues_case6_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW)']);

            eval(['Optimal_fixed_quantile_Revenues_case6=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Revenue_cdf_2018  ']);
            Annual_Optimal_fixed_quantile_Revenues_case6=zeros(1,99);
            for i=1:365
                Interm_revenue=Optimal_fixed_quantile_Revenues_case6{i, 1}/4; % /4 for MWh
                Annual_Optimal_fixed_quantile_Revenues_case6=Annual_Optimal_fixed_quantile_Revenues_case6+Interm_revenue;
            end
            Interm_idx=find(Annual_Optimal_fixed_quantile_Revenues_case6==max(Annual_Optimal_fixed_quantile_Revenues_case6));
            Interm_idx_min=find(Annual_Optimal_fixed_quantile_Revenues_case6==min(Annual_Optimal_fixed_quantile_Revenues_case6));
            MaxFixedQuantileCase6(capacity_idx)=Interm_idx;
            MinFixedQuantileCase6(capacity_idx)=Interm_idx_min;
            eval(['Annual_best_fixed_quantile_Revenues_case6_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=Annual_Optimal_fixed_quantile_Revenues_case6;']);
            eval(['Annual_Optimal_fixed_quantile_Revenues_case6_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=Annual_Optimal_fixed_quantile_Revenues_case6(Interm_idx)']);
            %             Annual_Optimal_fixed_quantile_Revenues_case6=Annual_Optimal_fixed_quantile_Revenues_case6(Interm_idx);

            eval(['Optimal_power_profile_Revenues_case6=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Revenue_opt']);
            Optimal_power_profile_Revenues_case6=cell2mat(Optimal_power_profile_Revenues_case6)/4;% /4 for MWh
            eval(['Annual_Optimal_power_profile_Revenues_case6_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(Optimal_power_profile_Revenues_case6)']);
            eval(['Actual_power_profile_Revenues_case6_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Revenue_best_guess_2018/4']); % /4 for MWh
            eval(['Annual_Actual_power_profile_Revenues_case6_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(Actual_power_profile_Revenues_case6_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW)']);
            Pactual=NaN;
            Pmeasur=NaN;
            for days=1:365
                eval(['Pactual_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Pactual(days);']);
                eval(['Pmeasur_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.P_measur_day(days);']);
                Pactual_interm=cell2mat(Pactual_interm);
                Pmeasur_interm=cell2mat(Pmeasur_interm);
                Pactual=[Pactual;Pactual_interm];
                Pmeasur=[Pmeasur;Pmeasur_interm];
            end
            Pactual(1,:)=[];
            Pmeasur(1,:)=[];
            E_measur=trapz(Pmeasur)/4; % /4 for MWh
            E_actual=trapz(Pactual)/4; % /4 for MWh
            E_curt=E_measur-E_actual;
            Curtailements=Pmeasur-Pactual;

            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Pmeasur=Pmeasur;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Pactual=Pactual;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Curtailements=Curtailements;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_measur=E_measur;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_actual=E_actual;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_curt=E_curt;'])
        end
    end
end % end of cases 1:6


% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Figure 10 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Figure 10 Impact of commitment strategies on the annual revenue of a non-overplanted OWF: Case 2

% create figure
figure('DefaultAxesFontSize',14,'InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized')

%  Prepare the data for ploting
dataCase5=[Annual_P50_Revenue_case5_339_MW/Annual_P50_Revenue_case5_339_MW*100,...
    Annual_Optimal_fixed_quantile_Revenues_case5_339_MW/Annual_P50_Revenue_case5_339_MW*100,...
    Annual_Actual_power_profile_Revenues_case5_339_MW/Annual_P50_Revenue_case5_339_MW*100,...
    Annual_Variable_quantile_Revenues_case5_339_MW/Annual_P50_Revenue_case5_339_MW*100,...
    Annual_Optimal_power_profile_Revenues_case5_339_MW/Annual_P50_Revenue_case5_339_MW*100];

% Plot the bars
bar(dataCase5)

% Plot the number over the bars
text(1:length(dataCase5),round(dataCase5),num2str(round(dataCase5)'),'vert','bottom','horiz','center');

% Set the xlabels
set(gca,'xticklabel',{'P50 (Reference)','P_F_i_x_e_d_Q_u_a_n_t_i_l_e','P_a_c_t_u_a_l','P_V_ar_Q_u_a_n_t_i_l_e','P_O_p_t_i_m_P_r_o_f_i_l_e'})
ylabel('Annual revenue,reference % ')
xlabel('Commitment strategy')

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Figure 11 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Figure 11 Annual revenue as a function of fixed quantile for a non-overplanted OWF
% Normalize the revenue over the reference (P50 quantile), in %
Normalized_bestfixed_revenue=Annual_best_fixed_quantile_Revenues_case5_339_MW/Annual_best_fixed_quantile_Revenues_case5_339_MW(50)*100;% in %

% Create figure
figure('DefaultAxesFontSize',14,'InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized')

% Create axes
axes1 = axes('Position',[0.13 0.112375296912114 0.775 0.815]);
hold(axes1,'on');

% Create plot
plot(Normalized_bestfixed_revenue,'DisplayName','Fixed Quantile','LineWidth',3,'Color',[0 0 1]);

% Create scatter
scatter(MaxFixedQuantileCase6(1),Normalized_bestfixed_revenue(MaxFixedQuantileCase6(1)),'DisplayName','Maximal revenue Fixed Quantile',...
    'MarkerFaceColor',[1 0 0],...
    'MarkerEdgeColor',[0 0 0]);

% Create scatter
scatter(50,Normalized_bestfixed_revenue(50),'DisplayName','P50',...
    'MarkerFaceColor',[0.800000011920929 0.800000011920929 0.800000011920929],...
    'MarkerEdgeColor',[0 0 0]);

% Create ylabel
ylabel('Annual Revenue, % of P50 revenue');

% Create xlabel
xlabel('Fixed quantiles: P1 - P99');

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 100]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[211000000 219000000]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(axes1,[-1 1]);
box(axes1,'off');

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Figure 12 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Figure 12 Scatter plot showing the price and power differences for a non-overplanted OWF

P50=[];
Pactual=[];
P_cdf=[];

% Extract the power profiles in each day
for i=1:365
    % Extract a power profile for each quantile at i day
    P_intermP_cdf=Overplanting_results.Capacity_339_MW.TempConstraints.P_cdf{i,1};

    % Extract a power profile for P50 at i day
    P_intermP50=Overplanting_results.Capacity_339_MW.TempConstraints.Pplan_P50_Elia{i,1};

    % Extract actual power profile at i day
    P_intermPactual=Overplanting_results.Capacity_339_MW.TempConstraints.Pactual{i,1};

    % Add the power profile at i day to one vector
    P50=cat(1,P50,P_intermP50);
    P_cdf=cat(1,P_cdf,P_intermP_cdf);
    Pactual=cat(1,Pactual,P_intermPactual);
end

% Find the difference between forecasted profiles and actual power delivery
Difference=P_cdf-Pactual;
DifferencePower=P50-Pactual;

% Find the difference between P50 and P01/P99
DifferenceP50_P01=P50-P_cdf(:,1);
DifferenceP50_P99=P50-P_cdf(:,99);


% create a figure
figure('DefaultAxesFontSize',14,'InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized')
hold on
xlabel('Difference: PXX-Pfact, MW')
ylabel('Difference: Day-ahead price-Imbalance price, €')
zlabel('Frequency (out of 35040)')

% Load day-ahead and imbalance prices
load('data_DA_IMB_years_2015 2020.mat')

% Find the difference between day-ahaed price and imbalance + price
Diff2018plus=DA_2018-Cb_plus_2018; % in 2018

% Plot the histogram for P01, P50 and P99
h_P01=histogram2(Difference(:,1),Diff2018plus,[461 length(unique(round(Diff2018plus)))])
h_P99=histogram2(Difference(:,99),Diff2018plus,[461 length(unique(round(Diff2018plus)))])
h_P50=histogram2(DifferencePower,Diff2018plus,[461 length(unique(round(Diff2018plus)))])

legend('P01','P99','P50')



% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Figure 13 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Figure 13 Imbalance+ prices versus DA prices during the considered period in France

% create a figure
figure('DefaultAxesFontSize',14,'InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized')

% Plot the histogram: Day-ahaed vs imbalance price
histogram2(DA_2018,Cb_plus_2018,'FaceColor','flat')
xlabel('Day-ahead prices,€/MWh')
ylabel('Imbalance prices,€/MWh')
view(2) % show X-Y view


% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Figure 14-15 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Figure 14-15 Example of the day when DA prices were below/higher the imbalance prices
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  YMATRIX2:  matrix of y data
load('Elia_Jan13_2016_Jan13_2017_powers.mat','t_year')

% Index of the VarQuantile (the best quantile to trade at the given day)
idx=Overplanting_results.Capacity_339_MW.TempConstraints.Quantile2018;

for i=1:365
    if ~(length(idx{i, 1}  )==1) % if there are few quantiles
        % Save the vector
        idx_interm1=cell2mat(idx(i));

        % Take the first value in the vector
        idx_interm(i)=idx_interm1(1);
    else % otherwise (one quantile)
        idx_interm(i)=cell2mat(idx(i));

    end % end if
end % end for cycle

% Save idx in double format
idx=idx_interm;

% Create a variable
PvarQuant=[];
P_optimProfil=[];

for i=1:365
    % Extact all quantiles at given day
    P_cdf_interm=Overplanting_results.Capacity_339_MW.TempConstraints.P_cdf{i, 1};
    % extract the best quantile at the given day
    PvarQuant_interm=P_cdf_interm(:,idx(i));
    % Save in one vector
    PvarQuant=cat(1,PvarQuant,PvarQuant_interm);

    % Extact optimal power profile  at given day
    PoptimProf_interm=Overplanting_results.Capacity_339_MW.TempConstraints.P_plan_opt{i, 1};
    % Save in one vector
    P_optimProfil=cat(1,P_optimProfil,PoptimProf_interm);
end


% Create figure
figure('DefaultAxesFontSize',14,'InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized')

% Create subplot
subplot1 = subplot(2,1,1);
hold(subplot1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(time_2018,[P_cdf(:,90),P_cdf(:,50),P_cdf(:,10),PvarQuant,Pactual,P_optimProfil],'LineWidth',2);
set(plot1(1),'DisplayName','P90');
set(plot1(2),'DisplayName','P50');
set(plot1(3),'DisplayName','P10');
set(plot1(4),'DisplayName','P_V_a_r_Q_u_a_n_t_i_l_e','LineWidth',1);
set(plot1(5),'DisplayName','P_a_c_t_u_a_l','LineWidth',3,'Color','c');
set(plot1(6),'DisplayName','P_O_p_t_i_m_P_r_o_f_i_l_e',...
    'Color',[0.466666668653488 0.674509823322296 0.18823529779911]);

% Create ylabel
ylabel('Power, MW','FontSize',11);

% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot1,[-18.1195999052204 352.06790737855]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot1,[-1 1]);
hold(subplot1,'off');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.195618864856791 0.760664085059999 0.124669138472695 0.100059384631431],...
    'NumColumns',2,...
    'EdgeColor',[1 1 1]);

% Create subplot
subplot2 = subplot(2,1,2);
hold(subplot2,'on');

% Create multiple lines using matrix input to plot
plot2 = plot(time_2018,[DA_2018,Cb_plus_2018,Cb_minus_2018],'LineWidth',2);
set(plot2(1),'DisplayName','Day-Ahead');
set(plot2(2),'DisplayName','Imbalance+');
set(plot2(3),'DisplayName','Imbalance-');

% Create ylabel
ylabel('Price, €/MWh');

% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot2,[0 80]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot2,[-1 1]);
hold(subplot2,'off');
% Create legend
legend2 = legend(subplot2,'show');
set(legend2,...
    'Position',[0.195304550186496 0.330661122608261 0.0750397045313204 0.0786817119127214],...
    'EdgeColor',[1 1 1]);

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Figure 16 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Figure 16 Annual revenue as a function of overplanting rate and commitment strategies, with STR (left) and DTR (right)

% create a figure
figure('DefaultAxesFontSize',14,'InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized')
subplot(1,2,1)
ylabel('Annual Revenu,€')
% case5: preparing the data
dataCase5=[Annual_P50_Revenue_case5_339_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Variable_quantile_Revenues_case5_339_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case5_339_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_power_profile_Revenues_case5_339_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Actual_power_profile_Revenues_case5_339_MW/Annual_P50_Revenue_case5_339_MW*100;...
    Annual_P50_Revenue_case5_373_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Variable_quantile_Revenues_case5_373_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case5_373_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_power_profile_Revenues_case5_373_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Actual_power_profile_Revenues_case5_373_MW/Annual_P50_Revenue_case5_339_MW*100;...
    Annual_P50_Revenue_case5_407_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Variable_quantile_Revenues_case5_407_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case5_407_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_power_profile_Revenues_case5_407_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Actual_power_profile_Revenues_case5_407_MW/Annual_P50_Revenue_case5_339_MW*100;...
    Annual_P50_Revenue_case5_441_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Variable_quantile_Revenues_case5_441_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case5_441_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_power_profile_Revenues_case5_441_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Actual_power_profile_Revenues_case5_441_MW/Annual_P50_Revenue_case5_339_MW*100;...
    Annual_P50_Revenue_case5_509_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Variable_quantile_Revenues_case5_509_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case5_509_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_power_profile_Revenues_case5_509_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Actual_power_profile_Revenues_case5_509_MW/Annual_P50_Revenue_case5_339_MW*100;...
    Annual_P50_Revenue_case5_577_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Variable_quantile_Revenues_case5_577_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case5_577_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_power_profile_Revenues_case5_577_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Actual_power_profile_Revenues_case5_577_MW/Annual_P50_Revenue_case5_339_MW*100;...
    Annual_P50_Revenue_case5_645_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Variable_quantile_Revenues_case5_645_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case5_645_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_power_profile_Revenues_case5_645_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Actual_power_profile_Revenues_case5_645_MW/Annual_P50_Revenue_case5_339_MW*100;...
    Annual_P50_Revenue_case5_679_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Variable_quantile_Revenues_case5_679_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case5_679_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Optimal_power_profile_Revenues_case5_679_MW/Annual_P50_Revenue_case5_339_MW*100,Annual_Actual_power_profile_Revenues_case5_679_MW/Annual_P50_Revenue_case5_339_MW*100];

% bar(dataCase5)

% Changing the column places for better view
dataCase5_dash=dataCase5;
dataCase5_dash(:,2)=dataCase5(:,3);
dataCase5_dash(:,3)=dataCase5(:,5);
dataCase5_dash(:,4)=dataCase5(:,2);
dataCase5_dash(:,5)=dataCase5(:,4);
dataCase5=dataCase5_dash;

% Plot bars
bar3(dataCase5)

% Display the bar values at their top
[X,Y] = meshgrid(1:size(dataCase5,2), 1:size(dataCase5,1));
text(X(:), Y(:), dataCase5(:), num2str(round(dataCase5(:))), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')

set(gca,'xticklabel',{'P50 strategy','Best fixed quantile', 'Actual power profile','Variable quantile','Optimal power profile'})
set(gca,'yticklabel',Overplanting_rate)
ylabel('Overplanting level,pu')
zlabel('Annual revenue,%  of reference')
xlabel('Market strategy')


% create a figure
subplot(1,2,2)
ylabel('Annual Revenue,€')

% Preparing the data
dataCase6=[Annual_P50_Revenue_case6_339_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Variable_quantile_Revenues_case6_339_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case6_339_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_power_profile_Revenues_case6_339_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Actual_power_profile_Revenues_case6_339_MW/Annual_P50_Revenue_case6_339_MW*100;...
    Annual_P50_Revenue_case6_373_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Variable_quantile_Revenues_case6_373_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case6_373_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_power_profile_Revenues_case6_373_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Actual_power_profile_Revenues_case6_373_MW/Annual_P50_Revenue_case6_339_MW*100;...
    Annual_P50_Revenue_case6_407_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Variable_quantile_Revenues_case6_407_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case6_407_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_power_profile_Revenues_case6_407_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Actual_power_profile_Revenues_case6_407_MW/Annual_P50_Revenue_case6_339_MW*100;...
    Annual_P50_Revenue_case6_441_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Variable_quantile_Revenues_case6_441_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case6_441_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_power_profile_Revenues_case6_441_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Actual_power_profile_Revenues_case6_441_MW/Annual_P50_Revenue_case6_339_MW*100;...
    Annual_P50_Revenue_case6_509_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Variable_quantile_Revenues_case6_509_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case6_509_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_power_profile_Revenues_case6_509_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Actual_power_profile_Revenues_case6_509_MW/Annual_P50_Revenue_case6_339_MW*100;...
    Annual_P50_Revenue_case6_577_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Variable_quantile_Revenues_case6_577_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case6_577_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_power_profile_Revenues_case6_577_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Actual_power_profile_Revenues_case6_577_MW/Annual_P50_Revenue_case6_339_MW*100;...
    Annual_P50_Revenue_case6_645_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Variable_quantile_Revenues_case6_645_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case6_645_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_power_profile_Revenues_case6_645_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Actual_power_profile_Revenues_case6_645_MW/Annual_P50_Revenue_case6_339_MW*100;...
    Annual_P50_Revenue_case6_679_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Variable_quantile_Revenues_case6_679_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_fixed_quantile_Revenues_case6_679_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Optimal_power_profile_Revenues_case6_679_MW/Annual_P50_Revenue_case6_339_MW*100,Annual_Actual_power_profile_Revenues_case6_679_MW/Annual_P50_Revenue_case6_339_MW*100];

% Changing the column place for better view
dataCase6_dash=dataCase6;
dataCase6_dash(:,2)=dataCase6(:,3);
dataCase6_dash(:,3)=dataCase6(:,5);
dataCase6_dash(:,4)=dataCase6(:,2);
dataCase6_dash(:,5)=dataCase6(:,4);
dataCase6=dataCase6_dash;

% Plot bars
bar3(dataCase6)

% Displaythe bar alues at the top
[X,Y] = meshgrid(1:size(dataCase6,2), 1:size(dataCase6,1));
text(X(:), Y(:), dataCase6(:), num2str(round(dataCase6(:))), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')

set(gca,'xticklabel',{'P50 strategy','Best fixed quantile', 'Actual power profile','Variable quantile','Optimal power profile'})
set(gca,'yticklabel',Overplanting_rate)
ylabel('Overplanting level,pu')
zlabel('Annual revenue,%  of reference')
xlabel('Market strategy')

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Figure 17 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Figure 17 Annual revenue for PFixedQuantile (top). Best fixed quantile as a function of the overplanting rate (middle). Revenue difference
% between the best-fixed quantile strategy and the P50 strategy (bottom).

% Create a figure
figure('DefaultAxesFontSize',14,'InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized')

subplot(3,1,1)
hold on

% Calculate mean, max and min revenue for the best fixed quantil strategy
for capacity_idx=1:length(Overplanting_capacity)
    % case 5 (current constraints)
    eval(['Mean_valueCase5(capacity_idx)=mean(Annual_best_fixed_quantile_Revenues_case5_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW);']);
    eval(['Max_valueCase5(capacity_idx)=max(Annual_best_fixed_quantile_Revenues_case5_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW);']);
    eval(['Min_valueCase5(capacity_idx)=min(Annual_best_fixed_quantile_Revenues_case5_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW);']);

    % case 6 (temperature constraints)
    eval(['Mean_valueCase6(capacity_idx)=mean(Annual_best_fixed_quantile_Revenues_case6_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW);']);
    eval(['Max_valueCase6(capacity_idx)=max(Annual_best_fixed_quantile_Revenues_case6_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW);']);
    eval(['Min_valueCase6(capacity_idx)=min(Annual_best_fixed_quantile_Revenues_case6_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW);']);

end

% Plot the bars with errors
errorbar(Overplanting_rate,(Mean_valueCase5/Annual_best_fixed_quantile_Revenues_case5_339_MW(50))*100,(Mean_valueCase5/Annual_best_fixed_quantile_Revenues_case5_339_MW(50))*100-(Min_valueCase5/Annual_best_fixed_quantile_Revenues_case5_339_MW(50))*100,(Max_valueCase5/Annual_best_fixed_quantile_Revenues_case5_339_MW(50))*100-(Mean_valueCase5/Annual_best_fixed_quantile_Revenues_case5_339_MW(50))*100,'x')
errorbar(Overplanting_rate,(Mean_valueCase6/Annual_best_fixed_quantile_Revenues_case6_339_MW(50))*100,(Mean_valueCase6/Annual_best_fixed_quantile_Revenues_case5_339_MW(50))*100-(Min_valueCase6/Annual_best_fixed_quantile_Revenues_case6_339_MW(50))*100,(Max_valueCase6/Annual_best_fixed_quantile_Revenues_case6_339_MW(50))*100-(Mean_valueCase6/Annual_best_fixed_quantile_Revenues_case6_339_MW(50))*100,'x')
ylabel('Annual Revenue, % of ref')
xlabel('Overplanting rate')
legend('Case3','Case4')

subplot(3,1,2)
hold on
% Plot the best quantile (maximizing the revenue over the year)
plot(Overplanting_rate,MaxFixedQuantileCase5,'-^','LineWidth',0.5)
plot(Overplanting_rate,MaxFixedQuantileCase6,'-^','LineWidth',0.5)
legend('Case3','Case4')
ylabel('Best Quantile')
xlabel('Overplanting rate')

subplot(3,1,3)
% prepare the data with difference against P50
a=[Max_valueCase5(1)-Annual_P50_Revenue_case5_339_MW,Max_valueCase5(2)-Annual_P50_Revenue_case5_373_MW,Max_valueCase5(3)-Annual_P50_Revenue_case5_407_MW,Max_valueCase5(4)-Annual_P50_Revenue_case5_441_MW,Max_valueCase5(5)-Annual_P50_Revenue_case5_509_MW,Max_valueCase5(6)-Annual_P50_Revenue_case5_577_MW,Max_valueCase5(7)-Annual_P50_Revenue_case5_645_MW,Max_valueCase5(8)-Annual_P50_Revenue_case5_679_MW]/Annual_best_fixed_quantile_Revenues_case5_339_MW(50)*100;
b=[Max_valueCase6(1)-Annual_P50_Revenue_case6_339_MW,Max_valueCase6(2)-Annual_P50_Revenue_case6_373_MW,Max_valueCase6(3)-Annual_P50_Revenue_case6_407_MW,Max_valueCase6(4)-Annual_P50_Revenue_case6_441_MW,Max_valueCase6(5)-Annual_P50_Revenue_case6_509_MW,Max_valueCase6(6)-Annual_P50_Revenue_case6_577_MW,Max_valueCase6(7)-Annual_P50_Revenue_case6_645_MW,Max_valueCase6(8)-Annual_P50_Revenue_case6_679_MW]/Annual_best_fixed_quantile_Revenues_case5_339_MW(50)*100;

% Plot the bar
hB=bar(Overplanting_rate,[a;b]);
legend(' Case3','Case4','Location','northwest'	)
xlabel('Overplanting rate')
ylabel('Difference with P50,%')

% Pltoing the bar values at the top
hAx=gca;            % get a variable for the current axes handle
hT=[];              % placeholder for text object handles
for i=1:length(hB)  % iterate over number of bar objects
    hT=[hT text(hB(i).XData+hB(i).XOffset,hB(i).YData,num2str(hB(i).YData.','%.1f'), ...
        'VerticalAlignment','bottom','horizontalalign','center')];
end

%% Table 3 
% To add the code later 


%% Figure 18 - 19 + Table 4

% Clear workspace
clear all

% Uncomment the neccesary cases and year:

% Choosing the DTR case
% cases=5; %  current constraint
cases=6; % temperature constraint

% Choosing the year for market prices
% year = 2018; % Uncomment
year = 2022; % Uncomment

% Note that some code below would still use 2018 even though 2022 is
% chosen as the year. This is not an error as simulations saved results for
% 2022 as 2018. The right year data would be defined by filename below

% Choose the name of precalculated data
if cases==5 % STR
    if year == 2018
        filename=sprintf('main_simulations_STR_2018.mat');
    elseif year == 2022
        filename=sprintf('main_simulations_STR_2022.mat');
    else
        error('Choose the year either 2018 or 2022')
    end
elseif cases == 6 % DTR
    if year == 2018
        filename=sprintf('main_simulations_DTR_2018.mat');

    elseif year == 2022
        filename=sprintf('main_simulations_DTR_2022.mat');
    else
        error('Choose the year either 2018 or 2022')
    end
else
    error('Choose the case either 5 or 6')
end


% Load data
load(filename)

%  Set of  OWF capacities
Overplanting_capacity=[339.4387  373.3825  407.3264  441.2703  509.1580...
    577.0457  644.9334  678.8773]; % MW

% Set a feed-in tariff of offshore wind farm
FIT=69.5513; % €/MWh % LCOE for 1.3 pu
FIT_Saint_Nazaire=143; % €/MWh % FIT for first offshore wind farm in France (at Saint Nazaire)

% discount rate in % for NPV calculations
discount_rate=2.5; %

% Project horizon for NPV calculations
Project_horizon=27; % years

% Share of cable export costs in all cable installation costs(for all cables:inner and export):
coef_export=0.5;

% GPB/EUR 1.1405 in 2019
exchange_rate_GBP_EUR=1.1405;

% Calculate the CAPEX and OPEX (€) for each overplanting rate
[CAPEX,OPEX]=CATAPULT_costs_2019(coef_export,exchange_rate_GBP_EUR);

% CATAPULT costs : https://guidetoanoffshorewindfarm.com/wind-farm-costs

% Extracting the annual revenue for overplanted OWF
% For cycle: extracting the annual revenues
E_measurs=[];
E_actuals=[];

for capacity_idx=1:length(Overplanting_capacity)
    % Extract the revenu for P50 strategy. Note that /4 is neccesary to
    % normalize € from MW15min (calculated initially) to € for MWh
    if cases==3 || cases==6
        eval(['P50_Revenues=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Revenue_P50_cdf_2018/4']);% €
    else  % other cases 1,2,4,5
        eval(['P50_Revenues=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Revenue_P50_cdf_2018/4']);% €
    end

    % Calculate the annual revenue
    eval(['Annual_P50_Revenue_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(P50_Revenues)']);

    % Extract the revenu for variable quantile strategy. Note that /4 is
    % neccesary to % normalize € from MW15min (calculated initially) to € for MWh
    if cases==3 || cases==6
        eval(['Variable_quantile_Revenues_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.max_Revenue_cdf_2018/4']);% MWh
    else % other cases 1,2,4,5
        eval(['Variable_quantile_Revenues_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.max_Revenue_cdf_2018/4']);% MWh
    end

    % Calculate the annual revenue
    eval(['Annual_Variable_quantile_Revenues_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(Variable_quantile_Revenues_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW)']);

    if cases==3 || cases==6
        % Extract the revenu for variable quantile strategy. Note that /4 is
        eval(['Optimal_fixed_quantile_Revenues=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Revenue_cdf_2018']);
    else % other cases 1,2,4,5
        eval(['Optimal_fixed_quantile_Revenues=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Revenue_cdf_2018']);

    end
    % Prepare the zeros vector
    Annual_Optimal_fixed_quantile_Revenues=zeros(1,99);
    for i=1:365
        % Extract the revenue for i day
        Interm_revenue=Optimal_fixed_quantile_Revenues{i, 1}/4; % /4 is neccesary to % normalize € from MW15min (calculated initially) to € for MWh

        % Update the accumulated revenue for the given day
        Annual_Optimal_fixed_quantile_Revenues=Annual_Optimal_fixed_quantile_Revenues+Interm_revenue;
    end

    % Find the index corresponding to the maximal annual revenue
    Interm_idx=find(Annual_Optimal_fixed_quantile_Revenues==max(Annual_Optimal_fixed_quantile_Revenues));

    % Find the index corresponding to the minimal annual revenue
    Interm_idx_min=find(Annual_Optimal_fixed_quantile_Revenues==min(Annual_Optimal_fixed_quantile_Revenues));

    %     MaxFixedQuantileCase6(capacity_idx)=Interm_idx;
    %     MinFixedQuantileCase6(capacity_idx)=Interm_idx_min;

    % Change variable
    eval(['Annual_best_fixed_quantile_Revenues_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=Annual_Optimal_fixed_quantile_Revenues;']);

    % Save the highest annual revenue of the fixed quantile strategy
    eval(['Annual_Optimal_fixed_quantile_Revenues_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=Annual_Optimal_fixed_quantile_Revenues(Interm_idx)']);

    if cases==3 || cases==6
        % Extract the revenu for optimal power profile strategy.
        eval(['Optimal_power_profile_Revenues=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Revenue_opt ']);
    else % other cases 1,2,4,5
        eval(['Optimal_power_profile_Revenues=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Revenue_opt ']);
    end

    % Convert Cell to double
    Optimal_power_profile_Revenues=cell2mat(Optimal_power_profile_Revenues)/4; % /4 is neccesary to % normalize € from MW15min (calculated initially) to € for MWh

    % Calculate the annual revenue
    eval(['Annual_Optimal_power_profile_Revenues_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(Optimal_power_profile_Revenues)']);

    if cases==3 || cases==6
        % Extract the revenu for actual power profile strategy. Note that /4 is
        % neccesary to % normalize € from MW15min (calculated initially) to € for MWh
        eval(['Actual_power_profile_Revenues_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Revenue_best_guess_2018/4  ']);
    else % other cases 1,2,4,5
        eval(['Actual_power_profile_Revenues_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Revenue_best_guess_2018/4  ']);
    end

    % Calculate the annual revenue
    eval(['Annual_Actual_power_profile_Revenues_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(Actual_power_profile_Revenues_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW)']);

    % Create NaN variable
    Pactual=NaN; %  power profile after curtailements
    Pmeasur=NaN; %  power profile without curtailement (Elia load factor  * overplanting capacity)
    P50=NaN;    % power profile for P50 without curtailement (Elia load factor P50 forecast * overplanting capacity)

    % For each day  extract power and price profiles as well as  energy values
    for days=1:365
        if cases==3 || cases==6

            % Extract power values for the given day
            eval(['Pactual_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Pactual(days);']);
            eval(['Pmeasur_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.P_measur_day(days);']);
            eval(['P50_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Pplan_P50_Elia(days);']);
        else % other cases
            % Extract power values for the given day
            eval(['Pactual_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Pactual(days);']);
            eval(['Pmeasur_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.P_measur_day(days);']);
            eval(['P50_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Pplan_P50_Elia(days);']);
        end
        % Convert cell data to double data
        Pactual_interm=cell2mat(Pactual_interm);
        Pmeasur_interm=cell2mat(Pmeasur_interm);
        P50_interm=cell2mat(P50_interm);

        % Combine the power profile for the given day with previous days
        Pactual=[Pactual;Pactual_interm];
        Pmeasur=[Pmeasur;Pmeasur_interm];
        P50=[P50;P50_interm];

        if cases==3 || cases==6
            % Extract prices for the given day
            eval(['Cb_minus_day_2018=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Cb_minus_day_',num2str(year),'(days);']);
            eval(['Cb_plus_day_2018=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Cb_plus_day_',num2str(year),'(days);']);
            eval(['DAh_day_2018=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.DA_day_',num2str(year),'(days);']);
        else % other cases
            eval(['Cb_minus_day_2018=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Cb_minus_day_',num2str(year),'(days);']);
            eval(['Cb_plus_day_2018=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Cb_plus_day_',num2str(year),'(days);']);
            eval(['DAh_day_2018=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.DA_day_',num2str(year),'(days);']);

        end
        % Convert cell data to double data
        Cb_minus_day_2018=cell2mat(Cb_minus_day_2018);
        Cb_plus_day_2018=cell2mat(Cb_plus_day_2018);
        DAh_day_2018=cell2mat(DAh_day_2018);

        % Create the vector of feed-in tariff
        FIT_LCOE=linspace(FIT,FIT,96)';
        FIT_SN=linspace(FIT_Saint_Nazaire,FIT_Saint_Nazaire,96)';

        % Caclulate the revenue for given FIT
        Revenue_FIT_SN_Pactual(days)=sum(Pactual_interm.*FIT_SN-Cb_minus_day_2018.*max(0,Pactual_interm-Pactual_interm)+Cb_plus_day_2018.*max(0,Pactual_interm-Pactual_interm));

        % Caclulate the revenue for given FIT
        Revenue_FIT_Pactual(days)=sum(Pactual_interm.*FIT_LCOE-Cb_minus_day_2018.*max(0,Pactual_interm-Pactual_interm)+Cb_plus_day_2018.*max(0,Pactual_interm-Pactual_interm));

    end

    % Delete NaN values
    Pactual(1,:)=[];
    Pmeasur(1,:)=[];
    P50(1,:)=[];

    % Calculate the curtailed energy
    E_measur=trapz(Pmeasur)/4;  % MWh without curtailement
    E_actual=trapz(Pactual)/4;  % MWh with curtailement to ensure a cable limit
    E_max=Overplanting_capacity(capacity_idx)*8760;
    E_curt=E_measur-E_actual; % MWh
    E_measurs(capacity_idx,1)=E_measur;
    E_actuals(capacity_idx,1)=E_actual;

    % Capacity_factor
    eval(['Capacity_factors.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=E_actual/E_max'])

    % Calculate curtailement profile
    Curtailements=Pmeasur-Pactual;

    % Save results: power and energy (PandE)
    eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Pmeasur=Pmeasur;'])
    eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Pactual=Pactual;'])
    eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Curtailements=Curtailements;'])
    eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_measur=E_measur;'])
    eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_actual=E_actual;'])
    eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_curt=E_curt;'])

    % Calculate the annual revenue for FIT
    eval(['Annual_FIT_SN_Actual_Revenues_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(Revenue_FIT_SN_Pactual/4);' ]);
    eval(['Annual_FIT_Actual_Revenues_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(Revenue_FIT_Pactual/4);' ]);

end

% Caclulating cash flows

% Preparing the data
Annual_revenues_P50=[Annual_P50_Revenue_339_MW Annual_P50_Revenue_373_MW...
    Annual_P50_Revenue_407_MW Annual_P50_Revenue_441_MW...
    Annual_P50_Revenue_509_MW Annual_P50_Revenue_577_MW...
    Annual_P50_Revenue_645_MW Annual_P50_Revenue_679_MW];

Annual_revenues_fixed_quantile=[Annual_Optimal_fixed_quantile_Revenues_339_MW Annual_Optimal_fixed_quantile_Revenues_373_MW...
    Annual_Optimal_fixed_quantile_Revenues_407_MW Annual_Optimal_fixed_quantile_Revenues_441_MW...
    Annual_Optimal_fixed_quantile_Revenues_509_MW Annual_Optimal_fixed_quantile_Revenues_577_MW...
    Annual_Optimal_fixed_quantile_Revenues_645_MW Annual_Optimal_fixed_quantile_Revenues_679_MW];

Annual_revenues_actual_power=[Annual_Actual_power_profile_Revenues_339_MW Annual_Actual_power_profile_Revenues_373_MW...
    Annual_Actual_power_profile_Revenues_407_MW Annual_Actual_power_profile_Revenues_441_MW...
    Annual_Actual_power_profile_Revenues_509_MW Annual_Actual_power_profile_Revenues_577_MW...
    Annual_Actual_power_profile_Revenues_645_MW Annual_Actual_power_profile_Revenues_679_MW];

Annual_revenues_variable_quantile=[Annual_Variable_quantile_Revenues_339_MW Annual_Variable_quantile_Revenues_373_MW...
    Annual_Variable_quantile_Revenues_407_MW Annual_Variable_quantile_Revenues_441_MW...
    Annual_Variable_quantile_Revenues_509_MW Annual_Variable_quantile_Revenues_577_MW...
    Annual_Variable_quantile_Revenues_645_MW Annual_Variable_quantile_Revenues_679_MW];

Annual_revenues_optimal_power=[Annual_Optimal_power_profile_Revenues_339_MW Annual_Optimal_power_profile_Revenues_373_MW...
    Annual_Optimal_power_profile_Revenues_407_MW Annual_Optimal_power_profile_Revenues_441_MW...
    Annual_Optimal_power_profile_Revenues_509_MW Annual_Optimal_power_profile_Revenues_577_MW...
    Annual_Optimal_power_profile_Revenues_645_MW Annual_Optimal_power_profile_Revenues_679_MW];

Annual_revenues_FIT_SN_Actual_Revenues=[Annual_FIT_SN_Actual_Revenues_339_MW Annual_FIT_SN_Actual_Revenues_373_MW...
    Annual_FIT_SN_Actual_Revenues_407_MW Annual_FIT_SN_Actual_Revenues_441_MW...
    Annual_FIT_SN_Actual_Revenues_509_MW Annual_FIT_SN_Actual_Revenues_577_MW...
    Annual_FIT_SN_Actual_Revenues_645_MW Annual_FIT_SN_Actual_Revenues_679_MW];

Annual_revenues_FIT_Actual=[Annual_FIT_Actual_Revenues_339_MW Annual_FIT_Actual_Revenues_373_MW...
    Annual_FIT_Actual_Revenues_407_MW Annual_FIT_Actual_Revenues_441_MW...
    Annual_FIT_Actual_Revenues_509_MW Annual_FIT_Actual_Revenues_577_MW...
    Annual_FIT_Actual_Revenues_645_MW Annual_FIT_Actual_Revenues_679_MW];

% Calculating the cash flows (not discounted yet)
CF_P50=Annual_revenues_P50-OPEX;
CF_fixed_quantile=Annual_revenues_fixed_quantile-OPEX;
CF_actual_power=Annual_revenues_actual_power-OPEX;
CF_variable_quantile=Annual_revenues_variable_quantile-OPEX;
CF_optimal_power=Annual_revenues_optimal_power-OPEX;
CF_FIT_SN_Actual=Annual_revenues_FIT_SN_Actual_Revenues-OPEX;
CF_FIT_Actual=Annual_revenues_FIT_Actual-OPEX;

% Calculate discounted cash flows
for t=1:Project_horizon % for 1 to n years
    CF_discounted_P50(t,:)=CF_P50/((1+discount_rate/100)^t);
    CF_discounted_BestFixed(t,:)=CF_fixed_quantile/((1+discount_rate/100)^t);
    CF_discounted_ActualP(t,:)=CF_actual_power/((1+discount_rate/100)^t);
    CF_discounted_VarQuant(t,:)=CF_variable_quantile/((1+discount_rate/100)^t);
    CF_discounted_Optimal(t,:)=CF_optimal_power/((1+discount_rate/100)^t);
    CF_discounted_FIT_SN_Actual(t,:)=CF_FIT_SN_Actual/((1+discount_rate/100)^t);
    CF_discounted_FIT_Actual(t,:)=CF_FIT_Actual/((1+discount_rate/100)^t);

end

% Calculate NPV
for i=1:length(Overplanting_capacity)
    NPV_P50(1,i)=-CAPEX(i)+sum(CF_discounted_P50(:,i));
    NPV_BestFixed(1,i)=-CAPEX(i)+sum(CF_discounted_BestFixed(:,i));
    NPV_ActualP(1,i)=-CAPEX(i)+sum(CF_discounted_ActualP(:,i));
    NPV_VarQuant(1,i)=-CAPEX(i)+sum(CF_discounted_VarQuant(:,i));
    NPV_Optimal(1,i)=-CAPEX(i)+sum(CF_discounted_Optimal(:,i));
    NPV_FIT_SN_Actual(1,i)=-CAPEX(i)+sum(CF_discounted_FIT_SN_Actual(:,i));
    NPV_FIT_Actual(1,i)=-CAPEX(i)+sum(CF_discounted_FIT_Actual(:,i));

end

% Calculate LCOE
for capacity_idx=1:length(Overplanting_capacity)
    OPEX_tot=0;
    E_tot=0;
    eval(['E_actual=PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_actual;'])
    for t=1:Project_horizon
        OPEX_tot=OPEX_tot+OPEX(capacity_idx)/((1+discount_rate/100)^t);
        E_tot=E_tot+E_actual/((1+discount_rate/100)^t);
    end
    LCOE(capacity_idx)=(CAPEX(capacity_idx)+OPEX_tot)/(E_tot);
end

% Calculating the payback period
for strategy=1:7
    for capacity_idx=1:length(Overplanting_capacity)
        t=0;
        Revenue_discount=0;
        if strategy==1
            CF=CF_P50(capacity_idx);
        elseif strategy==2
            CF=CF_fixed_quantile(capacity_idx);
        elseif strategy==3
            CF=CF_actual_power(capacity_idx);
        elseif strategy==4
            CF=CF_variable_quantile(capacity_idx);
        elseif strategy==5
            CF=CF_optimal_power(capacity_idx);
        elseif strategy==6
            CF=CF_FIT_Actual(capacity_idx);
        else
            CF=CF_FIT_SN_Actual(capacity_idx);
        end

        while CAPEX(capacity_idx)>Revenue_discount
            t=t+1;
            if Revenue_discount==Revenue_discount+CF/((1+discount_rate/100)^t)
                t=inf;
                break
            else
                Revenue_discount=Revenue_discount+CF/((1+discount_rate/100)^t);

            end
        end

        Payback_period(capacity_idx,strategy)=t;
    end
end


Table_outputs=table();

Table_outputs.Overplanting_rate=Overplanting_rate';
Table_outputs.Payback_P50=Payback_period(:,1);
Table_outputs.Payback_FixedQuantile=Payback_period(:,2);
Table_outputs.Payback_Actual_Power=Payback_period(:,3);
Table_outputs.Payback_VariableQuantile=Payback_period(:,4);
Table_outputs.Payback_OptimalPower=Payback_period(:,5);
Table_outputs.Payback_FIT_LCOE=Payback_period(:,6);
Table_outputs.Payback_FIT_SN=Payback_period(:,7);

Table_outputs

% Displaying the results
figure('DefaultAxesFontSize',14,'InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized')
hold on
plot(Overplanting_rate,NPV_P50/1e6,'LineWidth',2,DisplayName='P50') % /1e6 - transformation to millions
plot(Overplanting_rate,NPV_BestFixed/1e6,'LineWidth',2,DisplayName='P_F_i_x_e_d_Q_u_a_n_t_i_l_e')
plot(Overplanting_rate,NPV_ActualP/1e6,'LineWidth',2,DisplayName='P_a_c_t_u_a_l')
plot(Overplanting_rate,NPV_VarQuant/1e6,'LineWidth',2,DisplayName='P_V_a_r_Q_u_a_n_t_i_l_e')
plot(Overplanting_rate,NPV_Optimal/1e6,'LineWidth',2,DisplayName='P_O_p_t_i_m_P_r_o_f_i_l_e')
plot(Overplanting_rate,NPV_FIT_SN_Actual/1e6,'LineWidth',2,DisplayName='FIT 143 €/MWh')
plot(Overplanting_rate,NPV_FIT_Actual/1e6,'LineWidth',2,DisplayName='FIT 69.6 €/ MWh')
xlabel('Overplanting rate,pu')
ylabel('Net present value,M€')
legend show
ylim([-2000 10000])

figure('InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized');
plot(Overplanting_rate,LCOE,'LineWidth',2,DisplayName='LCOE') % LCOE
xlabel('Overplanting rate')
ylabel('LCOE,€/MWh')
legend show
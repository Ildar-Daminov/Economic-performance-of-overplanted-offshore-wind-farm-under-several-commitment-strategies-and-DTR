function [CAPEX,OPEX]=CATAPULT_costs_2019(coef_export,exchange_rate_GBP_EUR)
%% Puppose 
% This function calculates CAPEX and OPEX for the given OWF capacity

% Input
% coef_export - Share of cable export costs in all cable installation costs
%               (for all cables:inner and export):
% exchange_rate_GBP_EUR - exchange rate for GBP_EUR in 2019

% Output: 
% CAPEX - Capital expenditures for each overplanting rate, EUR 
% OPEX - Operational (annual) expenditures for each overplanting rate, EUR 
%% Economic data from CATAPULT
% Source https://guidetoanoffshorewindfarm.com/wind-farm-costs

%  Set of overplanting capacities
Overplanting_capacity=[339.4387  373.3825  407.3264  441.2703  509.1580...
    577.0457  644.9334  678.8773]; % MW

%% Category	Rounded cost (£/MW)

% Development and project management
C_develop_pm=120000;% £/MW

% Turbine
C_turbine=1000000;% £/MW

%% Balance of plant
% C_BoP=600000;

% Export cable
C_export_cable=130000; % £/MW always for 339 MW

% Array cable
C_array_cable=35000; % £/MW

% Cable protection
C_cable_protection=2000; % % £/MW Export cable?

% Turbine foundation
C_turbine_foundation=280000;% £/MW

% Offshore substation
C_offshore_substation=120000;% £/MW

%  Onshore substation
C_onshore_substatino=30000;% £/MW

%  Operations base
C_operation_base=3000;% £/MW

% Full balance of plant (except of export cable)
C_BoP=C_array_cable+C_cable_protection+C_turbine_foundation+...
    C_offshore_substation+C_onshore_substatino+C_operation_base;
% 130000 export cable costs

%% installation and commissioning
% C_install_commiss=650000; % £/MW total costs "Component costs are below"

% Foundation installation
C_foundation_install=100000;% £/MW

% Offshore substation installation
C_offshore_substation_install=35000;% £/MW

% Onshore substation construction
C_onshore_substation_construct=25000;% £/MW

% Onshore export cable installation
C_onshore_export_cable_install=5000; % £/MW

% Offshore cable installation
% For export cable
C_offshore_export_cable_install=coef_export*220000; % £/MW 0.5

% For array cable
C_offshore_array_cable_install=(1-coef_export)*220000; % £/MW

% Final cost of installation and commisioning (except offshore export cable)
C_install_commiss=C_offshore_array_cable_install+C_onshore_export_cable_install+...
    C_onshore_substation_construct+C_offshore_substation_install+C_foundation_install;

% Operation, maintenance and service (per annum)
C_opex=75000; % £/MW

% Decommissioning
% C_decommisiong=330000;

% Turbine decommissioning
C_turbine_decomiss=45000; % £/MW
% Foundation decommissioning
C_foundation_decomiss=75000; % £/MW

% Cable decommissioning	in total
C_cable_decomiss=140000; % £/MW

% Export cable decommisionig
C_export_cable_decomiss=coef_export*C_cable_decomiss; % £/MW

% Array cables decommisionig
C_array_cable_decomiss=(1-coef_export)*C_cable_decomiss; % £/MW

% Substation decommissioning
C_substation_decomiss=65000; % £/MW

% Cost of decommisioning without export cable
C_decommisiong=C_substation_decomiss+C_array_cable_decomiss+C_foundation_decomiss+...
    C_turbine_decomiss; % £/MW

%% CAPEX and OPEX calculations
% CAPEX of export cable
CAPEX_export_cable=Overplanting_capacity(1)*(C_export_cable+...
    C_offshore_export_cable_install+C_export_cable_decomiss);  %  £/MW

% CAPEX of OWF (except export cable)
CAPEX=(C_develop_pm+C_turbine+C_BoP+C_install_commiss+C_decommisiong);%£/MW

% Total OPEX
OPEX=C_opex; % £/MW

% Convert costs in £/MW into costs in EUR
CAPEX=exchange_rate_GBP_EUR*(CAPEX*Overplanting_capacity+CAPEX_export_cable); % GPB/EUR 1.1405 in 2019
OPEX=exchange_rate_GBP_EUR*OPEX*Overplanting_capacity;
end
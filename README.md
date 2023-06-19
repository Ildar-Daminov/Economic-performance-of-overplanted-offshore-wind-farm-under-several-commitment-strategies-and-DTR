Upd (19/06/2023): The first release of MATLAB code and data.

# Economic performance of overplanted offshore wind farm under several commitment strategies and dynamic thermal ratings of submarine export cable
<img align="left" alt="Coding" width="80" src="https://ars.els-cdn.com/content/image/X03062619.jpg">

  
This repository shares the MATLAB code and data for the journal article 📋:\
Ildar Daminov, Anne Blavette, Salvy Bourguet, Hamid Ben Ahmed, Thomas Soulard, Pierre Warlop, "Economic performance of overplanted offshore wind farm under several commitment strategies and dynamic thermal ratings of submarine export cable," in Applied Energy, 2023, https://doi.org/10.1016/j.apenergy.2023.121326
  
  
## Paper's abstract
The overplanting of an offshore wind farm (OWF) with dynamic thermal ratings (DTR) is a promising solution for enhancing OWF performance. As the overplanted OWF generates an additional energy and DTR ensures its better transfer, the research question is if final OWF profits would be higher than associated costs. While there has been growing attention on this subject in recent years, there is still no research investigating the commitment strategies of overplanted OWFs with DTR. This paper investigates how commitment strategies may affect the economic performance of an overplanted OWF with DTR. The results show that, depending on the day-ahead commitment strategy, the annual revenue of OWFs may theoretically increase by up to 21% even without overplanting nor DTR, and by up to 204 % with overplanting and DTR. However, although commitment strategies, overplanting and DTR may significantly increase annual revenues of an OWF, its net present value (NPV) still heavily depends on market prices. In the presence of low market prices as in 2018, overplanting actually reduces the NPV of an OWF and, under the conditions considered here, keeps its NPV always negative. On the contrary, in the presence of high market prices as in 2022 or feed-in tariffs, the overplanting increases the OWF NPV. In the latter case, the NPV of the overplanted OWF is estimated between 1.5 billion and 9.5 billion euros while the discounted payback period is 3–12 years. The economic benefits of DTR for the overplanted OWF are estimated between 0.7 and 1 billion euros. The paper also shows that the optimal overplanting rate is very different whether the NPV or LCOE are considered. Finally, the paper shows that committing to the actual OWF power production (i.e. assuming a perfect power forecast) does not necessarily result in the highest revenue.



![2022_poster_CEL4Wind_version_1 1 - Copie](https://user-images.githubusercontent.com/73365375/229062926-c752438a-c98d-44bd-8a4a-3cdede853c8e.jpg)
<p align="center">Figure 1 Illustration of the electrical infrastructure of OWF

## How to run a code 
There are two ways how you may run this code: I. Do simulations yourself (but note that it takes few days) or II. Use precalculated data 
  
I. Launching all calculations yourself. This will reproduce the data in the paper but it would take several days:
1. Copy this repository to your computer 
2. Open the script main.m
3. Launch the script "main.m" by clicking on the button "Run" (usually located at the top of MATLAB window).\
As alternative, you may type ```main``` 
in Command Window to launch the entire script. 


II. Using the precalculated data to reproduce the particular figure : 
1. Copy this repository to your computer 
2. Open the script the Creating_Figures.m
3. Find the section (Figure XX) corresponding to the Figure you would like to reproduce. 
4. Put the cursor at any place of this section and click on the button "Run Section" (usually located at the top of MATLAB window)


## Files description

Principal scripts:
* main.m - the script which launches all calculations at a computer. Note that the entire calculates may take several days (due to neccesity to launch heavy and numerous optimization)
* Creating_Figures.m - this script repoduces the figures from the conference paper by using the precalculated data. 

Additional functions: 
* cable_thermal_model_IEC_60853_2.m - a thermal model of submarine export cable 225 kV according to the standard IEC 60853-2 (15 min is a time resolution of input data)
* losses_partial_transients.m - this script estimates losses in a export cable. Also, it calculates partial transient (p 33 IEC 60853-2) of cable and its environment  
* losses_partial_transients_new.m - the same as losses_partial_transients.m but it requires 2 additional inputs (R_AC,Lambda1_new). These two parameters are changing as a function of temperature
* temperature_correction_on_losses_new.m -this function adjusts the final profile of temperature rises as a function of losses change due to a temperature variations.
* power_curtailement.m - this function estimates power curtailements of offshore wind farm if cable's limit of temperature (90 degC) is exceeded
* I_admissible.m - This function finds the admissible current in accordance with IEC 60287-1-1, p11
* I2P.m - this function converts I (in A) to power values (in MVA) for submarine cable 225 kV
* P2I.m - this function converts power S (in MVA) to current I (in A) for submarine cable 225 kV

Initial data:
* main_simulations_case1.mat - precalculated results for the case 1 (when Static Thermal Rating is used for cable constraints)
* main_simulations_case2.mat - precalculated results for the case 2 (when Dynamic Thermal Rating is used for cable constraints)
* data_DA_IMB_years_2015 2020.mat - day-ahead and imbalance prices for 2015-2020 (all at 15-min resolution). Source: [ENTSO-E. Transparency Platform](https://transparency.entsoe.eu/)
* Elia_Fev01_2012_Jan12_2017_preloads.mat - preload power profiles of offshore wind farm (2012-2017). Source: [Elia](https://www.elia.be/en/grid-data/power-generation/wind-power-generation?csrt=6075160236430889381)
* Elia_Fev01_2012_Jul28_2021_monitored_capacity_offshore.mat - monitored installed capacity of Belgian offshore wind farm (2012-2021). Source: [Elia](https://www.elia.be/en/grid-data/power-generation/wind-power-generation?csrt=6075160236430889381)
* Elia_Jan13_2016_Jan13_2017_powers.mat - Annual power profile of Belgian offshore wind farms. Both forecast and measurements. Source: [Elia](https://www.elia.be/en/grid-data/power-generation/wind-power-generation?csrt=6075160236430889381)
* ENTSO_E_Jan13_2016_Jan13_2017_prices.mat - Annual day-ahead and imbalance prices in France according to ENTSO-e. Source: [ENTSO-E. Transparency Platform](https://transparency.entsoe.eu/)
* theta_VECTOR_ISGT_main_cases_339MW_871A.mat - Temperature rises of export cable (after current 871 A was applied at preload intervals)


SOON, this repository will contain MATLAB code and data for our [research article](https://authors.elsevier.com/c/1hDc015eif8Gcw) 

 
So far, you may see our other repositories for already published articles: 
* Conference paper: Optimal energy management of offshore wind farms considering the combination of overplanting and dynamic rating 2022, [GitHub repository](https://github.com/Ildar-Daminov/Optimal_energy_management_of_offshore_wind_farms-considering_overplanting_and_DTR_of_export_cables)
* Article: Assessment of dynamic transformer rating, considering current and temperature limitations. [GitHub repository](https://github.com/Ildar-Daminov/Assessment_Dynamic_Thermal_Rating_of_Transformers)
* Article: Demand Response Coupled with Dynamic Thermal Rating for Increased Transformer Reserve and Lifetime. [GitHub repository](https://github.com/Ildar-Daminov/Demand-response-coupled-with-DTR-of-transformers)
* Article: Energy limit of oil-immersed transformers: A concept and its application in different climate conditions. [GitHub repository](https://github.com/Ildar-Daminov/Energy-limit-of-power-transformer)
* Conference paper: Optimal ageing limit of oil-immersed transformers in flexible power systems [GitHub repository](https://github.com/Ildar-Daminov/MATLAB-code-for-CIRED-paper)
* Conference paper: Receding horizon algorithm for dynamic transformer rating and its application for real-time economic dispatch. [GitHub repository](https://github.com/Ildar-Daminov/Receding-horizon-algorithm-for-dynamic-transformer-rating)
* Conference paper: Application of dynamic transformer ratings to increase the reserve of primary substations for new load interconnection. [GitHub repository](https://github.com/Ildar-Daminov/Reserve-capacity-of-transformer-for-load-connection)

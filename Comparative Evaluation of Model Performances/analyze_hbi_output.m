
%%%%%%%%%%%%%%%%%% these script calculated the liklihoods accros subjects and creates and
%saves tables of the hbi results (tables include: Liklihoods, model
%frequencies, PEP, and estimated group parameter means
clear all; close all; clc;
%addpath('/projects/crunchie/buidze/TCG_project/cbm/cbm/codes') % add cbm toolbox to the current path only if you want to visualize the data with the cbm_hbi_plot function below
lapFile='C:\Users\user\Desktop\TCGpaper\Comparative Evaluation of Model Performances\lap_data';%name lap data file
hbiFile='C:\Users\user\Desktop\TCGpaper\Comparative Evaluation of Model Performances\hbi_data';%name hbi data file
outFile='C:\Users\user\Desktop\TCGpaper\Comparative Evaluation of Model Performances';
%% load state model and get liklihood and mean parameter values
cd(lapFile);
load("lap_state_1st_data.mat"); %import the lap() output file for the state model (if you want to run the same script for the dataset 2 change 1st to 2nd in the name.
loglik_state=mean(cbm.output.loglik);
for i=1:length(cbm.output.parameters)
    cbm.output.parameters(i,:)=transform_parameters_state_model(cbm.output.parameters(i,:));%transform the parameters for state model
end
mean_s_model=mean(cbm.output.parameters); %find the mean parameter values for state model
%% load state - movement model and get liklihood and mean parameter values
load("lap_state_movement_1st_data.mat"); %load lap() for the state movement model
loglik_sm=mean(cbm.output.loglik);

for i=1:length(cbm.output.parameters)
    cbm.output.parameters(i,:)=transform_parameters_sm_model(cbm.output.parameters(i,:)); %transform parameters for the sm model
end
mean_sm_model=mean(cbm.output.parameters);
%% load movement model and get liklihood and mean parameter values
load("lap_movement_1st_data.mat");%load lap() for the movement model
loglik_m=mean(cbm.output.loglik);
for i=1:length(cbm.output.parameters)
    cbm.output.parameters(i,:)=transform_parameters_m_model(cbm.output.parameters(i,:));%transform parameters for m model
end
mean_m_model=mean(cbm.output.parameters);

% import hbi_data
cd(hbiFile)
load("hbi_results_1st_data.mat");

%% create a table to import model statistics
% Define the data for each column
Model = {'State Movement Model'; 'State Model'; 'Movement model'};%names for different models
Loglik= [loglik_sm; loglik_state; loglik_m]; %get loglik together
Model_frequency = cbm.output.model_frequency';%model-frequencies for all the models
Exceedence_prob=cbm.exceedance.pxp';%protected ep

% Create the table
T1 = table(Model, Loglik, Model_frequency,Exceedence_prob)
% Define the path and name of the Excel file you want to create
filename = 'hbi_result_table_1st_data.xlsx';  % Update this to your desired path and file name
% Write the table to an Excel file
cd(outFile);
writetable(T1, filename);
% Confirm the file has been created
disp(['Table written to ', filename]);

%% create a table to store the model parameters
param_names = {'tau';'epsilon';'lambda';'gamma';'alpha'};
state_movement_model=mean_sm_model';
state_model=[mean_s_model(1:2)';0;0;mean_s_model(end)];
movement_model=[mean_m_model';0];

T2 = table(param_names,state_movement_model,state_model,movement_model)
% Define the path and name of the Excel file you want to create
filename = 'mean_parameters_table_1st_data.xlsx';  % Update this to your desired path and file name
% Write the table to an Excel file
writetable(T2, filename);
% Confirm the file has been created
disp(['Table written to ', filename]);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sub-functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1) transform_parameters for state model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par=transform_parameters_state_model(par)

LB = 0.0001;
UB = 0.3;
for i=1:length(par(:,1))
    par(i,1)=abs(par(i,1))+1;%tau
    par(i,2)=1/(1+exp(- par(i,2)));%epsilon
    par(i,3)=abs(par(i,3))+1;%alpha
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2) transform_parameters for state-movement model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par=transform_parameters_sm_model(par)

LB = 0.0001;
UB = 0.3;
for i=1:length(par(:,1))
    par(i,1)=abs(par(i,1))+1;%tau
    par(i,2)=1/(1+exp(- par(i,2)));%epsilon
    par(i,3) =1/(1+exp(-par(i,3)));%lambda
    par(i,3)= (UB-LB)* par(i,3) + LB;%lambda
    par(i,4)=abs(par(i,4))+1;%gamma
    par(i,5)=abs(par(i,5))+1;%alpha
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3) transform_parameters for movement model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par=transform_parameters_m_model(par)

LB = 0.0001;
UB = 0.3;
for i=1:length(par(:,1))
    par(i,1)=abs(par(i,1))+1;%tau
    par(i,2)=1/(1+exp(- par(i,2)));%epsilon
    par(i,3) =1/(1+exp(-par(i,3)));%lambda
    par(i,3)= (UB-LB)* par(i,3) + LB;%lambda
    par(i,4)=abs(par(i,4))+1;%gamma
end
end
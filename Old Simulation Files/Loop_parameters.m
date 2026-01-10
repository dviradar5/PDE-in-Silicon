% create Loop parameter table 
%
clc;clear; close all

exp_name='';
dir_exp=['results\',exp_name];
mkdir(dir_exp);

%% Loop parameters

% Material parameters
LP.material_name={'Si','GaAs'};
LP.PDE_method={'Drude_model'}; %

% Pump parameters
LP.E=[42e-9]*0.72*(0.96)^2*(0.98)^3;
LP.ND=[0.2,0.7];
LP.l=2;
LP.lambda=775e-9;
LP.z0=[-20];
LP.O_NA=0.4;

% Pump parameters
LP.pr_z0=-100;
LP.Op_NA=0.015; % Probe objective NA

% Diffusion parameters
LP.method_BC={'Neumann'};
LP.dt=25e-12;
LP.Dap=30e-4; % Diffusion cof m^2/s
LP.tau_eff=2e-9; % effective carrier lifetime
LP.nt=round([0]*1e-12/LP.dt);

%% plan simulation loops
LP_names=fieldnames(LP);
% show loops
for LP_count=1:length(LP_names)
    SL.(LP_names{LP_count})=LP.(LP_names{LP_count})(1);
end

count=1; SL(count).count=count;
for LP_count=1:length(LP_names)
    for jj=2:length(LP.(LP_names{LP_count}))
        count=count+1;
        SL(count)=SL(count-1);
        SL(count).count=count;
        SL(count).(LP_names{LP_count})=LP.(LP_names{LP_count})(jj);

    end

end

SLT=struct2table(SL);
%% save as table
writetable(SLT,[dir_exp,'SLT.xlsx']);

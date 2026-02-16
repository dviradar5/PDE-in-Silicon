clc;clear;close all
% version 0.5 11/08/2023
to_save=true;
inspect_only = false;
exp_name='N1_single';
dir_exp=['results\',exp_name];
exptype='';




%% Geometric Parameters
L0 = 1e-6;  % length unit
wl = .775;  % wavelength
dL = .02;  % grid size
dpml=10;
Lpml = dpml*dL;  % PML thickness
xn = -5-Lpml; xp = 5+Lpml;  % x boundaries
yn = -5-Lpml; yp = 20+Lpml;  % y boundaries
zn = 0; zp = dL;  % z boundaries
h = 0.00;  % metal thickness


%% Source
y_src = -1;
polarization = Axis.x;
%src = PlaneSrc(Axis.y, y_src, polarization);
src =  RectSrc(Axis.y, y_src, [zn zp; -.1 .1], polarization);

%% Loop parameters
SP=build_PDE_Loop('sizes',[xn+Lpml, xp-Lpml, 0, 0, 0, yp-Lpml],'dl',[dL,dL,dL]);

% Material parameters
SP.material_name={'Si'};


% Pump parameters
Energy=16e-9; % [J]
SP.E=8e-9; %[J], delta_E. energy/delta_e=number of pulses
SP.ND=0;
SP.z0=-7;
SP.NA=.4;
SP.update_material;
% Diffusion parameters
%LP.method_BC={'Neumann'};

SP.tau_eff=1e-9; % effective carrier lifetime

nsim=round(Energy/SP.E); %it will be the number of iterations
shape1 = Box([xn xp; yn 0; zn zp]);
SC_Material = Material(SP.material_name, 'none', SP.material_eps);
count_photons=zeros(1,nsim);
%initializing the variables:
I=0;
old_eps=0;
FCC=0;
DP=0;
save('doping_parameters','SP','I','FCC','old_eps');
%% run simulations


for ii=1:nsim

    


    %% Solution
    tic;

    [E, H, obj_array, src_array,extra] = maxwell_run_Loop(...
        'OSC', L0, wl, ...
        'DOM', SC_Material , [xn xp; yn yp; zn zp], dL, BC.p, [Lpml Lpml 0],...
        'OBJ', {'vacuum', 'k', 1.0} , shape1, ...
        'SRCJ', src, ...
        inspect_only);


    disp(['running time: ',num2str(fix(toc/60)),'min ',num2str(round(mod(toc,60))),'sec']);
    %% Ploting
    EF=abs(E{1}.array(dpml:end-dpml,1:end-dpml,2)).^2;
    disp(max(EF(:)));
    EF=EF/max(EF(:));
    xx=linspace(xn+Lpml,xp-Lpml,size(EF,1));
    yy=linspace(yn+Lpml,yp-Lpml,size(EF,2));
    eps=extra.eps{1}.array;
    if ii==1
        mind=round(size(eps,1)/2);
        ind=find(abs(eps(mind,1:end,1))>1,2,'first');
        ind=ind(2);
        I_0=create_radial_image(EF(:,ind));
    end
    EF=EF(:,ind:end)/sum(I_0(:));
    yy=yy(ind:end);
    I=EF;

    old_eps=extra.m{1};
    FCC=extra.doping;
    %figure; imagesc(EF);

    figure; subplot(321);imagesc(xx,xx,I_0); axis image
    subplot(322);    
    II=sum(I);
    plot(II/II(1)); title('normalized I');
    subplot(3,2,[3,4]);imagesc(yy,xx,I); title('I');
    colorbar;
    subplot(3,2,[5,6]);imagesc(yy,xx,FCC); title('FCC');
    colorbar;
    count_photons(ii)=sum(FCC(:));
    save('doping_parameters','SP','I','FCC','old_eps');
    status = copyfile('doping_parameters.mat', [dir_exp,'_LN',num2str(ii),'.mat']);
end



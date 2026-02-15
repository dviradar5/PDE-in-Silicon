classdef build_PDE_Loop < handle
    % build_PDE ver 18/1/2022
    properties
        %% Beam parameter
        NA=0.4 % Beam NA
        E % [J] pulse energy ;
        ND=0; %ND filter
        lambda=775e-9;   % [m] wave length
        n_index; % material complex refractive index
        z0=0; % beam waist location

        material_eps
        gauss_src=true;


        %% Geometric Parameters
        dl=[.1 .1 .1] %[dx,dy,dz] grid size
        um=1e-6  %  length unit
        Lpml = 1;  % PML thickness
        sizes;
        h;  % metal thickness

        %%  Material parameters
        material_name='Si';
        T0=300; %  room temperature (300 K)
        mu_e=1400e-4; % electron mobility m^2/(V?s)
        mu_h=450e-4; % hole mobility m^2/(V?s)
        PDE_method='drude_model'; % plasma dispersion effect model
        m_ce % the conductivity effective mass of electron
        m_ch % the conductivity effective mass of holes
        m0=9.11e-31;
        tau_D=1e-14; % Drude damping time
        charge_type
        special_prop=''; % special property: [realonly\imagonly- only the real\imaginary part of the refractive index
        
        %%  Diffusion parameters
        tau_eff=2.5e-9; % effective carrier lifetime
        tau=50e-12; % hot election life time
        dt=20e-12; %Width of each time step
        HED=false; % D(t=0) m^2/s hot election Diffusion
        nt
        Dap=2.5e-4; % Diffusion cof m^2/s
        method_BC='Dirichlet';
        Diffusion_type='ambipolar';
        
    end
    properties (Constant)
        c0 = physconst('LightSpeed');   % [m/s] speed of light
        kB=physconst('Boltzmann');
        h_bar=4.135667516*10^-15;
        eV=1.602176634e-19; % Elementary charge [C]
        eps0= 8.8541878128e-12; % Vacuum permittivity [F/m]
        e_m0=9.1093837015e-31; % electron mass
    end
    methods
        %% main
        function this = build_PDE_Loop(varargin)
            for ii=1:nargin
                if ischar(varargin{ii})
                    this.(varargin{ii})=varargin{ii+1};
                end
            end
            % load material parameters
            m=load(['material_parameters_',this.material_name]);
            this.m_ce=m.m_ce;
            this.m_ch=m.m_ch;
            [~,ind]=min(abs(this.lambda-m.wvlen*1e-9));
            this.n_index=m.n(ind)+1i*m.k(ind);
            this.material_eps=real(m.eps(ind))-1i*imag(m.eps(ind));       
        end
        %% replace properties
        function replace(this,varargin)
            for ii=1:nargin
                if ischar(varargin{ii})
                    this.(varargin{ii})=varargin{ii+1};
                end
            end
            
        end
        %% calculate number of photons
        function n_p=n_photons(this)
            n_p=this.E/this.eV/this.h_bar*this.lambda/this.c0;
        end
        function w0=waist(this)
            w0=this.lambda/this.NA/pi;
        end
        %% calculate the plasma dispersion effect (modification of the permittivity)
        function new_eps=PDE_eps(this,dNe,dNh,eps)
            if strcmp(this.PDE_method,'empirical')
                nc= -8.8e-22*dNe-8.5e-18*dNh.^0.8+1i*(8.5e-18*dNe+6e-18*dNh)/4/pi*this.lambda*100;
                if strcmp(this.special_prop,'imagonly')
                    k=imag(nc)+imag(sqrt(eps));
                    nr=real(sqrt(eps));
                elseif strcmp(this.special_prop,'realonly')
                    k=imag(sqrt(eps));
                    nr=real(nc)+real(sqrt(eps));
                else
                    k=imag(nc)+imag(sqrt(eps));
                    nr=real(nc)+real(sqrt(eps));
                end
                
                new_eps=(nr - 1i*k).^2;
            else
                afreq=2*pi*this.c0/this.lambda;
                plasma_freq=sqrt((this.eV^2/this.eps0/this.e_m0)*(dNe/this.m_ce+dNh/this.m_ch)*1e6);
                new_eps=eps-(plasma_freq/afreq).^2*(1/(1-1i/this.tau_D/afreq));
                eps_pop=-(eps-1).*(dNh+dNe)/2e23;
                new_eps=new_eps+eps_pop;
            end
        end
        %% update material properties
        function update_material(this)
            % load material parameters
            this.material_name=char(this.material_name);           
            m=load(['material_parameters_',this.material_name]);
            this.m_ce=m.m_ce;
            this.m_ch=m.m_ch;
            [~,ind]=min(abs(this.lambda-m.wvlen*1e-9));
            this.n_index=m.n(ind)+1i*m.k(ind);
            this.material_eps=real(m.eps(ind))-1i*imag(m.eps(ind));  % save as ε'-iε''     
        end
        %%
        function update_damping_time(this,cmethod)
            mcc=1/(1/this.m_ce+1/this.m_ch);
            if strcmp(cmethod,'diffusion')
                this.tau_D=this.e_m0*mcc/this.kB/this.T0*this.Dap;
            elseif strcmp(cmethod,'mobility')
                this.tau_D=this.e_m0*mcc/(1/this.mu_e+1/this.mu_h)/this.eV;
            end
        end
        %% Diffusion parameters calculation
        function D=Df(this,t)
            
            if strcmp(this.Diffusion_type,'ambipolar')
                D=this.Dap;
            else
                if strcmp(this.charge_type,'electron')
                    mu=this.mu_e;
                else
                    mu=this.mu_h;
                end
                D=mu*this.T0*this.kB/this.eV;
                if this.HED
                    D=D+mu*this.h_bar*this.c0/this.lambda/1.1*exp(-t/this.tau);
                end
            end
        end
    end
end

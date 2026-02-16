function [E_cell, H_cell, obj_array, src_array, extra] = maxwell_run_Loop(varargin)
DEFAULT_METHOD = 'direct';  % 'direct', 'gpu', 'aws', 'inputfile'

% Set solver options.
iarg = nargin; arg = varargin{iarg};
inspect_only = false;
if istypesizeof(arg, 'logical')
    inspect_only = arg;
    iarg = iarg - 1; arg = varargin{iarg};
end

is_solveropts = false;
if istypesizeof(arg, 'struct')
    solveropts = arg;
    is_solveropts = true;
    iarg = iarg - 1;  % arg = varargin{iarg};
end

if ~is_solveropts || ~isfield(solveropts, 'showstruct')
    solveropts.showstruct = true;
end

if ~is_solveropts || ~isfield(solveropts, 'method')
    solveropts.method = DEFAULT_METHOD;
end

if is_solveropts && isequal(solveropts.method, 'inputfile')
    chkarg(isfield(solveropts, 'filenamebase'), '"solveropts" should have "filenamebase" field.');
end

if ~is_solveropts || ~isfield(solveropts, 'maxit')
    % 		solveropts.maxit = intmax;
    solveropts.maxit = 1e6;
else
    chkarg(istypesizeof(solveropts.maxit, 'real') && solveropts.maxit > 0, ...
        'solveropts.maxit should be positive.');
end

if ~is_solveropts || ~isfield(solveropts, 'tol')
    solveropts.tol = 1e-6;
else
    chkarg(istypesizeof(solveropts.tol, 'real') && solveropts.tol > 0, ...
        'solveropts.tol should be positive.');
end

if ~is_solveropts || ~isfield(solveropts, 'eqtype')
    solveropts.eqtype = EquationType(FT.e, GT.prim);
else
    chkarg(istypesizeof(solveropts.eqtype, 'EquationType'), ...
        'solveropts.eqtype should be instance of EquationType.');
end

if ~is_solveropts || ~isfield(solveropts, 'pml')
    solveropts.pml = PML.sc;
else
    chkarg(istypesizeof(solveropts.pml, 'PML'), ...
        'solveropts.pml should be instance of PML.');
end

if ~is_solveropts || ~isfield(solveropts, 'returnAandb')
    solveropts.returnAandb = false;
else
    chkarg(istypesizeof(solveropts.returnAandb, 'logical'), ...
        'solveropts.returnAandb should be logical.');
end

if ~is_solveropts || ~isfield(solveropts, 'returnDiv')
    solveropts.returnDiv = false;
else
    chkarg(istypesizeof(solveropts.returnDiv, 'logical'), ...
        'solveropts.returnDiv should be logical.');
end

chkarg(iarg > 0, 'first argument is not correct.');

if inspect_only
    fprintf('%s begins (inspection only).\n', mfilename);
else
    fprintf('%s begins.\n', mfilename);
end
ge = solveropts.eqtype.ge;
fprintf('E-field grid type: %s\n', char(ge));
pm = ProgMark();

% Build the system.
% Make sure to pass the first consecutive elements of varargin to
% build_system() for correct error reports.
[osc, grid3d, s_factor, eps, mu, J, M, obj_array, src_array, mat_array, eps_node, mu_node, isiso] = ...
    build_system(ge, solveropts.pml, varargin{1:iarg}, pm);
%%

[m,J,FCC] = add_new_eps(eps,osc,src_array,grid3d,obj_array,J);
extra.doping=FCC;
eps=m;

%%
if ~is_solveropts || ~isfield(solveropts, 'F0')
    solveropts.F0 = 'zero';  % 'rand' is the other choice
else
    chkarg(isequal(solveropts.F0, 'zero') || isequal(solveropts.F0, 'rand') ...
        || istypesizeof(solveropts.F0, 'complexcell', [1 Axis.count], grid3d.N), ...
        'solveropts.F0 should be length-%d cell array whose each element is %d-by-%d-by-%d array with complex numbers.', ...
        Axis.count, grid3d.Ncell{:});
end

if inspect_only  % inspect objects and sources
    if solveropts.showstruct
        figure;
        set(gcf, 'units','normalized','position',[0.5 0 0.5 0.5]);
        withpml = true;
        visobjsrc(grid3d, obj_array, src_array, withpml);
        drawnow
        pm.mark('domain visualization');
    end
    
    % Visualize modes.
    is_modalsrc = false;
    for src = src_array
        if istypesizeof(src, 'ModalSrc')
            is_modalsrc = true;
            modalsrc = src;
            opts.withabs = true;
            
            E2d = modalsrc.E2d;
            cmax = max(abs([E2d{Axis.x}.array(:); E2d{Axis.y}.array(:); E2d{Axis.z}.array(:)]));
            opts.cmax = cmax;
            for w = Axis.elems
                figure;
                set(gcf, 'units','normalized','position',[(int(w)-1)/3 1/2 1/3 1/3]);
                vis2d(E2d{w}, obj_array, opts);
                drawnow;
            end
            
            H2d = modalsrc.H2d;
            cmax = max(abs([H2d{Axis.x}.array(:); H2d{Axis.y}.array(:); H2d{Axis.z}.array(:)]));
            opts.cmax = cmax;
            for w = Axis.elems
                figure;
                set(gcf, 'units','normalized','position',[(int(w)-1)/3 0 1/3 1/3]);
                vis2d(H2d{w}, obj_array, opts);
                drawnow;
            end
        end
    end
    
    if is_modalsrc
        pm.mark('distributed source visualization');
    end
    fprintf('%s finishes (inspection only).\n\n', mfilename);
    E = {};
    H = {};
else  % inspect_only == false
    if isequal(solveropts.method, 'inputfile')
        % Define eps_node_array at vertices of the E-field edges.
        chkarg(isiso, 'anisotropic materials are not supported for solveropts.method == ''inputfile''.');
        
        if ge == GT.prim
            eps_node = eps_node{Axis.x};  % ad hoc solution
            
            % Old implementation (without anisotropy support) starts
            % here.
            eps_node_array = eps_node.data_expanded();  % (Nx+2) x (Ny+2) x (Nz+2)
            eps_node_array = eps_node_array(1:end-1, 1:end-1, 1:end-1);  % (Nx+1) x (Ny+1) x (Nz+1)
            
            % The below line is an ad hoc solution for the error.  It
            % is just to pass an array with correct size to
            % write_input(), but eps_node_array is not used in
            % write_input() currently.
            eps_node_array = eps_node_array(1:end-1, 1:end-1, 1:end-1);  % Nx x Ny x Nz
            %
            % 				% (Anisotropy support? begins)
            % % 				eps_node_cell = eps_cell;
            % % 				for w = Axis.elems
            % % 					ind_g = {':',':',':'};
            % % 					if grid3d.bc(w) == BC.p
            % % 						ind_g{w} = grid3d.N(w);
            % % 					else
            % % 						ind_g{w} = 1;
            % % 					end
            % % 					eps_node_cell{w} = cat(int(w), eps_node_cell{w}(ind_g{:}), eps_node_cell{w});
            % % 				end
            % % 				eps_node_cell{Axis.x} = 2./(1./eps_node_cell{Axis.x}(1:end-1, :, :) + 1./eps_node_cell{Axis.x}(2:end, :, :));
            % % 				eps_node_cell{Axis.y} = 2./(1./eps_node_cell{Axis.y}(:, 1:end-1, :) + 1./eps_node_cell{Axis.y}(:, 2:end, :));
            % % 				eps_node_cell{Axis.z} = 2./(1./eps_node_cell{Axis.z}(:, :, 1:end-1) + 1./eps_node_cell{Axis.z}(:, :, 2:end));
            % %
            % % 				eps_node_array = (eps_node_cell{Axis.x} + eps_node_cell{Axis.y} + eps_node_cell{Axis.z})./3;
            % %
            % % 				eps_node_array = (eps_node_array(1:end-1,1:end-1,1:end-1) + eps_node_array(2:end,1:end-1,1:end-1) ...
            % % 							+ eps_node_array(1:end-1,2:end,1:end-1) + eps_node_array(1:end-1,1:end-1,2:end) ...
            % % 							+ eps_node_array(1:end-1,2:end,2:end) + eps_node_array(2:end,1:end-1,2:end) ...
            % % 							+ eps_node_array(2:end,2:end,1:end-1) + eps_node_array(2:end,2:end,2:end))./8;
            % 				% (Anisotropy support? ends)
            %
            % 				eps_node_array = 8./(1./eps_node_array(1:end-1,1:end-1,1:end-1) + 1./eps_node_array(2:end,1:end-1,1:end-1) ...
            % 							+ 1./eps_node_array(1:end-1,2:end,1:end-1) + 1./eps_node_array(1:end-1,1:end-1,2:end) ...
            % 							+ 1./eps_node_array(1:end-1,2:end,2:end) + 1./eps_node_array(2:end,1:end-1,2:end) ...
            % 							+ 1./eps_node_array(2:end,2:end,1:end-1) + 1./eps_node_array(2:end,2:end,2:end));
        else
            eps_node_array = eps_node.data_original();  % Nx x Ny x Nz
        end
        
        write_input(solveropts.filenamebase, solveropts.eqtype, osc, grid3d, s_factor, ...
            eps_node_array, eps, mu, J, M, solveropts.F0, solveropts.tol, solveropts.maxit);
        
        pm.mark('input file creation');
        fprintf('%s finishes. (input file created)\n\n', mfilename);
        E = {};
        H = {};
    else  % solveropts.method ~= 'inputfile'
        if isequal(solveropts.method, 'direct')
            [E, H] = solve_eq_direct(solveropts.eqtype, solveropts.pml, osc.in_omega0(), eps, mu, s_factor, J, M, grid3d);
        elseif isequal(solveropts.method, 'iterative')
            [E, H, relres, iter, resvec] = solve_eq_iterative(solveropts.maxit, solveropts.tol, solveropts.F0, solveropts.eqtype, solveropts.pml, osc.in_omega0(), eps, mu, s_factor, J, M, grid3d);
        else  % for solvers using E-field based equation
            %% Apply spatial inversion.
            % 				d_prim = grid3d.dl(:, GT.prim);
            % 				d_dual = grid3d.dl(:, GT.dual);
            % 				s_prim = s_factor(:, GT.prim);
            % 				s_dual = s_factor(:, GT.dual);
            d_prim = flip_vec(grid3d.dl(:, GT.dual));  % GT.dual, not GT.prim
            d_dual = flip_vec(grid3d.dl(:, GT.prim));  % GT.prim, not GT.dual
            s_prim = flip_vec(s_factor(:, GT.dual));  % GT.dual, not GT.prim
            s_dual = flip_vec(s_factor(:, GT.prim));  % GT.prim, not GT.dual
            mu = flip_vec(mu);
            eps = flip_vec(eps);
            J = neg_vec(flip_vec(J));  % pseudovector
            
            if isequal(solveropts.method, 'gpu')
                ds_prim = mult_vec(d_prim, s_prim);
                ds_dual = mult_vec(d_dual, s_dual);
                figure;
                F0 = {zeros(grid3d.N), zeros(grid3d.N), zeros(grid3d.N)};
                [E, H] = fds(osc.in_omega0(), ...
                    ds_prim, ds_dual, ...
                    mu, eps, ...
                    F0, J, ...
                    solveropts.maxit, solveropts.tol, 'plot');
                %   norm(A2 * ((1./e) .* (A1 * y)) - omega^2 * m .* y - A2 * (b ./ (-i*omega*e))) / norm(b) % Error for H-field wave equation.
            elseif isequal(solveropts.method, 'aws')
                ds_prim = mult_vec(d_prim, s_prim);
                ds_dual = mult_vec(d_dual, s_dual);
                F0 = {zeros(grid3d.N), zeros(grid3d.N), zeros(grid3d.N)};
                % 					F0 = {rand(grid3d.N), rand(grid3d.N), rand(grid3d.N)};
                % 					F0 = {rand(1)*ones(grid3d.N), rand(1)*ones(grid3d.N), rand(1)*ones(grid3d.N)};
                callback = maxwell(osc.in_omega0(), ...
                    ds_prim, ds_dual, ...
                    mu, eps, ...
                    F0, J, ...
                    solveropts.maxit, solveropts.tol);
                while ~callback(); end
                [~, E, H] = callback();
            end
            
            E = neg_vec(flip_vec(E));  % pseudovector
            J = neg_vec(flip_vec(J));  % pseudovector
            H = flip_vec(H);
        end
        
        pm.mark('solution calculation');
        if isequal(solveropts.method, 'iterative')
            %% Report the iterative solution process.
            semilogy(1:length(resvec), resvec);
            fprintf('\t# of iteration steps = %d\n', iter);
            fprintf('\trelative residual error = %e\n', relres);
        end
        fprintf('solve for: %s-field\n', char(solveropts.eqtype.f));
        fprintf('%s finishes.\n\n', mfilename);
    end
end

if solveropts.returnAandb
    %		[A, b, g_from_f] = create_eq(solveropts.eqtype, solveropts.pml, osc.in_omega0(), eps, mu, s_factor, J, M, grid3d);
    eq = MatrixEquation(solveropts.eqtype, solveropts.pml, osc.in_omega0(), eps, mu, s_factor, J, M, grid3d);
    
    [extra.A, extra.b] = eq.matrix_op();
    [~, ~, extra.GfromF] = eq.matrixfree_op();
end

if solveropts.returnDiv
    [Dive, Divm] = create_divs(ge, s_factor, grid3d);
    extra.Dive = Dive;
    extra.Divm = Divm;
end

% Construct Scalar3d objects.
E_cell = cell(1, Axis.count);
H_cell = cell(1, Axis.count);
J_cell = cell(1, Axis.count);
M_cell = cell(1, Axis.count);
for w = Axis.elems
    if ~isempty(E)
        E_cell{w} = array2scalar(E{w}, grid3d, ge, FT.e, w, osc, PhysQ.E);
    end
    
    if ~isempty(H)
        H_cell{w} = array2scalar(H{w}, grid3d, ge, FT.h, w, osc, PhysQ.H);
    end
    
    J_cell{w} = array2scalar(J{w}, grid3d, ge, FT.e, w, osc, PhysQ.J);
    M_cell{w} = array2scalar(M{w}, grid3d, ge, FT.h, w, osc, PhysQ.M);
end

extra.grid3d = grid3d;
extra.J = J_cell;
extra.M = M_cell;
extra.eps = eps_node;
extra.mu = mu_node;
extra.m=m;
end

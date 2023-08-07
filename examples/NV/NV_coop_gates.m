function [dyn, U1, U2] = NV_coop_gates(d_extra)
% Example: Concurrent optimization method for co-operative gates.
%
% Assume we wish to perform a Ramsey sequence:
% Start in the |0> state, do a pi/2_y rotation to the Bloch sphere
% equator, wait for a time t while the state gains a relative phase,
% do a -pi/2_y rotation, measure |0><0|.
% If the drift Hamiltonian is H = \hbar \omega Z/2, the measured
% signal should be s(t) = (1+cos(\omega t))/2.
%
% For the sequence to work, the rotations do not have to be exact y rotations.
% It's enough if they "co-operate" by bringing the state to the
% equator and back with phases that cancel each other. This extra
% degree of freedom can make gates shorter and/or easier to optimize.
%
% We can accomplish this by optimizing a state transfer sequence
% from |0> to |0>, but in the middle of the sequence project the state to the equator.
% Unless the first half of the sequence has already brought the
% state to the equator, it will lose some purity and thus the
% optimization can no longer reach zero error.
% The result of this optimization is a sequence where the first and
% second halves represent the two co-operative pi/2 rotations.
%
%! M. Braun and S. Glaser, New J. Phys. 16, 115002 (2014).

% Ville Bergholm 2015-2016


% We may do the propagation either in Hilbert or Liouville space.
% The Liouville space solution is much slower for large systems.
use_liouville_space = 0;
% In Hilbert space we may either use a vec expansion or a lrmul
% construction to perform the projection.
use_vec_expansion = 1;

%% define experiment parameters

D = 2.871; %% ZFS GHz
transition = 2.512040802; % Transitoin fre GHz
drive = 2.522040802;
detuning = (drive-transition);%*1e-3;

rabi_fre_x = 5e-3; % in Ghzgamma*B_field
rabi_fre_y = 5e-3;

rabi_fre = max([rabi_fre_x,rabi_fre_y]);


mult = 1;
pi_pul = mult*0.5*(1/(rabi_fre)); %% ns;
n_bin = 2; %ns time bin


gamma_nv = 28; % GHz/T

B_x = 0; B_y = 0; B_z = abs(D-transition)/gamma_nv; % Gauss
% 
% 
% h_zfs = D*(Z^2 - 2/3); %% canceled in the rotating frame
% 
% H_b = gamma_nv*(B_x*X + B_y*Y + B_z*Z); %% canceled in the rotating frame


%% dimensions

if nargin < 1
    d_extra = 2;  % dimension of the extra subsystem
end
dim = [2, d_extra];  % one qubit plus something extra

D = prod(dim);

%% operators

% Z = diag([1, -1]); % Pauli Z
p0 = diag([1, 0]); % |0><0|
p1 = diag([0, 1]); % |1><1|
S = [1,0,0,0; 0,0,1,0; 0,1,0,0; 0,0,0,1]; % 2-qubit SWAP gate

%% Liouvillian superoperators

% projector to the equator of the Bloch sphere, a Hermitianness-preserving mapping.
P_xy = [1, 0, 0, 1; 0, 2, 0, 0; 0, 0, 2, 0; 1, 0, 0, 1] / 2;
% Corresponding projector for the adjoint system
P_xy_adj = S * P_xy.' * S;

% extend the projectors to the other subsystem
P_xy = extend_P(P_xy, d_extra);
P_xy_adj = extend_P(P_xy_adj, d_extra);


%% corresponding lrmul Hilbert space operators

A = {diag([0.5,1]), diag([1,0.5]), [0,0;0.5,0], [0,0.5;0,0]};
B = {p0, p1, [0,1;0,0], [0,0;1,0]};

% extend the lrmul matrices to the other subsystem
I = eye(d_extra);
for k=1:length(A)
    A{k} = kron(A{k}, I);
    B{k} = kron(B{k}, I);
end

%% Set up Dynamo

ini = zeros(D, 1);
ini(1) = 1/sqrt(2);
ini(3) = 1/sqrt(2);
fin = ini;%zeros(D, 1); fin(end) = 1;
% random, diagonal, traceless H_drift
%temp = randn(D, 1);
%H_drift = diag(temp-sum(temp)/D); %Z / 2;

H_drift = detuning+diag([1,-1,1,-1]);
[H_ctrl, c_labels] = control_ops(dim, 'xy');


% limit driving Rabi frequency to the interval [-1, 1]
control_type = 'mmmm';
par_lim_x = ((1/sqrt(2))*rabi_fre_x*(2*pi))*[-1, 2];
par_lim_y = ((1/sqrt(2))*rabi_fre_y*(2*pi))*[-1, 2]; %% MHz %%%%%%%%%%%%%%%%%% w/2pi(1/sqrt(2))*
control_par = {par_lim_x, par_lim_y, par_lim_x, par_lim_y};


T = ceil(pi_pul); % pulse duration

if use_liouville_space
    task = 'open state overlap'
else
    task = 'closed state'
    fin = fin * fin';
end

% random initial controls

n_int = round((T)/n_bin);

tau_fac = 0; % set to 1 to have variable time bin
% % Rabi Freq limit
x_con = zeros(1,n_int) + 0.01*2*pi*rabi_fre*(1/mult);
y_con = x_con;

InitialControlles=[x_con' y_con' x_con' y_con'];

dyn = dynamo(task, ini, fin, H_drift, H_ctrl);
dyn.system.set_labels('Co-operative gates optimization', dim, c_labels);
dyn.seq_init(n_int, [T, tau_fac], control_type, control_par);

% random, constant initial controls
dyn.set_controls(InitialControlles);

%% Set up the projection

% replace the propagator in the middle bin by P_xy
t = ceil(n_int/2);
if use_liouville_space
    %temp = dyn.set_state_transform(t, 1, P_xy);
    temp = dyn.set_state_transform(t, 1, @proj_liouville);
else
    temp = dyn.set_state_transform(t, 1, @proj_hilbert);
end


options = struct(...
    'error_goal',        0.5 * (1e-5)^2 / dyn.system.norm2,...
    'max_evals',         1000,...
    'max_walltime',      180000,...
    'max_cputime',       10e6,...
    'min_gradient_norm', 1e-7,...
    'plot_interval',     1, ...
    'OptimalityTolerance',1e-7);

mask = dyn.full_mask(true);



mask = mask & temp;  % avoid clobbering it during update


%% Optimize and test the result

dyn.ui_open();
dyn.search(mask,options);


% extract the resulting pair of co-operating gates
U1 = dyn.cache.combined_propagator(1, t-1);
U2 = dyn.cache.combined_propagator(t+1, n_int);


return


    function ret = proj_liouville(bin, rho, adjoint)
        % Apply the projection on vectorized density matrices in Liouville space.
        
        if adjoint
            ret = rho * P_xy;
        else
            ret = P_xy * rho;
        end
    end


    function ret = proj_hilbert(bin, rho, adjoint)
        % Apply the projection on density matrices in Hilbert space, either by
        % (0) vectorizing the density matrix, applying the map in
        %     Liouville space, and then de-vectorizing the result, or
        % (1) using a lrmul-type construction.
        
        if use_vec_expansion
            if adjoint
                temp = P_xy_adj;
            else
                temp = P_xy;
            end
            ret = inv_vec(temp * vec(rho));
        else
            ret = 0;
            for k=1:length(A)
                if adjoint
                    ret = ret +B{k} * rho * A{k};
                else
                    ret = ret +A{k} * rho * B{k};
                end
            end
        end
    end
end

function P = extend_P(P, d_extra)
% Extends the given 1-qubit Liouville space operator P to act as
% identity for the other subsystem.

temp = lmap(P, {[2,2], [2,2]});
temp = gate.two(temp, [1,3], [2,d_extra,2,d_extra]);
%P = full(temp.data);
P = temp.data;
end
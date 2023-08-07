% Example 3 for DPG 2012 Stuttgart

% Simulating exciton transport on a spin chain with Markovian noise
% and Hamiltonian controls.


randseed(78318);


SX = [0 1; 1 0];
SY = [0 -1i; 1i 0];
SZ = [1 0; 0 -1];
I = eye(2);
SP = (SX +1i*SY)/2;


%% the physics of the problem

% Exciton transport chain with the fiducial sink as the last node,
% sucking probability from the given chain site.

n_sites = 3; % actual chain sites
dim = 2 * ones(1, n_sites+1); % dimension vector: qubits, or spin-1/2 chain
D = prod(dim); % total system dimension
target_site = n_sites

% 0: ground, 1: excited state, SP lowers/annihilates
a = SP;
n_op = a' * a; % == (I-SZ) / 2; % number op


%% parameters

% energy splittings
omega = [ones(1, n_sites), 0];
omega(2) = 4;
% site-to-site coupling
v = 0.8
% dephasing
gamma = zeros(1, n_sites)
% relaxation
Gamma = 1e-2 * ones(1, n_sites)
% transfer rate from target_site to sink
transfer_rate = 6


desc = sprintf('%d-qubit exciton transport chain with XX+YY interaction, Z controls, relaxation.', n_sites);
% 'transfer rate = %g, split = %g, v = %g', n_sites,transfer_rate,omega(2),v);
fprintf('%s\n\n', desc);



%% subspace restriction

% The Hamiltonian preserves exciton number, whereas the noise
% processes we use here leave the 0-or-1 exciton manifold
% invariant. Hence we shall limit ourselves to it:

ddd = n_sites + 2; % zero- and single-exciton subspaces (sink included)
p = [1, 1+2.^(n_sites:-1:0)]; % states to keep: zero, single exciton at each site/sink
q = setdiff(1:D, p); % states to throw away

st_labels = {'loss', 'site 1', 'site 2*', 'site 3', 'sink'};


%% Lindblad ops

diss = cell(1, n_sites);
deph = cell(1, n_sites);
for k = 1:n_sites
    diss{k} = op_list({{sqrt(2 * Gamma(k)) * a,   k}}, dim);
    diss{k} = diss{k}(p,p);

    deph{k} = op_list({{sqrt(2 * gamma(k)) * n_op, k}}, dim);
    deph{k} = deph{k}(p,p);
end
sink = {op_list({{a, target_site; sqrt(2 * transfer_rate) * a', n_sites+1}}, dim)};
sink{1} = sink{1}(p,p);


%% drift Hamiltonian

J = v * 2 * [1 1 0]; % XY coupling
C = diag([ones(1,n_sites-1), 0], 1); % linear chain

H_drift = heisenberg(dim, @(s,a,b) J(s)*C(a,b))...
  +op_sum(dim, @(k) omega(k) * n_op);

H_drift = H_drift(p,p);


%% drift Liouvillian (noise / dissipation)

L_drift = superop_lindblad(sink, H_drift) +superop_lindblad(diss) +superop_lindblad(deph);


%% controls

% Control Hamiltonians
[H_ctrl, c_labels] = control_ops(dim, '1:3z')
for k = 1:n_sites
    H_ctrl{k} = full(H_ctrl{k}(p,p));
end

% transformed controls?
control_type = char('m' + zeros(1, length(H_ctrl)));
pp = [-1, 2] * 8;
control_par = {pp, pp, pp};


%% initial and final states

% for pure state transfer
initial = zeros(D,1); initial(prod(dim(2:end))+1) = 1; % '10..0'
final   = zeros(D,1); final(2) = 1;  % '0..01'

dyn = dynamo('open state overlap', initial(p), final(p), L_drift, H_ctrl);
dyn.system.set_labels(desc, st_labels, c_labels);


%% set up controls
T = 10;
dyn.seq_init(151, T * [0.5, 1.0], control_type, control_par);
%dyn.set_controls(0.1);
dyn.set_controls(0.05*[1, 2, 3]);
%dyn.set_controls([0, 2.6, 0]);


%% now do the actual search

dyn.ui_open();
dyn.search();
%dyn.analyze();

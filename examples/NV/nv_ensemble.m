clear all; close all;

%% Units for the variable 
% Frequency - GHz
% Time - ns
% B field - Gauss

%% Teenie meenie details

pulse_type = 'pi_pulse'; %%'pi_pulse'; % %'pi_pulse'; %'pi_half_pulse'; %% % 'three_pi_half_pulse' %%
NV_name = 'NV1_PP1_Y';
time = clock;
data_file = ['Pulses\' pulse_type '_X_' NV_name '_dyn_data_' date '_' num2str(time(4)) '_' num2str(time(5)) ];


%% define experiment parameters
tran = '0_1';
D = 2.871; %% ZFS GHz
transition = 3.28062914; % Transitoin fre GHz
drive = 3.28062914;
detuning = drive-transition;

rabi_fre_x = 5e-3; % in Ghzgamma*B_field
rabi_fre_y = 5e-3;
rabi_fre = max([rabi_fre_x,rabi_fre_y]);

rabi_fre_x_lim = 7.9e-3; % in Ghzgamma*B_field
rabi_fre_y_lim = 9.9e-3;



mult = 3;
pi_pul = 0.5*(1/(rabi_fre)); %% ns;
n_bin = 2; %ns time bin


%% define matices and other parameters , units GHz---ns
switch tran
    case '0_1'
        X=0.5*[0 1;1 0];
        Y=0.5*[0,-1i; 1i 0];
        Z =[1, 0; 0,-1];
        Z2 = [1, 0; 0,0];
        H_det=2*pi*detuning*Z;
    case '0_min1'
        X=1/2*[0 1;1 0];
        Y=1/2*[0,1i; -1i 0];
        Z =[0, 0; 0,-1];
        Z2=[0, 0; 0,1];
        H_det=2*pi*[detuning,            0;...
            0,  -detuning];
end


gamma_nv = 28; % GHz/T

B_x = 0; B_y = 0; B_z = abs(D-transition)/gamma_nv; % Gauss


h_zfs = D*(Z^2 - 2/3); %% canceled in the rotating frame

H_b = gamma_nv*(B_x*X + B_y*Y + B_z*Z); %% canceled in the rotating frame

H_drift =  H_det;


%% Define the physics of the problem

% dimension vector for the quantum system
n_qubits = 1;
dim = 2 * ones(1, n_qubits);
D = prod(dim);

% Control Hamiltonians / Liouvillians
[H_ctrl, c_labels] = control_ops(dim, 'xy');
n_controls = length(H_ctrl);

% limit driving Rabi frequency to the interval [-1, 1]
control_type = 'mm';
par_lim_x = ((1/sqrt(2))*rabi_fre_x*(2*pi))*[-1, 2];
par_lim_y = ((1/sqrt(2))*rabi_fre_y*(2*pi))*[-1, 2]; %% MHz %%%%%%%%%%%%%%%%%% w/2pi(1/sqrt(2))*
control_par = {par_lim_x, par_lim_y};

shake_par = 0.1*max([par_lim_x,par_lim_y]);

switch pulse_type
    case 'pi_pulse'
        % closed gate optmisation
                initial = [1 0;0 1];
                Sy=[0,-1i;1i,0];
                Sx=[0, 1; 1,0];
                final = Sy;%1/sqrt(2)*[1 1]';
        % %         state optmisation
        %         initial = [1 0;0 0];
        %         final = [0 0;0 1];
        % ket optimisation
%         initial = [1, 0]';
%         final = [0, 1]';
        T = mult*pi_pul; %% ns
    case 'pi_half_pulse'
        
        % closed gate optmisation
                initial = [1 0;0 1];
                Had = (1/sqrt(2))*[1, 1; 1, -1];
                S_xh= (1/sqrt(2))*[1, -1j; -1j, 1];
                S_yh= (1/sqrt(2))*[1, -1; 1, 1];
                final = S_xh;%1/sqrt(2)*[1 1]';
        %
        %                 initial = [1 0;0 0]';
        %                 final = (1/sqrt(2))*[0 1;1 0]';
        % % State transfer
%         initial = [1, 0]';
%         final = (1/sqrt(2))*[1, 1]';
        
        T = mult*0.5*pi_pul; %% ns
    case 'three_pi_half_pulse'
        
        % closed gate optmisation
        initial = [1 0;0 1];
        S_x3h= -(1/sqrt(2))*[1, 1j; 1j, 1];
        S_y3h= -(1/sqrt(2))*[1, 1; -1, 1];
        final = S_x3h;%1/sqrt(2)*[1 1]';
        
        %         initial = [1 0;0 0]';
        %         final = (1/sqrt(2))*[0 1;1 0]';
        T = mult*1.5*pi_pul; %% ns
end


%% define the ensemble

no_of_det = 3;
no_of_amp_var = 1e2;

detun_limit = 2.1e-3; %GHz, additive
amp_limit = 0.5; %% multiplicative

detun = linspace(-detun_limit,detun_limit,no_of_det); %[8e-3,10e-3,12e-3];%
amp_var = linspace(1-amp_limit,1+amp_limit,no_of_amp_var); 

%detun = detun_limit*rand(1,no_of_det) - detun_limit; 
% amp_var = (1-amp_limit)*(rand(1,no_of_amp_var)); 

ens_det = kron(ones(1,length(amp_var)),detun);
ens_amp = [];
amp = zeros(1,length(detun));

for elee = 1:1:length(amp_var)
    for ele = 1:length(amp_var):length(ens_det)
        ens_amp = [ens_amp,amp_var(elee)];
    end
end

% set(0,'DefaultFigureWindowStyle','normal')

x = ens_det;
y = ens_amp;
A_func = @(x)   x*Z + H_drift;  % vary drift Guess:: additional term to the drift
B_func = @(y,c) (1-y) * H_ctrl{c};

%distribution of weights for the ensemble

%sigma = 0.5;
%mu = 1;%0.5*(1+amp_limit);1;%
%weight =  exp(-(y-mu).^2 / (2*sigma^2));% 
weight = ones(1,length(y));%


%weight_det = normrnd(1,0.2,[1,reps]);
%weight_amp = normrnd(0.65,0.2,[1,reps]);
weight = weight/max(weight);

f = figure;
subplot(3,2,1);
plot(y ,'-*')
hold on;
plot(weight);
title('Ensemble Amplitude')
subplot(3,2,2);
plot(x,'-x');
title('Ensemble detuning');


task_='open gate'; % other options

% 'closed gate': Target operation: unitary gate (ignoring global phase) in a closed system.
% 'open gate': quantum gate in an open system under Markovian noise.
% 'closed state': Target operation: mixed state transfer in a closed system.

dyn = dynamo(task_, initial, final, @(k) A_func(x(k)), @(k,c) B_func(x(k),c), weight, n_controls);
dyn.system.set_labels('NV ensemble optimization', dim, c_labels);


%% Initial controls

% random initial controls
tau_fac = 0;

n_int = round((T)/n_bin);
dyn.seq_init(n_int, [T, tau_fac], control_type, control_par); %% first parameter is time slots




%set controls

% sinusoidal
% x_con = (1/(sqrt(2)))*rabi_fre*sin(1:1:n_int);
% y_con = (1/(sqrt(2)))*rabi_fre*cos(1:1:n_int);

%Gaussian
% sigma = n_int/6;
% gaus= @(x,sigma) exp(-((x-0.5*n_int).^2)./(2*sigma^2));

% for i = 1:1:n_int
% x_con(i) = (1/(sqrt(2)))*rabi_fre*gaus(i,sigma);
% y_con(i) = -(1/(sqrt(2)))*rabi_fre*gaus(i,sigma);
% end


% % Rabi Freq limit
x_con = zeros(1,n_int) + 0.1*2*pi*rabi_fre*(1/mult);
y_con = x_con;


% % random initial pulse
% InitialControlles = (2*rand(1, 2)-1);

InitialControlles=[x_con' y_con'];

dyn.set_controls(InitialControlles);
dyn.shake(shake_par,shake_par,false);
mask = dyn.full_mask(true);


pause(2);

%% Now do the actual search


dyn.ui_open();

options = struct(...
    'error_goal',        0.5 * (1e-5)^2 / dyn.system.norm2,...
    'max_evals',         1000,...
    'max_walltime',      180000,...
    'max_cputime',       10e6,...
    'min_gradient_norm', 1e-6,...
    'plot_interval',     1, ...
    'OptimalityTolerance',1e-6);
dyn.search(mask,options);

% figure;
% for i = 1:1:length(x)
%     subplot(round(length(x)/4),4,i);
%     dyn.plot_pop(1,2);
%     str = ['amp var: ' num2str(ens_amp(i)) '| det: ' num2str(ens_det(i))];
%     title(str);
%     drawnow;
% end
%
% %

dyn.analyze;

% subplot(3,2,[3 4]); 
figure;
% dyn.plot_pop();

set(0, 'currentfigure', f);
%% plot control field

X = dyn.seq.fields(:,1);
Y = dyn.seq.fields(:,2);

amp = sqrt(X.^2 + Y.^2);
phi = atan(Y./X);
t = 1:1:length(X);
t = n_bin*t';
theta = t*2*pi*transition + phi;
sig = amp.*sin(theta);
subplot(3,2,[5 6]);
plot(t,sig,'b',t,amp,'r','LineWidth',2)
title('Control field');
data.per = T; data.x = X;data.y = Y; ...
    data.rabi_freq_x = rabi_fre_x; data.rabi_freq_y = rabi_fre_y; ...
    data.dt = n_bin; data.transition =  transition; ...
    data.drive = drive;data.detuning = detuning; data.task = task_; ...
    data.ens_amp = ens_amp;data.ens_det = ens_det;
savename = 'Complete_data_NV1_TH_Ex';
save(savename);
save(data_file,'data')

%% Plot optimal error
fig_0 = figure;

axes_0 = axes('Parent',fig_0);

semilogy(dyn.stats{1, 1}.error,'LineWidth',1, 'Color', [1 0 0], 'LineStyle','-','Linewidth',3)
set(axes_0,'FontSize',30,'FontName','CMU Serif','Linewidth',3)
% title('FoM Evolution','FontSize',30,'FontName','CMU Serif')
ylabel('Optimization error','FontSize',30,'FontName','CMU Serif')
xlabel('Algorithm iteration no.','FontSize',30,'FontName','CMU Serif')



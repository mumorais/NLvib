clearvars;
% close all;
clc;
addpath('../../SRC');
%% Geometry and physical Parameters of the TLCD
epsilon = 3.542;       % Headloss [.]
g = 9.788 * 1000;      % Gravity [mm/s^2]
b = (224);             % Horizontal length [mm]
h = 55;                % vertical length [mm]
L = (2*h + b);         % Total length [mm]
%% Dynamic Parameters of TLCD
omega_w = (2*g/L)^(1/2); % Natural frequency TLCD [rad/s]
gamma = 1.0; % tuning ratio [gamma = omega_w / omega_s]
%% Parameters of the Structure
omega_s = omega_w/gamma; % Natural frequency Ideal Structure [rad/s]
xi_s = 0.02; % Damping coefficient [.]
%% Dimensionless Parameters
alpha = b/L; % Aspect Ratio 
mu = 0.10; % Mass Ratio [mu = M_w / M_s] 
%% Matrices
M = [1+mu alpha*mu; 
    alpha*mu mu];
D = [2*omega_s*xi_s 0; 
    0 0];
K = [omega_s^2 0; 
    0 mu*omega_w^2];
%% Elementos não-lineares
D_fluido = mu*(1/(2*L))*epsilon;

nonlinear_elements{1} = struct('type','quadraticdamper',...
    'damping',D_fluido,'force_direction',[0;1]);

% nonlinear_elements{2} = struct('type','cubicSpring',...
%     'stiffness',0.01*omega_s^2,'force_direction',[1;0]);

%% Harmonic excitation acting on the structure
f0 = 100; % Forcing coefficient [mm/s^2]
Fex1 = [f0; 0];
%% Defining the system as chain of oscillators
n = size(M,1); % number of degrees of freedom

oscillator.M = M; % Acceleration coefficient
oscillator.D = D; % Velocity coefficient
oscillator.K = K; % Displacement coefficient
oscillator.Fex1 = Fex1; % Excitation force
oscillator.nonlinear_elements = nonlinear_elements; % Elements NL
oscillator.epsilon = epsilon; % Headloss 
oscillator.L = L; % Total length
oscillator.n = n; % Horizontal length
%% Compute frequency response using harmonic balance
% Analysis parameters
H = 5;          % harmonic order
N = 4*H+1;      % number of time samples per period, 
           % cf. Appendix A in https://doi.org/10.1016/j.ymssp.2019.106503
Om_s = 0.5;     % start frequency
Om_e = 70.00;   % end frequency

% Initial guess (solution of underlying linear system)
% (linear response to initiate Newton HB's method) 
Q1 = (-Om_s^2*oscillator.M + 1i*Om_s*oscillator.D + oscillator.K)\Fex1; 
y0 = zeros((2*H+1)*length(Q1),1); 
y0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)]; 

% Solve and continue w.r.t. Om
ds = .01;
Sopt = struct('jac','none'); 
Sopt.stepmax = 1500000;        % Set the maximum number of steps


[X_HB] = solve_and_continue(y0,...
    @(X) HB_residual_TLCD_2GDL(X,oscillator,H,N),...
    Om_s,Om_e,ds,Sopt);

%% Interpret solver output
Om_HB = X_HB(end,:); % Freq. Excitation
Q_HB = X_HB(1:end-1,:); % Coef. of Fourier

% Define amplitude as magnitude of the fundamental harmonic of the second
% coordinate's displacement

% Amplitude considering only the first harmonic
w_HB = sqrt(Q_HB(n+2,:).^2 + Q_HB(2*n+2,:).^2); 
u_HB = sqrt(Q_HB(n+1,:).^2 + Q_HB(2*n+1,:).^2);

% Amplitude considering the RMS of the sum of harm. up to the truncation 
% order.
u_rms_HB = sqrt(sum(Q_HB(1:2:end,:).^2))/sqrt(2);
w_rms_HB = sqrt(sum(Q_HB(2:2:end,:).^2))/sqrt(2);

%% Illustrate results

fig_TLCD_Est_2GDL_NLin = figure('Name', 'Respostas RMS', 'Position', [100, 100, 800, 600]);
tiledlayout(2, 1);

% --- Resposta da Estrutura (U) ---
ax1 = nexttile;
plot(Om_HB, u_rms_HB, 'r-', 'LineWidth', 2); hold on;
set(ax1, 'xlim', [Om_s Om_e], 'TickLabelInterpreter', 'latex', 'FontSize', 12);
ylabel('Amplitude da Estrutura $U_{rms}$ [mm]', 'Interpreter', 'latex');
title(sprintf('FRF RMS da Estrutura ($f_0 = %.2f$ $mm/s^2$)', f0), 'Interpreter', 'latex');
%legend('Não-linear (HB)', 'Linearizado', 'Location', 'best', 'Interpreter', 'latex');
grid on; 

% --- Resposta do TLCD (W) ---
ax2 = nexttile;
plot(Om_HB, w_rms_HB, 'b-', 'LineWidth', 2); hold on;
set(ax2, 'xlim', [Om_s Om_e], 'TickLabelInterpreter', 'latex', 'FontSize', 12);
xlabel('Frequência de excitação Ω [rad/s]', 'Interpreter', 'latex');
ylabel('Amplitude do TLCD $W_{rms}$ [mm]', 'Interpreter', 'latex');
title(sprintf('FRF RMS do TLCD ($f_0 = %.2f$ $mm/s^2$)', f0), 'Interpreter', 'latex');
%legend('Não-linear (HB)', 'Linearizado', 'Location', 'best', 'Interpreter', 'latex');
grid on; 

% %% Plot results
% figure; hold on;
% plot(Om_HB,u_HB,'g-');
% plot(Om_HB,w_HB,'b-');
% set(gca,'xlim',[Om_s Om_e]);
% xlabel('excitation frequency'); ylabel('response amplitude');
% legend('u','w');
% grid on;
% 
% 
% figure; hold on;
% plot(Om_HB,u_rms_HB,'g-');
% plot(Om_HB,w_rms_HB,'b-');
% set(gca,'xlim',[Om_s Om_e]);
% xlabel('excitation frequency'); ylabel('response amplitude');
% legend('u_{rms}','w_{rms}');
% grid on;
% 
% %figure; hold on;
% %plot(Om_HB,u_HB,'g-');
% %plot(Om_HB,w_HB,'b-');
% %plot(Om_HB,u_rms_HB,'r-');
% %plot(Om_HB,w_rms_HB);
% %set(gca,'xlim',[Om_s Om_e]);
% %xlabel('excitation frequency'); ylabel('response amplitude');
% %legend('u','w','u_{rms}','w_{rms}');
% 
% %% Maxima amplitude da estrutura 
% 
% [max_amplitude, indice_pico] = max(u_rms_HB);
% 
% freq_ressonancia = Om_HB(indice_pico);
% 
% % printar o resultado
% fprintf('A amplitude máxima encontrada é: %.4f mm\n', max_amplitude);
% fprintf('A frequência de ressonância é: %.4f rad/s\n', freq_ressonancia);
% 
% % Gráfico marcando o pico 
% figure; hold on;
% plot(Om_HB, u_rms_HB, 'b-');
% plot(freq_ressonancia, max_amplitude, 'ro', 'MarkerFaceColor', 'r'); % Marca o pico
% xlabel('Frequência (rad/s)');
% ylabel('Amplitude');
% title('FRF com Pico Identificado');
% legend('u_{rms}','Pico');
% grid on;

%% Illustrate results

% FRF of the isolated structure
path1 = 'C:\Users\Guilherme\OneDrive\Área de Trabalho\Mapa de Resposta Lin e NLin\FRF Estrutura 1GDL\Coef. de Força = 100.00';

if isfolder(path1)
    Amp_struct = load(fullfile(path1, 'a_RMS.mat'));
    Om_struct  = load(fullfile(path1, 'Om.mat'));

    Amp_1GDL = Amp_struct.a_RMS;
    Om_1GDL  = Om_struct.Om;
end

fig_TLCD_Est_2GDL_NLin_comp = figure('Name', 'Respostas RMS', 'Position', [100, 100, 800, 600]);
tiledlayout(2, 1);

% --- Resposta da Estrutura (U) ---
ax1 = nexttile;
plot(Om_HB, u_rms_HB, 'r-', 'LineWidth', 1, 'DisplayName', 'Com TLCD'); hold on;
plot(Om_1GDL, Amp_1GDL, 'k-', 'LineWidth', 1, 'DisplayName', 'Sem TLCD');
set(ax1, 'xlim', [Om_s 15], 'TickLabelInterpreter', 'latex', 'FontSize', 12);
ylabel('Amplitude da Estrutura $U_{rms}$ [mm]', 'Interpreter', 'latex');
title(sprintf('FRF RMS da Estrutura ($f_0 = %.2f$ $mm/s^2$)', f0), 'Interpreter', 'latex');
legend('show', 'Location', 'best');
grid on; 

% --- Resposta do TLCD (W) ---
ax2 = nexttile;
plot(Om_HB, w_rms_HB, 'b-', 'LineWidth', 1); hold on;
set(ax2, 'xlim', [Om_s 15], 'TickLabelInterpreter', 'latex', 'FontSize', 12);
xlabel('Frequência de excitação Ω [rad/s]', 'Interpreter', 'latex');
ylabel('Amplitude do TLCD $W_{rms}$ [mm]', 'Interpreter', 'latex');
title(sprintf('FRF RMS do TLCD ($f_0 = %.2f$ $mm/s^2$)', f0), 'Interpreter', 'latex');
%legend('Não-linear (HB)', 'Linearizado', 'Location', 'best', 'Interpreter', 'latex');
grid on; 

%% Salvando arquivos 
pasta_base = 'C:\Users\Guilherme\OneDrive\Área de Trabalho\Mapa de Resposta Lin e NLin';

% Salvando Para Amortecimento Não-Linear 
if length(nonlinear_elements) == 1
    Forcing = f0;

    pasta_TLCD_Est_2GDL_NLin = fullfile(pasta_base, 'FRF TLCD_Est_2GDL_NLin', ...
            sprintf('Coef. de Força = %.2f _ Ordem de Truncamento H = %.0f', Forcing, H));

    % Criar se não existir
    if ~exist(pasta_TLCD_Est_2GDL_NLin, 'dir')
        mkdir(pasta_TLCD_Est_2GDL_NLin);
    end

    % Arquivos que serão salvos
    save(fullfile(pasta_TLCD_Est_2GDL_NLin, 'u_rms_HB.mat'), 'u_rms_HB');
    save(fullfile(pasta_TLCD_Est_2GDL_NLin, 'w_rms_HB.mat'), 'w_rms_HB');
    save(fullfile(pasta_TLCD_Est_2GDL_NLin, 'Om_HB.mat'), 'Om_HB');

    % Figura que serão salvos 
    exportgraphics(fig_TLCD_Est_2GDL_NLin, fullfile(pasta_TLCD_Est_2GDL_NLin, 'FRF do TLCD_Estrutura 2GDL.png'));
    exportgraphics(fig_TLCD_Est_2GDL_NLin_comp, fullfile(pasta_TLCD_Est_2GDL_NLin, 'FRF do TLCD_Estrutura 2GDL (comparando sem TLCD).png'));
    
% Com a Mola Cúbica
else
    Forcing = f0;

    pasta_TLCD_Est_2GDL_NLin = fullfile(pasta_base, 'FRF TLCD_Est_2GDL_NLin (Com Mola Cúbica)', ...
            sprintf('Coef. de Força = %.2f', Forcing));

    % Criar se não existir
    if ~exist(pasta_TLCD_Est_2GDL_NLin, 'dir')
        mkdir(pasta_TLCD_Est_2GDL_NLin);
    end

    % Arquivos que serão salvos
    save(fullfile(pasta_TLCD_Est_2GDL_NLin, 'u_rms_HB.mat'), 'u_rms_HB');
    save(fullfile(pasta_TLCD_Est_2GDL_NLin, 'w_rms_HB.mat'), 'w_rms_HB');
    save(fullfile(pasta_TLCD_Est_2GDL_NLin, 'Om_HB.mat'), 'Om_HB');

    % Figura que serão salvos 
    exportgraphics(fig_TLCD_Est_2GDL_NLin, fullfile(pasta_TLCD_Est_2GDL_NLin, 'FRF do TLCD_Estrutura 2GDL com Mola Cúbica.png'), ...
        'Position', [100, 100, 800, 600]);
    exportgraphics(fig_TLCD_Est_2GDL_NLin_comp, fullfile(pasta_TLCD_Est_2GDL_NLin, 'FRF do TLCD_Estrutura 2GDL (comparando sem TLCD).png'), ...
        'Position', [100, 100, 800, 600]);
end

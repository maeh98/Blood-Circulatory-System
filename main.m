%{
--------------------------- MAIN FILE -------------------------------------
With this project we have modelled and simulated the circulatory system
using 14 compartments: [large arteries, arterioles, capillaries, venules, 
large veins, left atrium, left ventricle, large arteries, arterioles capillaries, 
venules, large veins, right atrium, right ventricle].

You can choose between four numerical methods -- see below.
You can choose between 2 different implementations of the adaptive time
step, with a tradeoff between speed and accuracy -- see below.
We have non-differentiable points in the model. To address that you can
choose whether slightly modify the model, or use a root-finding method --
see below.

METHOD_num has options @RKF45, @RK_Bogacki_Shampine, @RK_Cash_Karp,
@RK_Dormand_Prince.
METHOD_adapative has options @adaptive_deprecated (faster but less accurate) and
@adaptive (slower but more accurate).
METHOD_diff has options @F (with root_finding) and @F_modified (using
modified RHS).
%}

% clear all
clc
format long

METHOD_num      = @RKF45;
METHOD_adaptive = @adaptive;
METHOD_diff     = @F;
tf = 4.0;   % Set finishing time

% Set all the parameters
T = 1.0;                            % Duration of a heartbeat
t_c = T*[0.8; 0; 0.8; 0];           % Start contraction, heart chambers [LA,LV,RA,RV]
T_cp = T*[0.17; 0.3; 0.17; 0.3];    % Duration of contraction, heart chambers [LA,LV,RA,RV]
t_r = T*[0.97; 0.3; 0.97; 0.3];     % Start relaxation, heart chambers [LA,LV,RA,RV]
T_rp = T*[0.17; 0.15; 0.17; 0.15];  % Duration of relaxation, heart chambers [LA,LV,RA,RV]

vol_d = [50; 30; 53; 75; 75; 30; 10; 580; 300; 600; 1900; 50; 30; 10];     % Dead volumes
L = [0.00002; 0.00002; 0.00002; 0.00002; 0.00002; 0.0002; 0.00001; 0.0005; 0.006; 0.002; 0.001; 0.001; 0.0002; 0.00001];% Inertial parameters
R = [0.0250; 0.03; 0.050; 0.03; 0.025; 0.001; 0.003; 1.025; 0.23; 0.03; 0.0075; 0.0003; 0.001; 0.003];   % Viscous parameters
C = [2.62; 1.98; 2.4750; 6.5; 4.6; 0.000016; 0.000025; 1.36; 0.04; 2; 90; 20; 0.000016; 0.000025];       % Compliance
B = [1.6e-5; 2.5e-5; 1.6e-5; 2.5e-5];   % Turbolence parameters, heart chambers [LA,LV,RA,RV]
Cd = [0.01; 0.01; 0.01; 0.01];          % Dynamic characteristsics, valves [MV,AV,TV,PV]
Cv = [2.50; 2.5; 2; 2.5];               % Ejection velocity, valves [MV,AV,TV,PV]
E_A = [0.26; 2.75; 0.17; 0.6];          % Amplitude of the elasticity function, heart chambers [LA,LV,RA,RV]
E_B = [0.3; 0.08; 0.22; 0.05];          % Basis value of the elasticiy function, heart chambers [LA,LV,RA,RV]

%%EE

% Initialize
x_new = 1.0e+03 * [0.0862;    0.0590;    0.0907;    0.1756;    0.1463;    0.0581;    0.1939;    0.6043;    0.3004;    0.6136;    2.4910;    0.1955;    0.0448;    0.1407;   -0.0333;   -0.0186;   -0.0050;   -0.0014;    0.0115;    0.0072;   -0.00001;    0.0089;    0.0091;    0.0116;    0.0071;   -0.0323;   -0.0001;   -0.00001;    0.0009;    0.00001;    0.0001;    0.00001];  
P = 1.0e+03 * [0.0138;    0.0147;    0.0152;    0.0155;    0.0155;    0.0152;    0.0148;    0.0179;    0.0088;    0.0068;    0.0066;    0.0073;  0.0055;    0.0065] ;
dt = 1e-4;  % First step
t0 = 0.0;   % Initial time

% Pre-allocation
lengt = 2e3;    % Initial sizes of the arrays
t_array = zeros(lengt, 1);      % Contains time values
sol = zeros(32, lengt);         % Contains solutions at all time values
P_saved = zeros(14, lengt);     % Contains pressures at all time values
E_saved = zeros(4, lengt);      % Contains elasticity values (of heart chambers) at all time values
dt_hist = zeros(lengt, 1);      % Contains (next) time step sizes for all time values
adap_step_c = zeros(lengt, 1);  % Records for all time values the number of times a step had to be repeated due to accuracy
bjoerck_c = zeros(lengt, 1);    % Records for all time values the number of iterations needed for root-findind with Anderson-bjoerck

% Set variables for first loop
i = 1;
t = t0;
dt_new = 1e-4;
counter = 0;
crossed = 0;
% sol(:, 1) = x_new;
% P_saved(:, 1) = P;
% dt_hist(1) = dt;
tic

while t <= tf
    if mod(i, 1) == 0%  % saving every 10/20 loops smooths out the plot (but doesn't in the actual solution.)
        i_save = i/1;%  %save all of them
        sol(:, i_save) = x_new;
        P_saved(:, i_save) = P;
        dt_hist(i_save)=dt;
        t_array(i_save) = t;
    end
    
    % Update variables
    x_old = x_new;  
    i = i+1;
    t_old = t;
    t = t+dt;
    
    % Perform a step of size dt_new with the chosen method
    dt = dt_new;
    [x_new, dt, dt_new, adap_step_count, bjoerck_count, crossed] = METHOD_adaptive(x_old, t_old, P, dt, L, R, B, Cd, Cv, crossed, METHOD_num, METHOD_diff);
    
    % Update counter
    adap_step_c(i-1) = adap_step_count;
    bjoerck_c(i-1) = bjoerck_count;

    % Postprocessing on valves
    x_new(20:21) = x_new(20:21).* x_new(29:30);          
    x_new(27:28) = x_new(27:28).* x_new(31:32);
    
    % Evaluation of the elasticity functions
    E = zeros(4, 1);
    E(1) = E_A(1).*e_a(t, T_cp(1), T_rp(1), t_c(1), t_r(1), T) + E_B(1);
    E(2) = E_A(2).*e_v(t, T_cp(2), T_rp(2), T) + E_B(2);
    E(3) = E_A(3).*e_a(t, T_cp(3), T_rp(3), t_c(3), t_r(3), T) + E_B(3);
    E(4) = E_A(4).*e_v(t, T_cp(4), T_rp(4), T) + E_B(4);
    E_saved(:, i)= E;

    % Pressure update
    P(1:5) = (1./C(1:5)) .* (x_new(1:5) - vol_d(1:5));
    P(6:7) = (E(1:2).*(x_new(6:7) - vol_d(6:7)))./(1-0.0005*(x_new(19:20)-x_new(20:21)));
    P(8:12) = (1./C(8:12)) .* (x_new(8:12) - vol_d(8:12));
    P(13:14) = (E(3:4).*(x_new(13:14) - vol_d(13:14)))./(1-0.0005*(x_new(26:27)-x_new(27:28)));    
    
    % Preallocation issues
    if length(t_array)-2 <= i
            t_array = [t_array; zeros(size(t_array))];
            sol = cat(2, sol, zeros(size(sol)));
            E_saved = cat(2, E_saved, zeros(size(E_saved)));
            dt_hist = cat(1, dt_hist, zeros(size(dt_hist)));  
    end
    
    
end       
loop = toc  

% Resize all arrays with the stored solutions approprietely
k  = size(nonzeros(t_array), 1)-1;      
E_saved = E_saved(:, 1:k);
t_array = t_array(1:k);
dt_hist = dt_hist(1:k);
sol = sol(:, 1:k);
P_saved = P_saved(:, 1:k);
adap_step_c = adap_step_c(1:k);
bjoerck_c = nonzeros(bjoerck_c(1:k));  % Only considering when we do use the Anderson-Bjoerck method.

dt_mean = mean(dt_hist) % Was our adaptive timestep good?
bjoerck_mean = mean(bjoerck_c)
adap_step_c = mean(adap_step_c)


%%%%%%%%%%%%%% End Code --- Start plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
figure(1)
plot(t_array, sol(6:7, :),t_array, sol(13:14, :))
title('Volumes heart chambers')
xlabel('t')
ylabel('V')
legend('V_{LA}', 'V_{LV}','V_{RA}', 'V_{RV}')


%
figure(20)
plot(t_array.', dt_hist.', 'LineWidth',2.0)
title('Adaptive Time-step size')
xlabel('time [s]')
ylabel('Step size dt')
set(gca,'FontSize',15)
% xlim([5,7])

figure(2)
plot(t_array, sol(20:21, :), t_array, sol(27:28, :),'LineWidth',2.0)
title('Fluxes heart chambers')
xlabel('t')
ylabel('Q')
legend('Q_{MV}', 'Q_{AV}', 'Q_{TV}', 'Q_{PV}')
set(gca,'FontSize',15)
% xlim([5,7])
% 
%

% figure(3)
% plot(t_array, sol(29:32, :))
% title('Opening coefficients')
% xlabel('t')
% ylabel('O')
% ylim([-0.5,1.5])
% legend('O_{MV}', 'O_{AV}', 'O_{TV}', 'O_{PV}')
% %
% figure(4)
% plot(t_array, sol(1:5, :), t_array, sol(8:12, :))
% title('Volumes outside heart')
% xlabel('t')
% ylabel('V')
% legend('V_{1}', 'V_2', 'V_3', 'V_4', 'V_5', 'V_8', 'V_9', 'V_{10}', 'V_{11}', 'V_{12}')
% 
% figure(5)
% plot(t_array, sol(15:19, :), t_array, sol(22:26, :))
% title('Fluxes outside heart')
% xlabel('t')
% ylabel('Q')
% legend('Q_{1}', 'Q_{2}', 'Q_{3}', 'Q_{4}', 'Q_{5}', 'Q_{8}', 'Q_{9}', 'Q_{10}', 'Q_{11}', 'Q_{12}')
% % 
% figure(7)
% plot(t_array, P_saved(1:5, :), t_array, P_saved(8:12, :))
% title('Pressure outside heart')
% xlabel('t')
% ylabel('P')
% legend('P_{1}', 'P_2', 'P_3', 'P_4', 'P_5', 'P_8', 'P_9', 'P_{10}', 'P_{11}', 'P_{12}')
% 
% x = linspace(0,3);
% y1 = sin(5*x);
% y2 = sin(15*x);

figure(17)
plot(t_array, P_saved(6:7, :), t_array, P_saved(13:14, :))
title('Pressure inside heart')
xlabel('t')
ylabel('P')
legend('P_{6}', 'P_7', 'P_13', 'P_14')
% figure(20)
% tiledlayout(1,2)
% 
% %Top plot
% ax1 = nexttile;
% plot(ax1,t_array, sol(20:21, :), t_array, sol(27:28, :))
% title(ax1,'Top Plot')
% xlim(ax1, [11,14]);
% ylabel(ax1,'flow rate')
% 
% 
% % Bottom plot
% ax2 = nexttile;
% plot(ax2,t_array, P_saved(6:7, :), t_array, P_saved(13:14, :))
% title(ax2,'Top Plot')
% ylabel(ax2,'pressure')
% xlim(ax2, [11, 14]);

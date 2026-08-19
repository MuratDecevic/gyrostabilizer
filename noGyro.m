% Roll motion of a boat without gyrostabilizer
% Townsend and Shenoi

%% boat
length = 4; % m
beam = 1.3; % m
m = 500; % displacement [kg]
k = 0.4 * beam; % radius of gyration [m]
w4 = 2.64; % nondimensional natural roll frequency
beta_d = 0.3; % damping ratio
%% system parameters
I44 = m*k^2; % system roll inertia
A44 = 0.3*I44; % added mass
C44 = (w4^2)*(I44+A44); % roll restoring
B44 = 2*beta_d*sqrt(C44*(I44+A44)); % roll damping
Isys = I44+A44;
%% time
tspan=[0 100];
%% initial conditions
% x0 = [theta;theta']
x0 = [0;0];
%% solve
p=[Isys,B44,C44];
[t,x] = ode45(@dynamics, tspan, x0,[], p);

theta=rad2deg(x(:,1));
thetap=x(:,2);
beta=rad2deg(x(:,3));
betap=x(:,4);
%% plot
plot(t,theta,'-','LineWidth',1.5)
xlabel('t')
ylabel('deg')
legend('roll')
grid on


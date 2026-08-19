% A gyroscopic wave energy recovery system for marine vessels

%% boat
length = 4; % m
beam = 1.3; % m
m = 500; % displacement [kg]
k = 0.4 * beam; % radius of gyration [m]
w4 = 2.64; % nondimensional natural roll frequency
beta_d = 0.3; % damping ratio
%% flywheel
mg = 6.5; % mass of flywheel [kg]
mm = 4.3; % mass of motor [kg]
ro = 0.5; % outer radius of flywheel [m]
ri = 0.3; % inner radius of flywheel [m]
thickness = 0.03; % axial thickness of flywheel [m]
spinRate = 6000; % angle ksi [r/min]
ws = (spinRate*2*pi)/60; % rad/s
%Is= 0.5*(mg+mm)*rg^2; %kgm^2 solid disk
Is = 0.8*mg*ro^2; %kgm^2 rim-weighted rotor
Im = 0.1; % transversal moment of inertia of a motor, approximated as a sylinder radius=0.075m
Itr = (mg/12)*(3*(ro^2+ri^2)+thickness^2); %kgm^2
It = Itr+Im; % transversal moment of inertia of flywheel + motor
%% system parameters
I44 = m*k^2; % system roll inertia
A44 = 0.3*I44; % added mass
C44 = (w4^2)*(I44+A44); % roll restoring
B44 = 2*beta_d*sqrt(C44*(I44+A44)); % roll damping
Isys = I44+A44;
% values of parameters tkaen from "A gyroscopic wave energy recovery system
% for marine vessels" by Townsend and Shenoi
bg = [10 20 50 100 500 1000 2000]; % Nms
cg = [40 60 80 100 500 1000 2000]; % Nm
[B,C]=ndgrid(bg,cg);
b_c_combinations=[B(:),C(:)];
%% time
tspan=[0 100];

%% initial conditions
% x0 = [theta;theta';beta;beta']
x0 = [0;0;0;0];

%% load previous solutions

filename = 'solutions.mat';

if isfile(filename)
    load(filename, 'solutions');
    fprintf('Loaded %d previously calculated solutions.\n', ...
        numel(solutions));
else
    solutions=struct( ...
        'b', {}, ...
        'c', {}, ...
        't', {}, ...
        'theta', {}, ...
        'beta', {} ...
        );
    fprintf('No previous solutions found. Starting from the beginning.\n');
end

%% solve for new combinations
for i=1:size(b_c_combinations,1)
    Bg=b_c_combinations(i,1);
    Cg=b_c_combinations(i,2);
    
    fprintf('\nChecking Bg=%g, Cg=%g ...\n',Bg,Cg);
    
    % Check whether this combination was used
    if isempty(solutions)
        idx=[];
    else
        idx=find(...
            [solutions.b]==Bg & ...
            [solutions.c]==Cg, ...
            1 ...
            );
    end
    
    %If solution exists, don't solve again
    if ~isempty(idx)
        fprintf('Already calculated.\n')
        continue;
    end
    
    % Solve
    p=[Isys,B44,C44,Is,ws,It,Bg,Cg];
    [t,x] = ode45(@dynamics, tspan, x0,[], p);
    
    % Store the solution
    new_index=numel(solutions)+1;
    solutions(new_index).b=Bg;
    solutions(new_index).c=Cg;
    solutions(new_index).t=t;
    solutions(new_index).theta=rad2deg(x(:,1));
    solutions(new_index).beta=rad2deg(x(:,3));
    
    % Save
    save(filename,'solutions');
    fprintf('Solution saved.\n');
    
end
% theta=rad2deg(x(:,1));
% thetap=x(:,2);
% beta=rad2deg(x(:,3));
% betap=x(:,4);

%% plot
for k=1:numel(solutions)
    figure;
    plot(solutions(k).t,solutions(k).theta,'-', ...
        solutions(k).t,solutions(k).beta,'--', ...
        'LineWidth', 1.5)
    xlabel('t')
    ylabel('deg')
    title(sprintf(...
        'Bg=%g, Cg=%g',solutions(k).b, solutions(k).c));
    legend('roll','precession')
    grid on
end
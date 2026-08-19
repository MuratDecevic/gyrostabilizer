function Mwave = waveForcing(t)
% ============================================================
% WAVEFORCING
%
% JONSWAP wave-induced roll moment for a 4 m fishing boat
%
% Input:
%     t      time [s]
%
% Output:
%     Mwave  wave excitation moment [N*m]
%
% Boat:
%     L = 4.0 m
%     B = 1.3 m
%     displacement = 500 kg
%
% Sea state:
%     Hs = 2.0 m
%     Tp = 11.0 s
%     h  = 50 m
%     gamma = 3.3
%
% 50 random wave components
% ============================================================

%% ------------------------------------------------------------
% Sea-state parameters
% ------------------------------------------------------------

g     = 9.81;
Hs    = 2.0;       % significant wave height [m]
Tp    = 11.0;      % peak period [s]
h     = 50.0;      % water depth [m]
gamma = 3.3;       % JONSWAP peak enhancement

N = 50;

fmin = 0.04;       % minimum frequency [Hz]
fmax = 0.30;       % maximum frequency [Hz]

%% ------------------------------------------------------------
% Wave excitation coefficient
% ------------------------------------------------------------

% Converts wave slope into roll excitation moment.
%
% Units: N*m
%
% This is an assumed value for a 4 m fishing boat.

Kw = 5000;         % [N*m]

%% ------------------------------------------------------------
% Frequency discretization
% ------------------------------------------------------------

f = linspace(fmin,fmax,N);

df = f(2)-f(1);

omega = 2*pi*f;

fp = 1/Tp;

%% ------------------------------------------------------------
% JONSWAP spectrum
% ------------------------------------------------------------

sigma = 0.09*ones(1,N);

sigma(f <= fp) = 0.07;

r = exp( ...
    -(f-fp).^2 ./ ...
    (2*sigma.^2*fp^2) );

% Unscaled JONSWAP spectrum

Sshape = ...
    (g^2/(2*pi)^4) .* ...
    f.^(-5) .* ...
    exp(-1.25*(fp./f).^4) .* ...
    gamma.^r;

%% ------------------------------------------------------------
% Normalize spectrum to Hs = 2 m
% ------------------------------------------------------------

m0_target = Hs^2/16;

alpha = m0_target/trapz(f,Sshape);

S = alpha*Sshape;

%% ------------------------------------------------------------
% Wave amplitudes
% ------------------------------------------------------------

a = sqrt(2*S*df);

%% ------------------------------------------------------------
% Random phases
% ------------------------------------------------------------

% Fixed seed -> same wave realization every simulation

persistent phi

if isempty(phi)

    rng(20260812);

    phi = 2*pi*rand(1,N);

end

%% ------------------------------------------------------------
% Finite-depth wave numbers
%
% omega^2 = g*k*tanh(k*h)
% ------------------------------------------------------------

k = zeros(1,N);

for n = 1:N

    % Initial guess
    kn = omega(n)^2/g;

    % Newton-Raphson
    for iter = 1:50

        kh = kn*h;

        F = g*kn*tanh(kh) - omega(n)^2;

        dF = g*tanh(kh) + ...
             g*kn*h*sech(kh)^2;

        kn_new = kn - F/dF;

        if abs(kn_new-kn) < 1e-12
            break;
        end

        kn = kn_new;

    end

    k(n) = kn_new;

end

%% ------------------------------------------------------------
% Wave-induced roll moment
%
% Wave slope:
%
%     theta_wave = a*k*cos(...)
%
% Roll forcing:
%
%     Mwave = Kw * theta_wave
%
% ------------------------------------------------------------

Mwave = 0;

for n = 1:N

    Mwave = Mwave + ...
        Kw*a(n)*k(n)* ...
        cos(omega(n)*t + phi(n));

end

end
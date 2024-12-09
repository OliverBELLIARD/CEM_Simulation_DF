function scriptFDTD05(time, alpha)
%% Calcul FDTD (Finite Difference in Time Domain) - Magic-time step
% Valeurs par défaut en cas de non arguments
if nargin < 1
    time = 1500;
end
if nargin < 2
    alpha = 1;
end

%% Définition des constantes
eps0 = 8.854e-12;  % Permittivité du vide en F/m
epsr = 4;
mu0 = 4*pi*1e-7;

%% Maillage : discretisation spatiale
L = 0.5; % longueur du domaine de calcul
max_space = 501; % nb de points spatiaux (nb de champ E)
dz = L/(max_space-1);

%% Discretisation temporelle;
max_time = time / alpha;
dt = alpha * sqrt(eps0*mu0)*dz;

%% Initialisation des champs E et H
% dimension en H = dimension en E - 1 du au schema
% conditions aux limites imposees ici pour E(1) et E(max_space)
E = zeros(max_space, 1);
H = zeros(max_space-1, 1);

%% Conditions aux limites
% Conditions sur E
E(1) = 0;
E(max_space) = 0;

%% Constantes
gamma = -1/eps0 * dt/dz;
tau = -1/mu0 * dt/dz;

%% Conditions absorbantes
temp(1) = 0;
temp(2) = 0;

%% source
start = 2;
t0 = 400*dt;
spread = 1.6e-11;

%% Discretisation suivant z
for k=1:max_space
    zE(k) = (k-1)*dz;
end

%% Diélectrique
dielec_deb = 200;
dielec_fin = 300;

for u=1:max_space
  if (u>=dielec_deb & u<=dielec_fin)
    alphaEdielec(u) = gamma./epsr;
  else
    alphaEdielec(u) = gamma;
  end
end

%% Boucle temporelle
figure
for n=1:max_time
    t = (n-1) * dt;

    % Equation : calcul du champ electrique
    for k = 2:max_space - 1
        E(k) = E(k) + alphaEdielec(k)* (H(k) - H(k-1));
    end

    % Soft source
    pulse = exp(-1*((t-t0)/spread)^2);
    E(start) = E(start) + pulse;

    % Absorption aux limites
    E(1) = temp(1);
    temp(1) = E(2);
    E(max_space) = temp(2);
    temp(2) = E(max_space-1);

    % Calcul du champ magnétique
    for k = 1:max_space-1
        H(k) = H(k) + tau * (E(k+1) - E(k));
    end

    % Visualisation des champs
    hold off;
    area(zE, alphaEdielec-min(alphaEdielec), 'faceColor','c')
    hold on;
    plot(zE, E)
    title(["Champ E", "\nalpha = " num2str(alpha) ", max time = " num2str(max_time)])
    ylabel("E [V/m]")
    xlabel("z (position dans l'espace) [m]")
    axis([0 L -1.1 1.1])
    pause(0.0001)
end

end

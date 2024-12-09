%% Calcul FDTD (Finite Difference in Time Domain) - Ondes planes
function scriptFDTD02(time, alpha)
% Valeurs par défaut en cas de non arguments
if nargin < 1
    time = 100;
end
if nargin < 2
    alpha = 1;
end

%% Définition des constantes
eps0 = 8.854e-12;  % Permittivité du vide en F/m
mu0 = 4*pi*1e-7;

%% Maillage : discretisation spatiale
L = 2; % longueur du domaine de calcul
max_space = 201; % nb de points spatiaux (nb de champ E)
dz = L/(max_space-1);

%% Discretisation temporelle;
%alpha = 1;
max_time = time / alpha;
dt = alpha * sqrt(eps0*mu0)*dz;

%source
center = 101;
t0 = 40*dt;
spread = 1.6e-10;

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

figure
for n=1:max_time
    t = (n-1) * dt;
    
    % Equation : calcul du champ electrique
    for k = 2:max_space - 1
        E(k) = E(k) + gamma * (H(k) - H(k-1));
    end
    
    % Soft source
    %if n<=60
    pulse = exp(-1*((t-t0)/spread)^2);
    E(center) = E(center) + pulse;
    %end
    % Calcul du champ magnétique
    for k = 1:max_space-1
        H(k) = H(k) + tau * (E(k+1) - E(k));
    end

    %% Visualisation des champs
    subplot(1, 2, 1);
    plot(E)
    title("Champ E", "alpha="+alpha+", time="+time)
    ylabel("E [V/m]")
    xlabel("t")
    
    subplot(1, 2, 2);
    plot(H)
    title("Champ H")
    ylabel("H [T]")
    xlabel("t")
    pause(0.1)
end

end
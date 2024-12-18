function scriptFDTD03(time, alpha)
%% Calcul FDTD (Finite Difference in Time Domain) - Ondes planes dans l'espace
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
max_time = 1000;
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

%% Source spatiale
spread = 1.6e-10;
z0 = 1; % position du centre du domaine
c0 = 1./sqrt(eps0*mu0); % vitesse de propagation
vecz = (1:max_space-2)*dz;
gamma2 = (vecz-z0)./c0;
pulse2 = exp(-1*(gamma2./spread).^2);
E(2:max_space-1) = pulse2;

%% Discretisation suivant z
for k=1:max_space
    zE(k) = (k-1)*dz;
end

%% Boucle temporelle
figure
for n=1:max_time
    t = (n-1) * dt;
    
    % Equation : calcul du champ electrique
    for k = 2:max_space - 1
        E(k) = E(k) + gamma * (H(k) - H(k-1));
    end

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
    plot(zE, E)
    title("Champ E", "alpha="+alpha+", time="+time+", t="+t+" s")
    ylabel("E [V/m]")
    xlabel("z (position dans l'espace) [m]")
    axis([0 L -1.1 1.1])
    pause(0.05)
end

end
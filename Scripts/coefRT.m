function [R,T] = coefRT(epsr0, epsr1)
%% Permet de calculer les constantes R et T
% comme spécifié précédemment en fixant les valeurs des permittivités relatives
R = (sqrt(epsr0) - sqrt(epsr1))/(sqrt(epsr0) + sqrt(epsr1))
T = (2 * sqrt(epsr0))/(sqrt(epsr0) + sqrt(epsr1))
end
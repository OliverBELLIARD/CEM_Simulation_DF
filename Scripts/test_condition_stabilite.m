clear; close all; clc

alphas = [0.5 0.95 0.99 1 1.01 1.1];

for i=1:length(alphas)
    scriptFDTD01(100, alphas(i))
end
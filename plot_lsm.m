% plot_lsm.m
clear; clc; close all;

filename = 'image_lsm.txt';

if ~exist(filename, 'file')
    error('Le fichier %s n''existe pas. Lancez le main C++ d''abord.', filename);
end

% Lecture des données
fid = fopen(filename, 'r');
header = fscanf(fid, '%d %d', 2);
nx = header(1);
ny = header(2);
data = fscanf(fid, '%f %f %f', [3, Inf]); % Lit [x, y, val]'
fclose(fid);

data = data'; % Transpose pour avoir colonnes x, y, val

X = reshape(data(:,1), [ny, nx]);
Y = reshape(data(:,2), [ny, nx]);
Z = reshape(data(:,3), [ny, nx]);

% Visualisation
figure('Position', [100, 100, 1000, 400]);

% On trace souvent le log de l'indicateur pour mieux voir le contraste
% Indicateur LSM = 1 / ||h||
% Echelle log = 20 * log10(Indicateur)
surf(X, Y, 20*log10(Z), 'EdgeColor', 'none');

view(2);
axis equal; axis tight;
colormap jet;
colorbar;
title('Image Linear Sampling (dB)');
xlabel('x');
ylabel('y');

% Inverser l'axe Y si besoin selon la convention du maillage
set(gca, 'YDir', 'normal');
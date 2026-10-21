% Sweep of the periodicity a and of the hole radius for the inverse design with codice/inverse_design.py.
% Runs simula_membrana (same simulation as main_apr2026.m) on a grid of (a, radius) in the measured window and saves
% sweep_a_r.mat, which the Python code reads. The thickness stays fixed.
% Partial results are saved every few geometries: if you stop it, run it again and it continues from the saved file.
clearvars; close all; clc;

a_values = 720:2:880;       % nm, nominal 800 +- 10%
r_values = 180:4:220;       % nm, nominal 200 +- 10%
dslab    = 300;             % nm, fixed
lambdav  = 1520:1:1570;     % nm, the measured window (use a finer step if the resonances are very narrow)
file_out = 'sweep_a_r.mat';

Na = numel(a_values); Nr = numel(r_values); Nl = numel(lambdav);
Txx = nan(Na,Nr,Nl); Tyy = Txx; Txy = Txx; Tyx = Txx;
if exist(file_out, 'file')                       % continue from a previous partial run if the grid is the same
    old = load(file_out);
    if isequal(old.a_values, a_values) && isequal(old.r_values, r_values) && isequal(old.lambdav, lambdav) && old.dslab == dslab
        Txx = old.Txx; Tyy = old.Tyy; Txy = old.Txy; Tyx = old.Tyx;
        fprintf('Continuing from %s\n', file_out);
    end
end

t0 = tic;
for ia = 1:Na
    for ir = 1:Nr
        if ~isnan(Txx(ia,ir,1)), continue; end   % already done
        [Txx(ia,ir,:), Tyy(ia,ir,:), Txy(ia,ir,:), Tyx(ia,ir,:)] = simula_membrana(a_values(ia), r_values(ir), dslab, lambdav);
    end
    fprintf('a = %g nm done (%d/%d), %.1f min elapsed\n', a_values(ia), ia, Na, toc(t0)/60);
    save(file_out, 'a_values', 'r_values', 'lambdav', 'dslab', 'Txx', 'Tyy', 'Txy', 'Tyx', '-v7');
end
fprintf('Saved %s\n', file_out);

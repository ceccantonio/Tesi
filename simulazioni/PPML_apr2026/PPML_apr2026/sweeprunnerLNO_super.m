clearvars
close all

flag.ploteps = 0;
flag.plotS = 0;

%% 1) Carico dati (co-pol e cross-pol) e referenza (spectra_nm: colonna 1 = nm, colonna 2 = volt)
d   = readmatrix('C:\Users\Mario\OneDrive - University of Pisa\File di Alessandro Pitanti - tesi_Mario\copolC4\run_009\spectra_nm\point_011_x110.509_y146.798.csv');
dc  = readmatrix('C:\Users\Mario\OneDrive - University of Pisa\File di Alessandro Pitanti - tesi_Mario\crosspolC4\run_010\spectra_nm\point_011_x110.509_y146.798.csv');   % lo stesso punto in cross-pol
ref = readmatrix('C:\Users\Mario\OneDrive - University of Pisa\File di Alessandro Pitanti - tesi_Mario\referenceC4\spectra_nm\copol.csv');   % la stessa referenza per co e cross
lam_d = d(:,1);    V_d = d(:,2);
lam_c = dc(:,1);   V_c = dc(:,2);
lam_r = ref(:,1);  V_r = ref(:,2);

%% 2) Divido per la referenza -> trasmissione T
ref_i = interp1(lam_r, V_r, lam_d, 'linear', 'extrap');
ref_i = max(ref_i, 0.02*max(ref_i));
T = V_d ./ ref_i;
ref_c = interp1(lam_r, V_r, lam_c, 'linear', 'extrap');
ref_c = max(ref_c, 0.02*max(ref_c));
T_c = V_c ./ ref_c;

%% 3) Media mobile
window = 1001;
T_avg  = movmean(T, window);
T_avg_c = movmean(T_c, window);
half   = floor(window/2);
valid  = (half+1):(numel(T)-half);
valid_c = (half+1):(numel(T_c)-half);
step = 200;                          % disegno un punto ogni 200
idx  = valid(1:step:end);
idx_c = valid_c(1:step:end);

%% 4) Parametri fissi, valori nominali e range da esplorare
colonna_co = 3;                      % colonna di out per la co-pol: 3 = Tyy (2 = Txx)
colonna_cross = 4;                   % colonna di out per la cross-pol: 4 = Txy (5 = Tyx)
peso_cross = 0;                      % peso della cross-pol nello scarto: 0 = lo scarto usa solo la co-pol (la cross si vede soltanto); T cross e' molto piu' piccola, servono pesi grandi (es. 100)
ds = 300;                            % spessore dello strato inciso (nm): 300 = buchi che attraversano la membrana
nsc = 2;                             % la supercella ha nsc x nsc buchi (10 = molto pesante)
lambdav = 1540:1:1564;               % lunghezze d'onda da simulare (nm)
a_nom = 1300;                        % periodo nominale (nm)
r_nom = 257;                         % raggio nominale (nm)
tol_a = 0;                        % tolleranza sul periodo: 0 = valore fisso; 0.03 = +-3%
tol_r = 0;                           % tolleranza sul raggio
step_a = 39;                         % passo sul periodo (nm)
step_r = 1;                         % passo sul raggio (nm)
sigma_nom = 60;                      % deviazione standard nominale dello shift dei buchi (nm)
tol_s = 0;                         % tolleranza su sigma: 0 = valore fisso; 0.1 = +-10%
step_s = 1;                          % passo su sigma (nm)
n_plot = 100;                         % quante simulazioni migliori disegnare a colori

a_values = round((1-tol_a)*a_nom) : step_a : round((1+tol_a)*a_nom);
r_values = round((1-tol_r)*r_nom) : step_r : round((1+tol_r)*r_nom);
sigma_values = round((1-tol_s)*sigma_nom) : step_s : round((1+tol_s)*sigma_nom);
fprintf('%d x %d x %d = %d simulazioni\n', numel(a_values), numel(r_values), numel(sigma_values), numel(a_values)*numel(r_values)*numel(sigma_values));

%% 5) Tutte le combinazioni (a, r, sigma); dopo ogni simulazione il grafico live si aggiorna
[IA, IR, IS] = ndgrid(1:numel(a_values), 1:numel(r_values), 1:numel(sigma_values));
n_sim = numel(IA);
res_v = cell(n_sim, 1);                        % l'uscita di ogni simulazione
err_v = nan(n_sim, 1);                         % lo scarto di ogni combinazione dai dati
usa_parallelo = exist('parpool') > 0;          % Parallel Computing Toolbox installato? Allora le simulazioni girano insieme sui core, altrimenti una dopo l'altra
if usa_parallelo
    pool = gcp;
    F = parallel.FevalFuture.empty;
    for k = 1:n_sim
        F(k) = parfeval(pool, @main_LNO_super, 1, a_values(IA(k)), r_values(IR(k)), ds, sigma_values(IS(k)), flag, nsc, lambdav);
    end
end

figure(4); clf                                 % grafico live: dati, ultima simulazione finita e migliore finora (sopra co-pol, sotto cross-pol)
subplot(2,1,1); hold on
plot(lam_d(idx), T_avg(idx), 'k', 'LineWidth', 2, 'DisplayName', 'dati')
h_ultima(1) = plot(nan, nan, 'b', 'LineWidth', 1.2, 'DisplayName', 'ultima simulazione');
h_migliore(1) = plot(nan, nan, 'r--', 'LineWidth', 1.5, 'DisplayName', 'migliore finora');
xlim([lambdav(1) lambdav(end)]); ylabel('T co-pol'); legend('Location', 'eastoutside')
subplot(2,1,2); hold on
plot(lam_c(idx_c), T_avg_c(idx_c), 'k', 'LineWidth', 2, 'DisplayName', 'dati')
h_ultima(2) = plot(nan, nan, 'b', 'LineWidth', 1.2, 'DisplayName', 'ultima simulazione');
h_migliore(2) = plot(nan, nan, 'r--', 'LineWidth', 1.5, 'DisplayName', 'migliore finora');
xlim([lambdav(1) lambdav(end)]); xlabel('lunghezza d''onda (nm)'); ylabel('T cross-pol'); legend('Location', 'eastoutside')
migliore = inf;
for n_fatte = 1:n_sim
    if usa_parallelo
        [k, out] = fetchNext(F);               % la prossima simulazione che finisce (non nell'ordine della griglia)
    else
        k = n_fatte;
        out = main_LNO_super(a_values(IA(k)), r_values(IR(k)), ds, sigma_values(IS(k)), flag, nsc, lambdav);
    end
    res_v{k} = out;
    T_i = interp1(lam_d(valid), T_avg(valid), out(:,1));                 % i dati alle lunghezze d'onda della simulazione (colonna 1)
    T_ic = interp1(lam_c(valid_c), T_avg_c(valid_c), out(:,1));
    err_v(k) = sum((T_i - out(:,colonna_co)).^2, 'omitnan') + peso_cross*sum((T_ic - out(:,colonna_cross)).^2, 'omitnan');   % scarto dai dati
    fprintf('(%d/%d) a = %g, r = %g, sigma = %g: scarto %.4g\n', n_fatte, n_sim, a_values(IA(k)), r_values(IR(k)), sigma_values(IS(k)), err_v(k));
    set(h_ultima(1), 'XData', out(:,1), 'YData', out(:,colonna_co));
    set(h_ultima(2), 'XData', out(:,1), 'YData', out(:,colonna_cross));
    if err_v(k) < migliore
        migliore = err_v(k); k_migliore = k;
        set(h_migliore(1), 'XData', out(:,1), 'YData', out(:,colonna_co));
        set(h_migliore(2), 'XData', out(:,1), 'YData', out(:,colonna_cross));
    end
    subplot(2,1,1)
    title(sprintf('%d/%d: a = %g, r = %g, sigma = %g (scarto %.4g)   |   migliore finora: a = %g, r = %g, sigma = %g (scarto %.4g)', n_fatte, n_sim, ...
        a_values(IA(k)), r_values(IR(k)), sigma_values(IS(k)), err_v(k), a_values(IA(k_migliore)), r_values(IR(k_migliore)), sigma_values(IS(k_migliore)), migliore), 'FontSize', 9)
    drawnow
end
res = reshape(res_v, size(IA));
err = reshape(err_v, size(IA));

[~, ordine] = sort(err(:));                    % combinazioni dalla migliore alla peggiore
[ia_b, ir_b, is_b] = ind2sub(size(err), ordine(1));
fprintf('\nMIGLIORE: a = %g, r = %g, sigma = %g (scarto %.4g)\n', a_values(ia_b), r_values(ir_b), sigma_values(is_b), err(ia_b, ir_b, is_b));

%% 6) Grafici: sopra la co-pol, sotto la cross-pol
figure(1); clf
n_curve = min(n_plot, numel(ordine));
colori = turbo(n_curve);                       % un colore diverso per ogni curva (dalla migliore alla peggiore)
stili = {'-', '--', ':'};                      % e lo stile cambia ogni curva, per distinguere quelle con colori vicini
colonne = [colonna_co colonna_cross];
for p = 1:2
    subplot(2,1,p); hold on
    for j = 1:n_curve                          % le migliori a colori
        [ia, ir, is] = ind2sub(size(err), ordine(j));
        plot(res{ia, ir, is}(:,1), res{ia, ir, is}(:,colonne(p)), stili{mod(j-1, 3) + 1}, 'Color', colori(j,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('a = %g, r = %g, sigma = %g', a_values(ia), r_values(ir), sigma_values(is)))
    end
end
subplot(2,1,1)
plot(lam_d(idx), T_avg(idx), 'k', 'LineWidth', 2, 'DisplayName', 'dati')       % i dati sopra tutto
xlim([lambdav(1) lambdav(end)]); ylabel('T co-pol'); legend('Location', 'eastoutside')   % asse x = intervallo simulato
subplot(2,1,2)
plot(lam_c(idx_c), T_avg_c(idx_c), 'k', 'LineWidth', 2, 'DisplayName', 'dati')
xlim([lambdav(1) lambdav(end)]); xlabel('lunghezza d''onda (nm)'); ylabel('T cross-pol'); legend('Location', 'eastoutside')

if numel(a_values) > 1 && numel(r_values) > 1  % la mappa serve solo con piu' valori di a e r
    figure(2)
    imagesc(r_values, a_values, err(:,:,is_b)); colorbar
    xlabel('r (nm)'); ylabel('a (nm)'); title(sprintf('scarto tra dati e simulazione (sigma = %g nm)', sigma_values(is_b)))
end

if numel(sigma_values) > 1                     % scarto migliore (sugli altri due parametri) per ogni sigma
    figure(3)
    err_s = squeeze(min(min(err, [], 1), [], 2));
    plot(sigma_values, err_s, 'o-', 'LineWidth', 1.2)
    xlabel('sigma (nm)'); ylabel('scarto minimo su a e r'); title('scarto in funzione dello shift dei buchi')
end

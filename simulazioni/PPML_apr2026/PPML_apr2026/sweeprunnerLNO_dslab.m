clearvars
close all

flag.ploteps = 0;
flag.plotS = 0;

%% 1) Carico dato e referenza (spectra_nm: colonna 1 = nm, colonna 2 = volt)
d   = readmatrix('C:\Users\Mario\OneDrive - University of Pisa\File di Alessandro Pitanti - tesi_Mario\copolC4\run_009\spectra_nm\point_027_x109.010_y146.493.csv');
ref = readmatrix('C:\Users\Mario\OneDrive - University of Pisa\File di Alessandro Pitanti - tesi_Mario\referenceC4\spectra_nm\copol.csv');   % sempre la co-pol
lam_d = d(:,1);    V_d = d(:,2);
lam_r = ref(:,1);  V_r = ref(:,2);

%% 2) Divido per la referenza -> trasmissione T
ref_i = interp1(lam_r, V_r, lam_d, 'linear', 'extrap');
ref_i = max(ref_i, 0.02*max(ref_i));
T = V_d ./ ref_i;

%% 3) Media mobile
window = 10001;
T_avg  = movmean(T, window);
half   = floor(window/2);
valid  = (half+1):(numel(T)-half);
step = 200;                          % disegno un punto ogni 200
idx  = valid(1:step:end);

%% 4) Valori nominali, range da esplorare e componente da confrontare
colonna = 2;                         % colonna di out: 2 = Txx, 3 = Tyy (co-pol), 4 = Txy, 5 = Tyx (cross-pol)
a_nom = 1240;                        % periodo nominale (nm)
r_nom = 229;                         % raggio nominale (nm)
d_nom = 30;                          % spessore dello strato inciso nominale (nm); lo spessore totale resta 300 nm
tol_a = 0.0;                        % tolleranza sul periodo: 0 = valore fisso; 0.05 = +-5%
tol_r = 0;                        % tolleranza sul raggio
tol_d = 0;                        % tolleranza su dslab
step_a = 50;                         % passo sul periodo (nm)
step_r = 10;                         % passo sul raggio (nm)
step_d = 10;                          % passo sullo spessore dslab (nm)
n_plot = 100;                         % quante simulazioni migliori disegnare a colori

a_values = round((1-tol_a)*a_nom) : step_a : round((1+tol_a)*a_nom);
r_values = round((1-tol_r)*r_nom) : step_r : round((1+tol_r)*r_nom);
d_values = round((1-tol_d)*d_nom) : step_d : round((1+tol_d)*d_nom);
fprintf('%d x %d x %d = %d simulazioni\n', numel(a_values), numel(r_values), numel(d_values), numel(a_values)*numel(r_values)*numel(d_values));

%% 5) Tutte le combinazioni (a, r, dslab); dopo ogni simulazione il grafico live si aggiorna
[IA, IR, ID] = ndgrid(1:numel(a_values), 1:numel(r_values), 1:numel(d_values));
n_sim = numel(IA);
res_v = cell(n_sim, 1);                        % l'uscita di ogni simulazione
err_v = nan(n_sim, 1);                         % lo scarto di ogni combinazione dai dati
usa_parallelo = exist('parpool') > 0;          % Parallel Computing Toolbox installato? Allora le simulazioni girano insieme sui core, altrimenti una dopo l'altra
if usa_parallelo
    pool = gcp;
    F = parallel.FevalFuture.empty;
    for k = 1:n_sim
        F(k) = parfeval(pool, @main_LNO_func, 1, a_values(IA(k)), r_values(IR(k)), d_values(ID(k)), flag);
    end
end

figure(4); clf; hold on                        % grafico live: dati, ultima simulazione finita e migliore finora
plot(lam_d(idx), T_avg(idx), 'k', 'LineWidth', 2, 'DisplayName', 'dati')
h_ultima = plot(nan, nan, 'b', 'LineWidth', 1.2, 'DisplayName', 'ultima simulazione');
h_migliore = plot(nan, nan, 'r--', 'LineWidth', 1.5, 'DisplayName', 'migliore finora');
xlabel('lunghezza d''onda (nm)'); ylabel('trasmissione T'); legend('Location', 'eastoutside')
migliore = inf;
for n_fatte = 1:n_sim
    if usa_parallelo
        [k, out] = fetchNext(F);               % la prossima simulazione che finisce (non nell'ordine della griglia)
    else
        k = n_fatte;
        out = main_LNO_func(a_values(IA(k)), r_values(IR(k)), d_values(ID(k)), flag);
    end
    res_v{k} = out;
    if n_fatte == 1, xlim([out(1,1) out(end,1)]); end               % asse x = intervallo di lunghezze d'onda simulato
    T_i = interp1(lam_d(valid), T_avg(valid), out(:,1));                 % i dati alle lunghezze d'onda della simulazione (colonna 1)
    err_v(k) = sum((T_i - out(:,colonna)).^2, 'omitnan');                % scarto dalla componente scelta
    fprintf('(%d/%d) a = %g, r = %g, dslab = %g: scarto %.4g\n', n_fatte, n_sim, a_values(IA(k)), r_values(IR(k)), d_values(ID(k)), err_v(k));
    set(h_ultima, 'XData', out(:,1), 'YData', out(:,colonna));
    if err_v(k) < migliore
        migliore = err_v(k); k_migliore = k;
        set(h_migliore, 'XData', out(:,1), 'YData', out(:,colonna));
    end
    title(sprintf('%d/%d: a = %g, r = %g, dslab = %g (scarto %.4g)   |   migliore finora: a = %g, r = %g, dslab = %g (scarto %.4g)', n_fatte, n_sim, ...
        a_values(IA(k)), r_values(IR(k)), d_values(ID(k)), err_v(k), a_values(IA(k_migliore)), r_values(IR(k_migliore)), d_values(ID(k_migliore)), migliore), 'FontSize', 9)
    drawnow
end
res = reshape(res_v, size(IA));
err = reshape(err_v, size(IA));

[~, ordine] = sort(err(:));                    % combinazioni dalla migliore alla peggiore
[ia_b, ir_b, id_b] = ind2sub(size(err), ordine(1));
fprintf('\nMIGLIORE: a = %g, r = %g, dslab = %g (scarto %.4g)\n', a_values(ia_b), r_values(ir_b), d_values(id_b), err(ia_b, ir_b, id_b));

%% 6) Grafici
figure(1); clf; hold on
n_curve = min(n_plot, numel(ordine));
colori = turbo(n_curve);                       % un colore diverso per ogni curva (dalla migliore alla peggiore)
stili = {'-', '--', ':'};                      % e lo stile cambia ogni curva, per distinguere quelle con colori vicini
for j = 1:n_curve                              % le migliori a colori
    [ia, ir, id] = ind2sub(size(err), ordine(j));
    plot(res{ia, ir, id}(:,1), res{ia, ir, id}(:,colonna), stili{mod(j-1, 3) + 1}, 'Color', colori(j,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('a = %g, r = %g, dslab = %g', a_values(ia), r_values(ir), d_values(id)))
end
plot(lam_d(idx), T_avg(idx), 'k', 'LineWidth', 2, 'DisplayName', 'dati')     % i dati sopra tutto
xlim([res{1}(1,1) res{1}(end,1)]); xlabel('lunghezza d''onda (nm)'); ylabel('trasmissione T'); legend('Location', 'eastoutside')   % asse x = intervallo simulato

if numel(a_values) > 1 && numel(r_values) > 1  % la mappa serve solo con piu' valori di a e r
    figure(2)
    imagesc(r_values, a_values, err(:,:,id_b)); colorbar
    xlabel('r (nm)'); ylabel('a (nm)'); title(sprintf('scarto tra dati e simulazione (dslab = %g nm)', d_values(id_b)))
end

if numel(d_values) > 1                         % scarto migliore (sugli altri due parametri) per ogni dslab
    figure(3)
    err_d = squeeze(min(min(err, [], 1), [], 2));
    plot(d_values, err_d, 'o-', 'LineWidth', 1.2)
    xlabel('dslab (nm)'); ylabel('scarto minimo su a e r'); title('scarto in funzione dello spessore dslab')
end

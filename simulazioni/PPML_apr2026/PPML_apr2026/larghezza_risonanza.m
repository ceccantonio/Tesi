clearvars
close all

flag.ploteps = 0;
flag.plotS = 0;

%% 1) Dati: file del campione e della referenza (spectra_nm: colonna 1 = nm, colonna 2 = volt)
file_dati = 'C:\Users\Mario\OneDrive - University of Pisa\File di Alessandro Pitanti - tesi_Mario\copolC4\run_009\spectra_nm\point_027_x109.010_y146.493.csv';
file_ref  = 'C:\Users\Mario\OneDrive - University of Pisa\File di Alessandro Pitanti - tesi_Mario\referenceC4\spectra_nm\copol.csv';   % sempre la co-pol

%% 2) Cosa simulare e cosa confrontare
colonna = 3;                         % colonna di out: 2 = Txx, 3 = Tyy (co-pol), 4 = Txy, 5 = Tyx (cross-pol)
a_values = 1225;                     % periodi da simulare (nm): un valore o un elenco, es. 1200:5:1250
r_values = 261;                      % raggi (nm)
d_values = 20:10:150;                % spessori dslab (nm)
finestra = 1001;                     % media mobile dei dati (campioni; 0.8 nm circa): piu' stretta di quella dello sweep, cosi' non allarga il picco e non perde i bordi
prominenza = 0.03;                   % un dip conta se e' profondo almeno cosi' (in T)
peso_larghezza = 1;                  % importanza della larghezza nel compromesso (0 = ignorala)
peso_profondita = 1;                 % importanza della profondita' (0 = ignorala)

%% 3) Dati -> T e larghezza della risonanza
d   = readmatrix(file_dati);
ref = readmatrix(file_ref);
lam_d = d(:,1);    V_d = d(:,2);
ref_i = interp1(ref(:,1), ref(:,2), lam_d, 'linear', 'extrap');
ref_i = max(ref_i, 0.02*max(ref_i));
T = V_d ./ ref_i;
T_avg = movmean(T, finestra);
half  = floor(finestra/2);
valid = (half+1):(numel(T)-half);
lam_v = lam_d(valid);   T_v = T_avg(valid);
[w_dati, lam0_dati, prof_dati] = larghezza_dip(lam_v, T_v, prominenza);
fprintf('DATI: risonanza a %.2f nm, larghezza (FWHM) %.2f nm, profondita'' %.3f\n', lam0_dati, w_dati, prof_dati);

%% 4) Simulazioni: larghezza e profondita' del dip di ciascuna
[A, R, D] = ndgrid(a_values, r_values, d_values);
P = [A(:) R(:) D(:)];                          % una riga per simulazione: a, r, dslab
n_sim = size(P, 1);
res = cell(n_sim, 1);
w_sim = nan(n_sim, 1);   lam0_sim = nan(n_sim, 1);   prof_sim = nan(n_sim, 1);
for k = 1:n_sim
    out = main_LNO_func(P(k,1), P(k,2), P(k,3), flag);
    res{k} = out;
    [~, l0] = larghezza_dip(out(:,1), out(:,colonna), prominenza);          % dove cade il dip simulato
    if ~isnan(l0)                                                          % larghezza e profondita' nella stessa finestra dei dati, centrata sul dip: confrontabili
        m = out(:,1) >= l0 - (lam0_dati - lam_v(1)) & out(:,1) <= l0 + (lam_v(end) - lam0_dati);
        [w_sim(k), lam0_sim(k), prof_sim(k)] = larghezza_dip(out(m,1), out(m,colonna), prominenza);
    end
    fprintf('(%d/%d) a = %g, r = %g, dslab = %g: risonanza a %.2f nm, larghezza %.2f nm, profondita'' %.3f\n', k, n_sim, P(k,:), lam0_sim(k), w_sim(k), prof_sim(k));
end
passo = median(diff(res{1}(:,1)));             % passo in lunghezza d'onda della simulazione (nm)
if any(w_sim < 2*passo)
    fprintf('ATTENZIONE: dove la larghezza e'' sotto %.1f nm (due passi della simulazione) il dip non e'' risolto: il valore non e'' affidabile\n', 2*passo);
end

%% 5) Compromesso tra larghezza e profondita': scarto relativo dai dati di entrambe, dal migliore al peggiore
costo = peso_larghezza*((w_sim - w_dati)/w_dati).^2 + peso_profondita*((prof_sim - prof_dati)/prof_dati).^2;
[~, ordine] = sort(costo);                     % le simulazioni senza dip (NaN) vanno in fondo
tabella = table(P(ordine,1), P(ordine,2), P(ordine,3), lam0_sim(ordine), w_sim(ordine), w_sim(ordine)/w_dati, prof_sim(ordine), prof_sim(ordine)/prof_dati, costo(ordine), ...
    'VariableNames', {'a', 'r', 'dslab', 'lambda0_nm', 'FWHM_nm', 'FWHM_su_dati', 'profondita', 'prof_su_dati', 'costo'});
disp(tabella)
k = ordine(1);
fprintf('MIGLIOR COMPROMESSO: a = %g, r = %g, dslab = %g (FWHM %.2f nm contro %.2f nm, profondita'' %.3f contro %.3f)\n', P(k,:), w_sim(k), w_dati, prof_sim(k), prof_dati);

%% 6) Grafici
figure(1); clf; hold on                        % il miglior compromesso, con i dati e le larghezze segnate
plot(lam_v, T_v, 'k', 'LineWidth', 2, 'DisplayName', sprintf('dati (FWHM %.1f nm, prof. %.2f)', w_dati, prof_dati))
plot(res{k}(:,1), res{k}(:,colonna), 'r', 'LineWidth', 1.5, 'DisplayName', sprintf('a = %g, r = %g, dslab = %g (FWHM %.1f nm, prof. %.2f)', P(k,:), w_sim(k), prof_sim(k)))
xline(lam0_dati - w_dati/2, 'k:', 'HandleVisibility', 'off');   xline(lam0_dati + w_dati/2, 'k:', 'HandleVisibility', 'off')
xline(lam0_sim(k) - w_sim(k)/2, 'r:', 'HandleVisibility', 'off'); xline(lam0_sim(k) + w_sim(k)/2, 'r:', 'HandleVisibility', 'off')
xlim([res{1}(1,1) res{1}(end,1)]); xlabel('lunghezza d''onda (nm)'); ylabel('trasmissione T'); legend('Location', 'eastoutside')

variato = find([numel(a_values) numel(r_values) numel(d_values)] > 1);     % larghezza e profondita' contro il parametro che varia (se ne varia uno solo)
nomi = {'a (nm)', 'r (nm)', 'dslab (nm)'};
if isscalar(variato), x = P(:, variato); etichetta = nomi{variato}; else, x = (1:n_sim)'; etichetta = 'simulazione'; end
figure(2); clf
subplot(2,1,1); hold on
plot(x, w_sim, 'o-', 'LineWidth', 1.5, 'DisplayName', 'simulazione'); yline(w_dati, 'k', 'LineWidth', 2, 'DisplayName', 'dati')
plot(x(k), w_sim(k), 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'DisplayName', 'miglior compromesso')
ylabel('larghezza FWHM (nm)'); legend('Location', 'best'); grid on
subplot(2,1,2); hold on
plot(x, prof_sim, 'o-', 'LineWidth', 1.5, 'DisplayName', 'simulazione'); yline(prof_dati, 'k', 'LineWidth', 2, 'DisplayName', 'dati')
plot(x(k), prof_sim(k), 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'DisplayName', 'miglior compromesso')
xlabel(etichetta); ylabel('profondita'' del dip (T)'); grid on

figure(3); clf; hold on                        % le due grandezze insieme: piu' vicino alla stella dei dati = meglio
scatter(w_sim, prof_sim, 60, x, 'filled', 'DisplayName', 'simulazioni'); colorbar; ylabel(colorbar, etichetta)
plot(w_dati, prof_dati, 'kp', 'MarkerSize', 18, 'MarkerFaceColor', 'k', 'DisplayName', 'dati')
plot(w_sim(k), prof_sim(k), 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'DisplayName', 'miglior compromesso')
xlabel('larghezza FWHM (nm)'); ylabel('profondita'' del dip (T)'); legend('Location', 'best'); grid on

function [w, lam0, prof] = larghezza_dip(lam, T, prominenza)   % FWHM (nm) del dip piu' profondo di T(lam), a meta' della sua prominenza
[picchi, posizioni, larghezze, prom] = findpeaks(-T(:), lam(:), 'MinPeakProminence', prominenza, 'WidthReference', 'halfprom');
if isempty(picchi)
    w = NaN; lam0 = NaN; prof = NaN;
    return
end
[prof, i] = max(prom);
w = larghezze(i);
lam0 = posizioni(i);
end

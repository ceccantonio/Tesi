clearvars
close all

flag.ploteps = 0;
flag.plotS = 0;

%% 1) Simulazione con supercella e difetti
a = 1180;                            % periodo (nm)
r = 280;                             % raggio dei buchi (nm)
ds = 300;                             % spessore dello strato inciso (nm)
sigma = 0;                          % deviazione standard dello shift dei buchi dalla posizione di equilibrio (nm); 0 = reticolo perfetto
nsc = 2 
% la supercella ha nsc x nsc buchi (10 = molto pesante)
lambdav = 1540:4:1566;               % lunghezze d'onda da simulare (nm)
colonna_co = 3;                      % colonna di out per la co-pol: 3 = Tyy (2 = Txx)
colonna_cross = 4;                   % colonna di out per la cross-pol: 4 = Txy (5 = Tyx)

%% 2) Misure e referenza (spectra_nm: colonna 1 = nm, colonna 2 = volt)
file_co    = "C:\Users\Mario\OneDrive - University of Pisa\File di Alessandro Pitanti - tesi_Mario\copolC4\run_009\spectra_nm\point_011_x110.509_y146.798.csv";
file_cross = "C:\Users\Mario\OneDrive - University of Pisa\File di Alessandro Pitanti - tesi_Mario\crosspolC4\run_010\spectra_nm\point_011_x110.509_y146.798.csv";   % lo stesso punto misurato in cross-pol
file_ref   = 'C:\Users\Mario\OneDrive - University of Pisa\File di Alessandro Pitanti - tesi_Mario\referenceC4\spectra_nm\copol.csv';   % la stessa referenza per co e cross
finestra = 10001;                    % media mobile delle misure (campioni)

[lam_co, T_co] = trasmissione(file_co, file_ref, finestra);
[lam_cr, T_cr] = trasmissione(file_cross, file_ref, finestra);

%% 3) Simulazione (una sola: contiene sia la co-pol sia la cross-pol)
out = main_LNO_super(a, r, ds, sigma, flag, nsc, lambdav);

%% 4) Confronto: sopra la co-pol, sotto la cross-pol
etichetta = sprintf('simulazione: a = %g, r = %g, dslab = %g, sigma = %g nm, %dx%d buchi', a, r, ds, sigma, nsc, nsc);
figure(1); clf
subplot(2,1,1); hold on
plot(lam_co, T_co, 'k', 'LineWidth', 2, 'DisplayName', 'misura co-pol')
plot(out(:,1), out(:,colonna_co), 'r', 'LineWidth', 1.5, 'DisplayName', etichetta)
xlim([lambdav(1) lambdav(end)]); ylabel('T co-pol'); legend('Location', 'eastoutside')
subplot(2,1,2); hold on
plot(lam_cr, T_cr, 'k', 'LineWidth', 2, 'DisplayName', 'misura cross-pol')
plot(out(:,1), out(:,colonna_cross), 'r', 'LineWidth', 1.5, 'DisplayName', etichetta)
xlim([lambdav(1) lambdav(end)]); xlabel('lunghezza d''onda (nm)'); ylabel('T cross-pol'); legend('Location', 'eastoutside')

function [lam, T] = trasmissione(file_dati, file_ref, finestra)   % T = misura / referenza, con media mobile, senza i bordi (circa 500 punti)
d   = readmatrix(file_dati);
ref = readmatrix(file_ref);
ref_i = interp1(ref(:,1), ref(:,2), d(:,1), 'linear', 'extrap');
ref_i = max(ref_i, 0.02*max(ref_i));
T_avg = movmean(d(:,2) ./ ref_i, finestra);
half  = floor(finestra/2);
valid = (half+1):(size(d,1)-half);                       % senza i bordi, dove la media non e' completa
valid = valid(1:max(1, round(numel(valid)/500)):end);
lam = d(valid,1);
T = T_avg(valid);
end

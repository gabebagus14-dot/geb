%% ========================================================================
%                          PETROPHYSICAL ANALYSIS 
%                                TUGAS AKHIR 
%  ========================================================================
clc; clear; close all;

%% ========================================================================
%               1. READ LAS FILE (DYNAMIC & UNIVERSAL READER)
% =========================================================================
% A. Pilih File secara Interaktif (Pop-up GUI)
[file, path] = uigetfile('*.las', 'Pilih File Data LAS');
if isequal(file,0), error('ANALISIS BERHENTI: Tidak ada file yang dipilih.'); end
filename = fullfile(path, file);
fprintf('\n=> Membaca file: %s\n', file);

% B. Membaca Header untuk Mencari Urutan Kolom (Curve Information ~C)
fid = fopen(filename,'r');
curve_names = {};
in_curve = false;

while ~feof(fid)
    line = strtrim(fgetl(fid));
    if startsWith(upper(line), '~C')
        in_curve = true; continue;
    elseif startsWith(upper(line), '~') && in_curve
        in_curve = false;
    end
    if startsWith(upper(line), '~A')
        break; % Berhenti membaca header, masuk ke data angka
    end
    
    % Mengekstrak singkatan nama log (Mnemonic)
    if in_curve && ~isempty(line) && ~startsWith(line, '#')
        tokens = regexp(line, '^([^.\s]+)', 'tokens');
        if ~isempty(tokens)
            curve_names{end+1} = upper(tokens{1}{1});
        end
    end
end

% C. Membaca Data Angka (ASCII ~A)
ascii_data = [];
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    nums = str2num(line);
    if ~isempty(nums), ascii_data = [ascii_data; nums]; end
end
fclose(fid);

% D. Membersihkan Data NULL (-999.25 atau sejenisnya)
NULL = -999.25; ascii_data(ascii_data <= NULL) = NaN;

% E. Identifikasi Kolom Otomatis (Dynamic Column Mapping)
get_col = @(aliases) find(ismember(curve_names, aliases), 1);

% Mencari kolom Depth
idx_dept = get_col({'DEPT', 'DEPTH', 'DEP'});
if isempty(idx_dept), error('ANALISIS BERHENTI: Kolom Depth tidak ditemukan!'); end
DEPT = ascii_data(:, idx_dept);

% Mendefinisikan kolom dengan alias umum dari berbagai perusahaan (Schlumberger, Halliburton, dll)
idx_gr   = get_col({'GR', 'GAMMA', 'CGR', 'SGR'});
idx_lld  = get_col({'LLD', 'RT', 'ILD', 'RD', 'AHT90', 'RDEP'});
idx_lls  = get_col({'LLS', 'RXO', 'ILM', 'RS', 'AHT10', 'RMIC'});
idx_rhob = get_col({'RHOB', 'ZDEN', 'DEN', 'RHOZ', 'DENB'});
idx_nphi = get_col({'NPHI', 'CNC', 'PHIN', 'NPOR', 'HNPO'});
idx_dt   = get_col({'DT', 'AC', 'DTC', 'DTCO', 'ACCO'});

% F. Ekstraksi Data (Otomatis mengisi NaN jika log tidak tersedia di file tersebut)
n_rows = length(DEPT);
if ~isempty(idx_gr),   GR = ascii_data(:, idx_gr);     else, GR = NaN(n_rows, 1);   fprintf('   [!] Peringatan: Log GR tidak ditemukan.\n'); end
if ~isempty(idx_lld),  LLD = ascii_data(:, idx_lld);   else, LLD = NaN(n_rows, 1);  fprintf('   [!] Peringatan: Log Resistivitas Dalam (RT) tidak ditemukan.\n'); end
if ~isempty(idx_lls),  LLS = ascii_data(:, idx_lls);   else, LLS = NaN(n_rows, 1);  fprintf('   [!] Peringatan: Log Resistivitas Dangkal (RXO) tidak ditemukan.\n'); end
if ~isempty(idx_rhob), RHOB = ascii_data(:, idx_rhob); else, RHOB = NaN(n_rows, 1); fprintf('   [!] Peringatan: Log Densitas (RHOB) tidak ditemukan.\n'); end
if ~isempty(idx_nphi), NPHI = ascii_data(:, idx_nphi); else, NPHI = NaN(n_rows, 1); fprintf('   [!] Peringatan: Log Neutron (NPHI) tidak ditemukan.\n'); end
if ~isempty(idx_dt),   DT = ascii_data(:, idx_dt);     else, DT = NaN(n_rows, 1);   fprintf('   [!] Peringatan: Log Sonic (DT) tidak ditemukan.\n'); end

% G. Kalkulasi Step Kedalaman Dinamis
dz = median(diff(DEPT), 'omitnan'); if isnan(dz) || dz == 0, dz = 0.5; end
fprintf('=> Selesai! Data berhasil diekstrak. Resolusi kedalaman (dz) = %g ft\n', dz);

% NILAI CUT-OFF MANUAL (SESUAI GAMBAR JURNAL)
CO_Phi = 0.065; 
CO_Vsh = 0.37;  
CO_Sw  = 0.55;

%% ========================================================================
%                        FORMATION BOUNDARY  
% =========================================================================
Top_Formasi  = 5500; 
Base_Formasi = 7400; 
flag_reservoir = (DEPT >= Top_Formasi) & (DEPT <= Base_Formasi);
zone = flag_reservoir & ~isnan(GR) & ~isnan(LLD) & ~isnan(RHOB) & ~isnan(NPHI);

%% ======================================================================== 
%                  PERHITUNGAN PETROFISIK (VSH & PHI) - EVALUASI MURNI
% =========================================================================

% -------------------------------------------------------------------------
% BAGIAN 1: PENCARIAN BASELINE & VSHALE (TETAP DIPERTAHANKAN)
% -------------------------------------------------------------------------
GR_valid_all = GR(~isnan(GR));
if ~isempty(GR_valid_all)
    GR_shale_cutoff = prctile(GR_valid_all, 95);
    idx_true_shale = (GR >= GR_shale_cutoff) & ~isnan(RHOB) & ~isnan(NPHI) & ~isnan(LLD);
    if sum(idx_true_shale) < 5
        NPHI_valid = NPHI(~isnan(NPHI)); NPHI_shale_cutoff = prctile(NPHI_valid, 95);
        idx_true_shale = (NPHI >= NPHI_shale_cutoff) & ~isnan(RHOB) & ~isnan(LLD);
    end
end

rho_sh  = median(RHOB(idx_true_shale), 'omitnan');
phi_Nsh = median(NPHI(idx_true_shale), 'omitnan');
Rsh     = median(LLD(idx_true_shale), 'omitnan');

if rho_sh > 2.71 || isnan(rho_sh), rho_sh = prctile(RHOB(~isnan(RHOB)), 90); end
if isnan(phi_Nsh), phi_Nsh = prctile(NPHI(~isnan(NPHI)), 90); end 
if isnan(Rsh), Rsh = median(LLD(~isnan(LLD)), 'omitnan'); end

fprintf('=> Rho_sh = %.4f g/cc, NPHI_sh = %.4f v/v, Rsh = %.4f ohm.m\n', rho_sh, phi_Nsh, Rsh);

% Perhitungan Vshale (Linier MURNI TANPA BATAS)
GR_valid = GR(~isnan(GR) & GR > 0);
GR_clean = prctile(GR_valid, 5); 
GR_shale = prctile(GR_valid, 95); 
fprintf('=> GR Clean: %.2f API | GR Shale: %.2f API\n', GR_clean, GR_shale);

Vsh = (GR - GR_clean) ./ (GR_shale - GR_clean); 
% Keterangan: Batas bawah 0 dan batas atas 1 telah DIHAPUS.

% -------------------------------------------------------------------------
% BAGIAN 2: PERHITUNGAN POROSITAS EFEKTIF (SESUAI LAPORAN & EXCEL)
% -------------------------------------------------------------------------
rho_ma = 2.71; % Matriks Karbonat (Limestone)
rho_f = 1.0;   % Fluida Fresh Water

% Menghitung Porositas Log Murni (Sebelum Koreksi)
phi_D = (rho_ma - RHOB) ./ (rho_ma - rho_f); 

% Menghitung Baseline Porositas Density Shale (phi_Dsh)
phi_Dsh = (rho_ma - rho_sh) ./ (rho_ma - rho_f);

% --- KODE UNTUK MENCETAK ACUAN EXCEL KE COMMAND WINDOW ---
fprintf('\n========================================================\n');
fprintf('   ACUAN HITUNGAN MANUAL EXCEL (POROSITAS & SHALE)\n');
fprintf('========================================================\n');
fprintf(' Gunakan nilai ini sebagai konstanta di Excel Anda:\n');
fprintf(' - Densitas Shale (pb_sh / rho_sh)   : %.4f g/cc\n', rho_sh);
fprintf(' - Porositas Neutron Shale (Por_Nsh) : %.4f v/v\n', phi_Nsh);
fprintf(' - Porositas Density Shale (phi_Dsh) : %.4f v/v\n', phi_Dsh);
fprintf('========================================================\n\n');

% Mengaplikasikan Koreksi Volume Shale (Sesuai rumus laporan Anda)
phi_Dc = phi_D - (Vsh .* phi_Dsh); 
phi_Nc = NPHI - (Vsh .* phi_Nsh); 

% Menghitung Porositas Efektif Total (Arithmetic Mean)
phi_eff = (phi_Nc + phi_Dc) ./ 2;

%% ========================================================================
%               WATER SATURATION & PICKETT PLOT - EVALUASI MURNI
% =========================================================================
% Parameter Archie
a = 1.0; m = 2.0; n = 2.0;   

% Penentuan Rw (MURNI TANPA FILTER POROSITAS)
z_top_water = 7200; z_bot_water = 7220; 
idx_water = (DEPT >= z_top_water) & (DEPT <= z_bot_water) & ~isnan(LLD) & ~isnan(phi_eff);

% Pastikan rumus Rwa di Excel Anda sama persis dengan ini:
Rwa_water_zone = (LLD(idx_water) .* (phi_eff(idx_water).^m)) ./ a;
Rw = median(Rwa_water_zone, 'omitnan');

if isnan(Rw)
    error('ANALISIS BERHENTI: Gagal menghitung Rw. Tidak ada data log yang valid di zona air.');
end
fprintf('=> Nilai Rw terhitung murni dari zona air: %.6f ohm.m\n', Rw);

% Perhitungan Saturasi Air (MURNI TANPA BATAS)
Sw = ( (a .* Rw) ./ ((phi_eff.^m) .* LLD) ) .^ (1/n);
% Keterangan: Batas Sw < 0.05 dan Sw > 1.0 telah DIHAPUS.

% --- KODE PEMBUATAN FIGURE PICKETT PLOT (REVISI JURNAL) ---
fig_pickett = figure('Name', 'Pickett Plot - Final Validation', 'Color', 'w', 'Position', [200 200 900 800]);
h_all = loglog(LLD(zone), phi_eff(zone), 'o', 'MarkerSize', 4, 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', 'none'); 
hold on; grid on;

phi_line_range = logspace(log10(0.01), log10(1.0), 100); 
sw_lines = [1.0, 0.5, 0.25]; 
colors = {'b', 'g', 'r'};   

for i = 1:length(sw_lines)
    sw_val = sw_lines(i);
    rt_calc = (a * Rw) ./ ((phi_line_range.^m) .* (sw_val^n));
    plot(rt_calc, phi_line_range, 'Color', colors{i}, 'LineWidth', 2.5); % Garis dipertebal
    
    text_phi = 0.05; 
    text_rt = (a * Rw) / (text_phi^m * sw_val^n);
    if text_rt < 2000 
        % Font persentase Sw dinaikkan dari 9 ke 13
        text(text_rt, text_phi, [' Sw ' num2str(sw_val*100) '%'], ...
             'Color', colors{i}, 'FontWeight', 'bold', 'FontSize', 13, 'BackgroundColor', 'w');
    end
end

if sum(idx_water) > 0
    h_water = loglog(LLD(idx_water), phi_eff(idx_water), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k'); 
end
plot(Rw, 1, 'rp', 'MarkerSize', 16, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'r'); 

% Rw tidak dibulatkan (menggunakan 6 desimal) dan font judul diperbesar menjadi 16
title(['Pickett Plot | Rw = ' num2str(Rw, '%.6f') ' ohm.m'], 'FontSize', 16, 'FontWeight', 'bold');

% Label X dan Y diperbesar menjadi 15
xlabel('Deep Resistivity, R_t (ohm.m)', 'FontWeight', 'bold', 'FontSize', 15);
ylabel('Effective Porosity, \phi_e (v/v)', 'FontWeight', 'bold', 'FontSize', 15);

% Angka skala (tick marks) diperbesar menjadi 13
set(gca, 'XLim', [0.1 2000], 'YLim', [0.01 1.0], 'FontSize', 13);
grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);

if sum(idx_water) > 0
    % Font legend diperbesar menjadi 12
    legend([h_all, h_water], {'All Data', 'Water Zone (Verified)'}, 'Location', 'southwest', 'FontSize', 12);
else
    legend(h_all, {'All Data'}, 'Location', 'southwest', 'FontSize', 12);
end
hold off;

% =========================================================================
% EXPORT PICKETT PLOT KE RESOLUSI TINGGI
% =========================================================================
exportgraphics(fig_pickett, 'Pickett_Plot_600DPI.png', 'Resolution', 600);
fprintf('=> Sukses! Gambar Pickett Plot (High-Res 600 DPI) telah disimpan.\n');
%% ========================================================================
%       VISUALISASI CUT-OFF (RESOLUSI TINGGI UNTUK JURNAL)
% =========================================================================
% Definisikan variabel data formasi terlebih dahulu
vsh_form = Vsh(zone); 
phi_form = phi_eff(zone); 
sw_form  = Sw(zone);

% --- PENGATURAN FONT GLOBAL UNTUK CUT-OFF ---
Title_FS = 20;   % Ukuran Font Judul
Label_FS = 16;   % Ukuran Font Sumbu (X dan Y Label)
Tick_FS  = 14;   % Ukuran Font Angka Skala
Text_FS  = 15;   % Ukuran Font Teks Keterangan

% Perbesar ukuran figure agar 3 subplot tidak berdesakan
fig_cutoff = figure('Name', 'Analisis Cut-Off Statistik', 'Color', 'w', 'Position', [50 100 1600 550], 'NumberTitle', 'off');

% ---------------------------------------------------------
% Subplot 1: Histogram Porositas
% ---------------------------------------------------------
subplot(1,3,1);
histogram(phi_form(vsh_form <= CO_Vsh), 30, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k', 'LineWidth', 1.2); hold on; grid on;
xline(CO_Phi, 'r-', 'LineWidth', 3.5);
title('Histogram PHIE', 'FontWeight', 'bold', 'FontSize', Title_FS);
xlabel('Effective Porosity (PHIE)', 'FontWeight', 'bold', 'FontSize', Label_FS); 
ylabel('Frekuensi', 'FontWeight', 'bold', 'FontSize', Label_FS);

% KUNCI: Kunci sumbu X dan Y DULU sebelum meletakkan teks!
xlim([0 max(0.4, max(phi_form)+0.05)]);

% Teks diletakkan tepat di tengah secara vertikal (0.5) agar aman dari judul
text(CO_Phi + 0.015, max(ylim)*0.50, sprintf('PHIE \\geq %.3f', CO_Phi), 'Color', 'r', 'FontWeight', 'bold', 'FontSize', Text_FS, 'BackgroundColor', 'w', 'Margin', 2);
set(gca, 'FontSize', Tick_FS, 'LineWidth', 1.5, 'FontWeight', 'bold', 'GridLineStyle', ':', 'GridAlpha', 0.6);

% ---------------------------------------------------------
% Subplot 2: Crossplot VSH vs PHIE
% ---------------------------------------------------------
subplot(1,3,2);
scatter(vsh_form, phi_form, 45, DEPT(zone), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.8); 
colormap(jet); hold on; grid on;
xline(CO_Vsh, 'r-', 'LineWidth', 3.5); 
yline(CO_Phi, 'r-', 'LineWidth', 3.5);
title('Crossplot VSH vs PHIE', 'FontWeight', 'bold', 'FontSize', Title_FS);
xlabel('Volume Shale (VSH)', 'FontWeight', 'bold', 'FontSize', Label_FS); 
ylabel('Effective Porosity (PHIE)', 'FontWeight', 'bold', 'FontSize', Label_FS);

% KUNCI: Kunci sumbu X dan Y DULU sebelum meletakkan teks!
xlim([0 1.0]); 
ylim([0 max(0.4, max(phi_form)+0.05)]);

% Teks diparkir jauh dari judul dan keramaian titik
text(CO_Vsh + 0.03, max(ylim)*0.50, sprintf('VSH \\leq %.2f', CO_Vsh), 'FontWeight', 'bold', 'FontSize', Text_FS, 'Color', 'r', 'BackgroundColor', 'w', 'Margin', 2);
text(max(xlim)*0.70, CO_Phi + 0.02, sprintf('PHIE \\geq %.3f', CO_Phi), 'FontWeight', 'bold', 'FontSize', Text_FS, 'Color', 'r', 'BackgroundColor', 'w', 'Margin', 2, 'HorizontalAlignment', 'center');

set(gca, 'FontSize', Tick_FS, 'LineWidth', 1.5, 'FontWeight', 'bold', 'GridLineStyle', ':', 'GridAlpha', 0.6);

% ---------------------------------------------------------
% Subplot 3: Crossplot PHIE vs SW
% ---------------------------------------------------------
subplot(1,3,3);
scatter(phi_form, sw_form, 45, DEPT(zone), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.8); 
hold on; grid on;
xline(CO_Phi, 'r-', 'LineWidth', 3.5); 
yline(CO_Sw, 'r-', 'LineWidth', 3.5); 
title('Crossplot PHIE vs SW', 'FontWeight', 'bold', 'FontSize', Title_FS);
xlabel('Effective Porosity (PHIE)', 'FontWeight', 'bold', 'FontSize', Label_FS); 
ylabel('Water Saturation (SW)', 'FontWeight', 'bold', 'FontSize', Label_FS);

% KUNCI: Kunci sumbu X dan Y DULU sebelum meletakkan teks!
xlim([0 max(0.4, max(phi_form)+0.05)]); 
ylim([0 1.0]);

% Teks diparkir di tempat kosong
text(CO_Phi + 0.02, 0.85, sprintf('PHIE \\geq %.3f', CO_Phi), 'FontWeight', 'bold', 'FontSize', Text_FS, 'Color', 'r', 'BackgroundColor', 'w', 'Margin', 2);
text(max(xlim)*0.70, CO_Sw + 0.04, sprintf('SW \\leq %.2f', CO_Sw), 'FontWeight', 'bold', 'FontSize', Text_FS, 'Color', 'r', 'BackgroundColor', 'w', 'Margin', 2, 'HorizontalAlignment', 'center');

set(gca, 'FontSize', Tick_FS, 'LineWidth', 1.5, 'FontWeight', 'bold', 'GridLineStyle', ':', 'GridAlpha', 0.6);

% =========================================================================
% EXPORT GAMBAR CUT-OFF KE RESOLUSI TINGGI
% =========================================================================
exportgraphics(fig_cutoff, 'Analisis_Cutoff_Manual_600DPI.png', 'Resolution', 600);
fprintf('=> Sukses! Gambar Analisis Cut-Off (Manual 600 DPI) telah disimpan.\n');
%% ========================================================================
%                               PAY SUMMARY
% =========================================================================
% Membuat filter untuk MENGECUALIKAN zona data rusak (6700 - 6800 ft)
% CATATAN: Pastikan rentang ini sama persis dengan yang Anda hapus di Excel!
valid_log = ~(DEPT >= 6600 & DEPT <= 6800); 

% Perhitungan ketebalan & rata-rata keseluruhan formasi (Hanya Data Valid)
idx_form = find(DEPT >= Top_Formasi & DEPT <= Base_Formasi & valid_log); 
gross_thick = length(idx_form) * dz;

idx_net_res = find(DEPT >= Top_Formasi & DEPT <= Base_Formasi & Vsh <= CO_Vsh & phi_eff >= CO_Phi & valid_log); 
net_res_thick = length(idx_net_res) * dz;

idx_net_pay = find(DEPT >= Top_Formasi & DEPT <= Base_Formasi & Vsh <= CO_Vsh & phi_eff >= CO_Phi & Sw <= CO_Sw & valid_log); 
net_pay_thick = length(idx_net_pay) * dz;

if net_pay_thick > 0
    % Menghitung rata-rata HANYA pada zona yang lolos Cut-Off
    avg_vsh_pay = mean(Vsh(idx_net_pay), 'omitnan'); 
    avg_phi_pay = mean(phi_eff(idx_net_pay), 'omitnan'); 
    avg_sw_pay  = mean(Sw(idx_net_pay), 'omitnan');
else
    avg_vsh_pay = 0; avg_phi_pay = 0; avg_sw_pay = 0;
end

fprintf('\n=========================================================================\n');
fprintf('                 RESERVOIR PAY SUMMARY (FORMASI KUJUNG)\n');
fprintf('=========================================================================\n');
fprintf(' Interval Kedalaman : %g ft - %g ft\n', Top_Formasi, Base_Formasi);
fprintf(' *Catatan: Interval Washout/Failure dieksklusi (Data Tidak Valid)\n');
fprintf('-------------------------------------------------------------------------\n');
fprintf(' 1. Gross Rock Thickness    : %7.1f ft (Valid)\n', gross_thick);
fprintf(' 2. Net Reservoir Thickness : %7.1f ft\n', net_res_thick);
fprintf(' 3. Net Pay Thickness       : %7.1f ft\n', net_pay_thick);
fprintf(' 4. Net-to-Gross (N/G) Ratio: %7.3f\n', net_pay_thick / max(gross_thick, 0.001));
fprintf(' 5. Average Vshale (Pay)    : %7.5f (%.2f %%)\n', avg_vsh_pay, avg_vsh_pay*100);
fprintf(' 6. Average Porosity (Pay)  : %7.5f (%.2f %%)\n', avg_phi_pay, avg_phi_pay*100);
fprintf(' 7. Average Sw (Pay)        : %7.5f (%.2f %%)\n', avg_sw_pay, avg_sw_pay*100);
fprintf('=========================================================================\n\n');

%% ========================================================================
%             ADDITIONAL CALCULATION: LOGICAL SHADING FLAGS
% ========================================================================
% Pastikan variabel cut-off dinamis (CO_Vsh, CO_Phi, CO_Sw) sudah dihitung sebelumnya.

% 1. Reservoir Flag (Menerapkan valid_log agar tidak ada shading di zona rusak)
netResFlag = (Vsh <= CO_Vsh) & (phi_eff >= CO_Phi) & valid_log;

% 2. Net Pay Flag 
netPayFlag = netResFlag & (Sw <= CO_Sw) & valid_log;

% Menangani nilai NaN agar tidak ikut di-shading
netResFlag(isnan(Vsh) | isnan(phi_eff)) = 0;
netPayFlag(isnan(Vsh) | isnan(phi_eff) | isnan(Sw)) = 0;

%% ========================================================================
%               7. ZONAL SUMMARY & FLUID IDENTIFICATION
% =========================================================================
% 3 Zona Utama
Zona_Kualitatif = [
    5500, 6000;
    6000, 7000;
    7000, 7200
];

fprintf('\n===================================================================================================\n');
fprintf('                        RINCIAN PAY SUMMARY & FLUIDA PER ZONA (LUMPING)\n');
fprintf('===================================================================================================\n');
fprintf(' Zona | Interval (ft) | Gross(ft)*| Net Pay(ft) | Vsh(%%) | Phi(%%) | Sw(%%)  | Sep. D-N | FLUID TYPE\n');
fprintf('---------------------------------------------------------------------------------------------------\n');

for i = 1:size(Zona_Kualitatif, 1)
    z_top = Zona_Kualitatif(i, 1);
    z_bot = Zona_Kualitatif(i, 2);
    
    % Menerapkan filter valid_log di dalam setiap pencarian indeks
    idx_z = find(DEPT >= z_top & DEPT <= z_bot & valid_log);
    z_gross = length(idx_z) * dz; 
    
    idx_z_res = find(DEPT >= z_top & DEPT <= z_bot & Vsh <= CO_Vsh & phi_eff >= CO_Phi & valid_log);
    idx_z_pay = find(DEPT >= z_top & DEPT <= z_bot & Vsh <= CO_Vsh & phi_eff >= CO_Phi & Sw <= CO_Sw & valid_log);
    z_net_pay = length(idx_z_pay) * dz;
    
    fluid_type = 'UNKNOWN';
    z_avg_vsh = 0; z_avg_phi = 0; z_avg_sw = 0; z_crossover = 0;
    
    if isempty(idx_z_res)
        fluid_type = 'TIGHT/SHALE';
    else
        avg_sw_res = mean(Sw(idx_z_res), 'omitnan');
        
        if avg_sw_res >= CO_Sw
            fluid_type = 'WATER';
            z_avg_vsh = mean(Vsh(idx_z_res), 'omitnan') * 100;
            z_avg_phi = mean(phi_eff(idx_z_res), 'omitnan') * 100;
            z_avg_sw  = avg_sw_res * 100;
        else
            if ~isempty(idx_z_pay)
                z_avg_vsh = mean(Vsh(idx_z_pay), 'omitnan') * 100;
                z_avg_phi = mean(phi_eff(idx_z_pay), 'omitnan') * 100;
                z_avg_sw  = mean(Sw(idx_z_pay), 'omitnan') * 100;
                
                % Crossover murni tanpa batasan
                avg_phiD = mean(phi_Dc(idx_z_pay), 'omitnan');
                avg_phiN = mean(phi_Nc(idx_z_pay), 'omitnan');
                z_crossover = avg_phiD - avg_phiN;
                
                batas_gas_minyak = 0.05; % Diubah menjadi 5% agar lebih sensitif terhadap gas effect ringan
                
                if z_crossover >= batas_gas_minyak 
                    fluid_type = 'GAS';
                elseif z_crossover > 0 && z_crossover < batas_gas_minyak
                    fluid_type = 'OIL / WET GAS';
                else
                    fluid_type = 'OIL / TIGHT'; 
                end
            end
        end
    end
    
    fprintf('  %d   | %4d - %4d | %8.1f  | %9.1f   | %5.2f  | %5.2f  | %5.2f  |  %6.3f  | %s\n', ...
            i, z_top, z_bot, z_gross, z_net_pay, z_avg_vsh, z_avg_phi, z_avg_sw, z_crossover, fluid_type);
end
fprintf('---------------------------------------------------------------------------------------------------\n');
fprintf('===================================================================================================\n\n');

%% ========================================================================
%             LITHOLOGY & FLUID VALIDATION CROSSPLOTS (JURNAL SINTA)
% =========================================================================
% Definisi Zona Kualitatif
Zona_Kualitatif = [
    5500, 6000;
    6000, 7000;
    7000, 7200
];
ZoneColors = [
    1.0 0.0 0.0; % Z1 (Merah)
    0.0 0.0 1.0; % Z2 (Biru)
    0.0 0.8 0.0; % Z3 (Hijau)
];
idx_formasi = (DEPT >= Top_Formasi & DEPT <= Base_Formasi & valid_log);

% --- PENGATURAN FONT GLOBAL UNTUK JURNAL ---
Title_FS = 22;   % Ukuran Font Judul Gambar
Label_FS = 18;   % Ukuran Font Sumbu (X dan Y Label)
Tick_FS  = 14;   % Ukuran Font Angka Skala (Tick Marks)
Text_FS  = 14;   % Ukuran Font Teks di dalam plot (Nama Matriks)

% -------------------------------------------------------------------------
% 8.1 DENSITY-NEUTRON CROSSPLOT
% -------------------------------------------------------------------------
fig_dn = figure('Name','Density-Neutron Crossplot','Color','w','Position',[100 50 1000 950]);
hold on; box on; xlim([-0.05 0.45]); ylim([1.9 3.0]);

% UPGRADE JURNAL: Font angka skala diperbesar
set(gca, 'XDir', 'normal', 'YDir', 'reverse', 'FontSize', Tick_FS, 'LineWidth', 2.0, 'FontWeight', 'bold'); 
grid on; set(gca, 'GridLineStyle', ':', 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.6);

% UPGRADE JURNAL: Label diperbesar ke 18, Title ke 22
xlabel('Neutron Porosity (v/v) - Limestone Matrix', 'FontWeight','bold', 'FontSize',Label_FS);
ylabel('Bulk Density, \rho_b (g/cc)', 'FontWeight','bold', 'FontSize',Label_FS);
title({'Density-Neutron Crossplot'}, 'FontWeight','bold', 'FontSize',Title_FS);

phi_ref = 0:0.01:0.50; 
pts_n_ss=[-0.015,0.04,0.10,0.16,0.22,0.28,0.34]; pts_d_ss=[2.65,2.56,2.46,2.36,2.26,2.16,2.06]; 
n_line_ss=interp1(pts_n_ss,pts_n_ss,phi_ref,'linear','extrap'); d_line_ss=interp1(pts_n_ss,pts_d_ss,phi_ref,'pchip','extrap');
n_line_ls=phi_ref; d_line_ls=(1-phi_ref)*2.71+phi_ref*1.0; 
pts_n_dl=[0.01,0.07,0.13,0.19,0.25,0.31,0.38]; pts_d_dl=[2.87,2.76,2.65,2.54,2.43,2.32,2.20];
n_line_dl=interp1(pts_n_dl,pts_n_dl,phi_ref,'linear','extrap'); d_line_dl=interp1(pts_n_dl,pts_d_dl,phi_ref,'pchip','extrap');

% Garis dipertebal
lw=2.5; 
plot(n_line_ss, d_line_ss, 'k-', 'LineWidth', lw); plot(n_line_ls, d_line_ls, 'k-', 'LineWidth', lw); plot(n_line_dl, d_line_dl, 'k-', 'LineWidth', lw); 
text(-0.02,2.65,'Quartz','FontWeight','bold','FontSize',Text_FS,'Rotation',-25); 
text(0.00,2.72,'Limestone','FontWeight','bold','FontSize',Text_FS,'Rotation',-28); 
text(0.02,2.88,'Dolomite','FontWeight','bold','FontSize',Text_FS,'Rotation',-30);

for t = 0:0.05:0.40
    d_ss_t=(1-t)*2.65+t*1.0; n_ss_t=interp1(d_line_ss,n_line_ss,d_ss_t); 
    d_ls_t=(1-t)*2.71+t*1.0; n_ls_t=t; 
    d_dl_t=(1-t)*2.87+t*1.0; n_dl_t=interp1(d_line_dl,n_line_dl,d_dl_t);
    plot([n_ss_t n_dl_t],[d_ss_t d_dl_t],'k:','LineWidth',1.5);
    plot(n_ss_t,d_ss_t,'k_','MarkerSize',10,'LineWidth',2.5); plot(n_ls_t,d_ls_t,'k_','MarkerSize',10,'LineWidth',2.5); plot(n_dl_t,d_dl_t,'k_','MarkerSize',10,'LineWidth',2.5);
    if t>0, text(n_dl_t+0.01,d_dl_t,num2str(t*100),'FontSize',12,'FontWeight','bold','Color','b'); else, text(n_dl_t+0.01,d_dl_t,'0','FontSize',12,'FontWeight','bold','Color','b'); end
end

idx_bg = idx_formasi & ~isnan(NPHI) & ~isnan(RHOB);
scatter(NPHI(idx_bg), RHOB(idx_bg), 25, [0.6 0.6 0.6], 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5);

h_zones = zeros(1, size(Zona_Kualitatif, 1)); 
for i = 1:size(Zona_Kualitatif, 1)
    idx_z = (DEPT >= Zona_Kualitatif(i,1) & DEPT <= Zona_Kualitatif(i,2)) & ~isnan(NPHI) & ~isnan(RHOB);
    h_zones(i) = scatter(NPHI(idx_z), RHOB(idx_z), 55, ZoneColors(i,:), 'filled', 'MarkerEdgeColor','k', 'LineWidth',1.0, 'MarkerFaceAlpha',0.9);
end
legend_names = cell(1, size(Zona_Kualitatif, 1));
for j = 1:size(Zona_Kualitatif, 1)
    legend_names{j} = sprintf('Z%d (%d-%d)', j, Zona_Kualitatif(j,1), Zona_Kualitatif(j,2));
end
legend(h_zones, legend_names, 'Location','southwest', 'FontSize', Text_FS, 'FontWeight','bold');
plot(-0.01, 2.98, 'ko', 'MarkerFaceColor','w', 'LineWidth',2.0, 'MarkerSize',10); text(0.01, 2.98, 'Anhydrite', 'FontSize',Text_FS, 'FontWeight','bold');
plot(1.0, 1.0, 'bo', 'MarkerFaceColor','b', 'MarkerSize',10); text(0.40, 1.95, '\rho_f = 1.0 g/cc', 'FontSize',Text_FS, 'FontWeight','bold', 'Color','b'); hold off;

exportgraphics(fig_dn, '1_Density_Neutron_Crossplot_600DPI.png', 'Resolution', 600);

% -------------------------------------------------------------------------
% 8.2 MATRIX IDENTIFICATION PLOT (MID PLOT)
% -------------------------------------------------------------------------
fig_mid = figure('Name','MID Plot','Color','w','Position',[50 50 1100 1000]); 
axMain = axes('Position', [0.12 0.12 0.82 0.78]); hold(axMain, 'on'); box(axMain, 'on');
xlim(axMain, [35 75]); ylim(axMain, [2.4 3.1]); 
set(axMain, 'YDir', 'reverse', 'FontSize', Tick_FS, 'LineWidth', 2.0, 'FontWeight', 'bold');
set(axMain, 'XTick', 35:5:75, 'YTick', 2.4:0.1:3.1);  
grid(axMain, 'on'); set(axMain, 'XMinorGrid', 'on', 'YMinorGrid', 'on'); 
set(axMain, 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.6, 'GridLineStyle', ':');

rho_f_mid=1.00; dt_f_mid=189;
rho_maa = (RHOB - NPHI.*rho_f_mid) ./ (1 - NPHI);
dt_maa  = (DT - NPHI.*dt_f_mid) ./ (1 - NPHI);
QZ=[55.5, 2.65]; LS=[47.6, 2.71]; DL=[43.5, 2.87]; AN=[50.0, 2.98];
plot(axMain, [QZ(1) LS(1) DL(1) QZ(1)], [QZ(2) LS(2) DL(2) QZ(2)], 'k-', 'LineWidth', 2.5);

ms_main=10; 
plot(axMain, QZ(1), QZ(2), 'ro', 'MarkerFaceColor','w', 'LineWidth',2.0, 'MarkerSize',ms_main); 
plot(axMain, LS(1), LS(2), 'ro', 'MarkerFaceColor','w', 'LineWidth',2.0, 'MarkerSize',ms_main); 
plot(axMain, DL(1), DL(2), 'ro', 'MarkerFaceColor','w', 'LineWidth',2.0, 'MarkerSize',ms_main); 
plot(axMain, AN(1), AN(2), 'ko', 'MarkerFaceColor','w', 'LineWidth',1.5, 'MarkerSize',ms_main);

text(axMain, QZ(1)+1, QZ(2), 'Quartz', 'Color','r', 'FontWeight','bold', 'FontSize',Text_FS); 
text(axMain, LS(1)-5, LS(2), 'Calcite', 'Color','r', 'FontWeight','bold', 'FontSize',Text_FS); 
text(axMain, DL(1)-5, DL(2), 'Dolomite', 'Color','r', 'FontWeight','bold', 'FontSize',Text_FS); 
text(axMain, AN(1)+1, AN(2), 'Anhydrite', 'Color','k', 'FontWeight','bold', 'FontSize',Text_FS);

idx_bg_mid = idx_formasi & ~isnan(rho_maa) & ~isnan(dt_maa);
scatter(axMain, dt_maa(idx_bg_mid), rho_maa(idx_bg_mid), 25, [0.6 0.6 0.6], 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5);
for i = 1:size(Zona_Kualitatif, 1)
    idx_z_mid = (DEPT >= Zona_Kualitatif(i,1) & DEPT <= Zona_Kualitatif(i,2)) & ~isnan(rho_maa) & ~isnan(dt_maa);
    scatter(axMain, dt_maa(idx_z_mid), rho_maa(idx_z_mid), 55, ZoneColors(i,:), 'filled', 'MarkerEdgeColor','k', 'LineWidth',1.0, 'MarkerFaceAlpha',0.9);
end

xlabel(axMain, '\Delta t_{maa}, Apparent Matrix Transit Time (\mus/ft)', 'FontWeight','bold', 'FontSize',Label_FS);
ylabel(axMain, '\rho_{maa}, Apparent Matrix Density (g/cc)', 'FontWeight','bold', 'FontSize',Label_FS);
title(axMain, {'Matrix Identification Plot (MID)'}, 'FontWeight','bold', 'FontSize',Title_FS);
legend(legend_names, 'Location','southwest', 'FontSize', Text_FS, 'FontWeight','bold'); hold(axMain, 'off');

exportgraphics(fig_mid, '2_MID_Plot_600DPI.png', 'Resolution', 600);

% -------------------------------------------------------------------------
% 8.3 SONIC–DENSITY CROSSPLOT (PERBAIKAN SKALA & SPASI ANGKA)
% -------------------------------------------------------------------------
fig_sd = figure('Name','Sonic-Density Crossplot','Color','w','Position',[50 50 1100 1000]);
hold on; box on;

% PERBAIKAN SPASI ANGKA SUMBU (TICK MARKS DIJARANGKAN)
set(gca, 'XLim', [40 120], 'YLim', [1.9 3.0], 'YDir', 'reverse', 'FontSize', Tick_FS, 'LineWidth', 2.0, 'FontWeight', 'bold');
% XTicks diubah dari 2 menjadi 10 agar tidak bertumpuk rapat
set(gca, 'XTick', 40:10:120); 
% YTicks diubah dari 0.02 menjadi 0.1 agar lega
set(gca, 'YTick', 1.90:0.1:3.00); 

% Minor tick tetap ada agar grid kecil tetap rapi
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; set(gca, 'GridLineStyle', '-', 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.6);

xlabel('Interval Transit Time, \Delta t (\mus/ft)', 'FontSize',Label_FS, 'FontWeight','bold');
ylabel('Bulk Density, \rho_b (g/cc)', 'FontSize',Label_FS, 'FontWeight','bold');
title({'Sonic–Density Crossplot'}, 'FontSize',Title_FS, 'FontWeight','bold');

rho_f=1.00; dt_f=189; QZ=[2.65 55.5]; LS=[2.71 47.6]; DL=[2.87 43.5]; AN=[2.98 50.0]; HL=[2.03 67.0];
phi_fine=(0:0.1:45)/100; phi_tick=(0:1:40)/100;  
CalcRho = @(rho_ma, phi) (1 - phi) .* rho_ma + phi .* rho_f;
CalcDt_Wyllie = @(dt_ma, phi) (1 - phi) .* dt_ma + phi .* dt_f;
CalcDt_RHG = @(dt_ma, phi) dt_ma ./ (1 - 1.6 * phi);
Matrices = {'Quartz', QZ; 'Calcite', LS; 'Dolomite', DL}; LabelRot = [45, 50, 52];

for i = 1:size(Matrices,1)
    name=Matrices{i,1}; rho_ma=Matrices{i,2}(1); dt_ma=Matrices{i,2}(2);
    rho_c=CalcRho(rho_ma, phi_fine); dt_w=CalcDt_Wyllie(dt_ma, phi_fine); dt_rhg=CalcDt_RHG(dt_ma, phi_fine);
    plot(dt_w, rho_c, 'r-', 'LineWidth', 2.5); plot(dt_rhg, rho_c, 'k-', 'LineWidth', 3.0); 
    rho_t=CalcRho(rho_ma, phi_tick); dt_wt=CalcDt_Wyllie(dt_ma, phi_tick); dt_rt=CalcDt_RHG(dt_ma, phi_tick);
    for j = 1:length(phi_tick)
        p_val = phi_tick(j)*100; 
        plot(dt_wt(j), rho_t(j), 'r|', 'MarkerSize',10, 'LineWidth',2.0); plot(dt_rt(j), rho_t(j), 'k|', 'MarkerSize',10, 'LineWidth',2.0);
        if mod(p_val, 10) == 0
            text(dt_wt(j)+0.5, rho_t(j)+0.01, num2str(p_val), 'Color','r', 'FontSize',12, 'FontWeight','bold', 'Rotation', LabelRot(i)-5, 'HorizontalAlignment','left');
            text(dt_rt(j)-0.5, rho_t(j)-0.01, num2str(p_val), 'Color','k', 'FontSize',12, 'FontWeight','bold', 'Rotation', LabelRot(i)+5, 'HorizontalAlignment','right');
        end
    end
    plot(dt_ma, rho_ma, 'ko', 'MarkerFaceColor','w', 'MarkerSize', 12, 'LineWidth', 1.5); 
    text(dt_ma-2.0, rho_ma-0.03, name, 'FontWeight','bold', 'FontSize',Text_FS);
end
plot(AN(2), AN(1), 'ko', 'MarkerFaceColor','w', 'MarkerSize',12, 'LineWidth', 1.5); text(AN(2)+2.0, AN(1), 'Anhydrite', 'FontSize',Text_FS, 'FontWeight','bold');
plot(HL(2), HL(1), 'ko', 'MarkerFaceColor','w', 'MarkerSize',12, 'LineWidth', 1.5); text(HL(2)+2.0, HL(1), 'Halite', 'FontSize',Text_FS, 'FontWeight','bold');

idx_bg_sd = idx_formasi & ~isnan(DT) & ~isnan(RHOB);
scatter(DT(idx_bg_sd), RHOB(idx_bg_sd), 25, [0.6 0.6 0.6], 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5);
h_zones_sd = zeros(1, size(Zona_Kualitatif, 1));
for i = 1:size(Zona_Kualitatif, 1)
    idx_z = (DEPT >= Zona_Kualitatif(i,1) & DEPT <= Zona_Kualitatif(i,2)) & ~isnan(DT) & ~isnan(RHOB);
    h_zones_sd(i) = scatter(DT(idx_z), RHOB(idx_z), 55, ZoneColors(i,:), 'filled', 'MarkerEdgeColor','k', 'LineWidth',1.0, 'MarkerFaceAlpha',0.9);
end

h1 = plot(nan,nan, 'k-', 'LineWidth',3.0); h2 = plot(nan,nan, 'r-', 'LineWidth',2.5);
legend([h1 h2 h_zones_sd(1)], {'Empirical (Raymer-Hunt-Gardner)', 'Time Average (Wyllie)', 'Zona Potensi'}, 'Location', 'northeast', 'FontSize',Text_FS, 'FontWeight','bold');
hold off;

exportgraphics(fig_sd, '3_Sonic_Density_Crossplot_600DPI.png', 'Resolution', 600);

% -------------------------------------------------------------------------
% 8.4 NEUTRON-SONIC CROSSPLOT 
% -------------------------------------------------------------------------
fig_ns = figure('Name','Neutron-Sonic Crossplot','Color','w','Position',[50 20 1100 1100]); 
ax = axes('Position',[0.12 0.1 0.8 0.8]); 
hold on; box on; 
set(ax, 'XTick', [], 'YTick', []); 
title({'Neutron-Sonic Crossplot'}, 'FontWeight','bold', 'FontSize',Title_FS);
xlim([-5 45]); ylim([40 120]); 

col_major = [0.0 0.3 0.4]; col_minor = [0.0 0.6 0.7]; 
for x = -5:1:45
    if mod(x,5)==0
        plot([x x], [40 120], 'Color', col_major, 'LineWidth', 2.0); 
        text(x, 38.0, num2str(x), 'Horiz','center', 'FontSize',Tick_FS, 'FontWeight','bold');
    else
        plot([x x], [40 120], 'Color', col_minor, 'LineWidth', 1.0, 'LineStyle', '-'); 
    end
end
for y = 40:2:120
    if mod(y,10)==0
        plot([-5 45], [y y], 'Color', col_major, 'LineWidth', 2.0); 
        text(-5.5, y, num2str(y), 'Horiz','right', 'FontSize',Tick_FS, 'FontWeight','bold');
    else
        plot([-5 45], [y y], 'Color', col_minor, 'LineWidth', 1.0, 'LineStyle', '-'); 
    end
end
for ym = 140:20:380
    y_ft = ym / 3.28084;
    if y_ft >= 40 && y_ft <= 120
        plot([45 45.8], [y_ft y_ft], 'Color','k', 'LineWidth',2.0); 
        text(46.5, y_ft, num2str(ym), 'Horiz','left', 'FontSize',12, 'Color','k', 'FontWeight','bold'); 
    end
end

dt_qz=55.5; dt_ls=47.6; dt_dl=43.5; dt_fl=189; phi=0:0.1:46; phi_pct=phi;
w_qz=(1-phi/100)*dt_qz+(phi/100)*dt_fl; w_ls=(1-phi/100)*dt_ls+(phi/100)*dt_fl; w_dl=(1-phi/100)*dt_dl+(phi/100)*dt_fl;
pts_n=[0,10,20,30,40,45]; pts_d_qz=[55.5,67.5,84.0,106.0,138.0,158.0]; emp_qz=spline(pts_n,pts_d_qz,phi_pct);
pts_d_ls=[47.6,58.5,73.0,93.0,120.0,138.0]; emp_ls=spline(pts_n,pts_d_ls,phi_pct); pts_d_dl=[43.5,53.5,66.5,84.0,108.0,124.0]; emp_dl=spline(pts_n,pts_d_dl,phi_pct);

plot(phi_pct, w_qz, 'r-', 'LineWidth', 2.5); plot(phi_pct, w_ls, 'r-', 'LineWidth', 2.5); plot(phi_pct, w_dl, 'r-', 'LineWidth', 2.5); 
plot(phi_pct, emp_qz, 'k-', 'LineWidth', 3.0); plot(phi_pct, emp_ls, 'k-', 'LineWidth', 3.0); plot(phi_pct, emp_dl, 'k-', 'LineWidth', 3.0);

text(11, 55.5 + 11*1.335, '  Quartz  ', 'Rotation',45, 'FontSize',Text_FS, 'FontWeight','bold', 'BackgroundColor','w', 'EdgeColor','none', 'HorizontalAlignment', 'center'); 
text(16, 47.6 + 16*1.414, '  Calcite  ', 'Rotation',47, 'FontSize',Text_FS, 'FontWeight','bold', 'BackgroundColor','w', 'EdgeColor','none', 'HorizontalAlignment', 'center'); 
text(21, 43.5 + 21*1.455, '  Dolomite  ', 'Rotation',48, 'FontSize',Text_FS, 'FontWeight','bold', 'BackgroundColor','w', 'EdgeColor','none', 'HorizontalAlignment', 'center');

plot(0, dt_qz, 'ko','MarkerFaceColor','w', 'MarkerSize', 10, 'LineWidth', 1.5); 
plot(0, dt_ls, 'ko','MarkerFaceColor','w', 'MarkerSize', 10, 'LineWidth', 1.5); 
plot(0, dt_dl, 'ko','MarkerFaceColor','w', 'MarkerSize', 10, 'LineWidth', 1.5);

rectangle('Position', [1 103 26 15], 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 2.0); 
plot([3 7], [113 113], 'r-', 'LineWidth', 3.0); text(9, 113, 'Time Average', 'FontSize',Text_FS, 'FontWeight','bold'); 
plot([3 7], [107 107], 'k-', 'LineWidth', 3.0); text(9, 107, 'Empirical', 'FontSize',Text_FS, 'FontWeight','bold');

rectangle('Position', [22 43 23 7], 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5); 
text(33.5, 46.5, ['\Delta t_f = ' num2str(dt_fl) ' \mu s/ft'], 'Horiz','center', 'FontSize',Text_FS, 'FontWeight','bold');

NPHI_pct = NPHI * 100; 
idx_bg_ns = idx_formasi & ~isnan(DT) & ~isnan(NPHI);
scatter(NPHI_pct(idx_bg_ns), DT(idx_bg_ns), 25, [0.6 0.6 0.6], 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5);

for i = 1:size(Zona_Kualitatif, 1)
    idx_z = (DEPT >= Zona_Kualitatif(i,1) & DEPT <= Zona_Kualitatif(i,2)) & ~isnan(DT) & ~isnan(NPHI);
    scatter(NPHI_pct(idx_z), DT(idx_z), 55, ZoneColors(i,:), 'filled', 'MarkerEdgeColor','k', 'LineWidth',1.0, 'MarkerFaceAlpha',0.9);
end
hold off;

exportgraphics(fig_ns, '4_Neutron_Sonic_Crossplot_600DPI.png', 'Resolution', 600);
%% ========================================================================
%           8.5 TABEL KESIMPULAN CROSSPLOT (LITHOLOGY SUMMARY TABLE)
% =========================================================================
calc_rho_from_dt = @(dt_log, rho_m, dt_m) ...
    (1 - (0.625 * (1 - dt_m./dt_log))) * rho_m + (0.625 * (1 - dt_m./dt_log)) * 1.0;
fprintf('\n==========================================================================================================================\n');
fprintf('                             TABLE: COMPARISON OF LITHOLOGY ESTIMATION FROM VARIOUS CROSSPLOTS\n');
fprintf('==========================================================================================================================\n');
fprintf('        |                 |                                       Lithology                                      |\n');
fprintf('        |                 |--------------------------------------------------------------------------------------|\n');
fprintf('  Zone  |   Zone Range    |   Neutron-Density   |    Sonic-Density    |    Neutron-Sonic    |      MID Plot      |\n');
fprintf('        |      (ft)       |      Crossplot      |      Crossplot      |      Crossplot      |                    |\n');
fprintf('--------------------------------------------------------------------------------------------------------------------------\n');
for i = 1:size(Zona_Kualitatif, 1)
    z_top = Zona_Kualitatif(i,1);
    z_bot = Zona_Kualitatif(i,2);
    
    % Menggunakan filter valid_log (mengabaikan data washout)
    idx = (DEPT >= z_top & DEPT <= z_bot) & valid_log & ~isnan(RHOB) & ~isnan(NPHI) & ~isnan(DT);
    
    if sum(idx) == 0
        fprintf('   %d    |  %4d - %4d   |       NO DATA       |       NO DATA       |       NO DATA       |      NO DATA       |\n', i, z_top, z_bot);
        continue;
    end
    
    % Titik tengah data per zona
    N_m = mean(NPHI(idx));
    R_m = mean(RHOB(idx));
    D_m = mean(DT(idx));
    
    % --- 1. LITHOLOGY DARI NPHI-RHOB ---
    rho_maa_val = (R_m - N_m*1.0) / (1 - N_m);
    if rho_maa_val < 2.68 
        lit_ND = 'Gas Effect'; 
    elseif rho_maa_val < 2.80
        lit_ND = 'Limestone';          
    else
        lit_ND = 'Dolomite';          
    end
    
    % --- 2. LITHOLOGY DARI DT-RHOB ---
    rho_QZ = calc_rho_from_dt(D_m, 2.65, 55.5);
    rho_LS = calc_rho_from_dt(D_m, 2.71, 47.6);
    rho_DL = calc_rho_from_dt(D_m, 2.87, 43.5);
    [~, min_id_SD] = min(abs(R_m - [rho_QZ, rho_LS, rho_DL]));
    nama_lit_SD = {'Gas Effect', 'Limestone', 'Dolomite'}; 
    if R_m < (rho_QZ - 0.05)
        lit_SD = 'Gas Effect';     
    else
        lit_SD = nama_lit_SD{min_id_SD};
    end
    
    % --- 3. LITHOLOGY DARI NPHI-DT ---
    dt_maa_val = (D_m - N_m*189) / (1 - N_m);
    if dt_maa_val > 51 
        lit_NS = 'Gas Effect'; 
    elseif dt_maa_val > 44
        lit_NS = 'Limestone';
    else
        lit_NS = 'Dolomite';
    end
    
    % --- 4. LITHOLOGY DARI MID PLOT ---
    d1 = sqrt((dt_maa_val - 55.5)^2 + (50*(rho_maa_val - 2.65))^2);
    d2 = sqrt((dt_maa_val - 47.6)^2 + (50*(rho_maa_val - 2.71))^2);
    d3 = sqrt((dt_maa_val - 43.5)^2 + (50*(rho_maa_val - 2.87))^2);
    [~, min_id_MID] = min([d1, d2, d3]);
    if rho_maa_val < 2.60 || dt_maa_val > 58
        lit_MID = 'Gas Effect';
    else
        lit_MID = nama_lit_SD{min_id_MID}; 
    end
    
    % Print tabel
    fprintf('   %d    |  %4d - %4d   | %-19s | %-19s | %-19s | %-18s |\n', ...
            i, z_top, z_bot, lit_ND, lit_SD, lit_NS, lit_MID);
end
fprintf('==========================================================================================================================\n\n');
%% ========================================================================
%             SAFETY CLAMP (PEMBATAS FISIKA UNTUK HASIL AKHIR)
% =========================================================================
% Membuat variabel final agar plot dan Excel tidak meledak dari skala 0-1
Vsh_Final = Vsh; 
Vsh_Final(Vsh_Final < 0) = 0; 
Vsh_Final(Vsh_Final > 1) = 1; 

Sw_Final = Sw;
Sw_Final(phi_eff <= 0.01) = 1; 
Sw_Final(Sw_Final > 1) = 1;    
Sw_Final(Sw_Final < 0) = 0;    

%% ========================================================================
%        9. WELL LOG VISUALIZATION – JOURNAL PUBLICATION READY
% =========================================================================
% 9.1 Definisi Logika Flags & Lithology
valid_log = ~(DEPT >= 6600 & DEPT <= 6800);
LithCode = NaN(size(DEPT));
LithCode(Vsh_Final > CO_Vsh) = 1;                                      
LithCode(Vsh_Final <= CO_Vsh & phi_eff < CO_Phi) = 2;                  
LithCode(Vsh_Final <= CO_Vsh & phi_eff >= CO_Phi & Sw_Final > CO_Sw) = 3;    
LithCode(Vsh_Final <= CO_Vsh & phi_eff >= CO_Phi & Sw_Final <= CO_Sw) = 4;   
GR_cutoff = GR_clean + (CO_Vsh * (GR_shale - GR_clean));

netResFlag = (Vsh_Final <= CO_Vsh) & (phi_eff >= CO_Phi) & valid_log;
netPayFlag = netResFlag & (Sw_Final <= CO_Sw) & valid_log;
netResFlag(isnan(Vsh_Final) | isnan(phi_eff)) = 0;
netPayFlag(isnan(Vsh_Final) | isnan(phi_eff) | isnan(Sw_Final)) = 0;

% --- PENGATURAN UMUM VISUALISASI JURNAL ---
font_name = 'Arial';
grid_color = [0.75 0.75 0.75]; 
axis_size = 9;
margin_left = 0.08; margin_right = 0.05; margin_bottom = 0.05; 
h_header = 0.10; h_track = 0.82; gap = 0.005;
y_header = margin_bottom + h_track; y_track = margin_bottom;

disp('=> Sedang menggambar dan mengekspor visualisasi ke resolusi 600 DPI...');

%% ========================================================================
%  FIGURE 1: TRIPLE COMBO LOG (GR, RESISTIVITY, DENSITY-NEUTRON)
% =========================================================================
fig_tc = figure('Name','Triple Combo Log','Color','w','Position',[50 50 1000 800], 'NumberTitle','off');
total_width_tc = 1.0 - margin_left - margin_right;
w_tc = (total_width_tc - (2 * gap)) / 3; 

lbl_tc = {
    {'GAMMA RAY', '(API)', '0', '150'}, ...          
    {'RESISTIVITY', '(OHM.M)', '0.2', '2000'}, ...   
    {'RHOB (Red) / NPHI (Blue)', '(G/CC) / (V/V)', '1.7 / 0.6', '2.7 / 0'} ...
};

ax_tc = gobjects(1, 3);
curr_x = margin_left;
for i = 1:3
    % Header
    axes('Position',[curr_x, y_header, w_tc, h_header]); axis off; 
    rectangle('Position',[0 0 1 1],'LineWidth',1); 
    text(0.5, 0.75, lbl_tc{i}{1}, 'Horiz','center','FontWeight','bold','FontSize',10);
    text(0.5, 0.45, lbl_tc{i}{2}, 'Horiz','center','FontSize',9,'FontAngle','italic');
    text(0.02, 0.1, lbl_tc{i}{3}, 'Horiz','left','FontSize',9,'FontWeight','bold');
    text(0.98, 0.1, lbl_tc{i}{4}, 'Horiz','right','FontSize',9,'FontWeight','bold');
    
    % Track
    ax_tc(i) = axes('Position', [curr_x, y_track, w_tc, h_track]);
    set(ax_tc(i), 'FontName', font_name, 'FontSize', axis_size, 'Box', 'on', 'LineWidth', 1.0, ...
        'TickDir', 'out', 'XColor', 'k', 'YColor', 'k', 'YDir', 'reverse', 'XAxisLocation', 'bottom', ...
        'XTickLabel', [], 'XGrid', 'on', 'YGrid', 'on', 'GridColor', grid_color, 'GridAlpha', 0.8, ...
        'XMinorGrid', 'on', 'YMinorGrid', 'on', 'MinorGridColor', grid_color, 'MinorGridAlpha', 0.4);
    hold(ax_tc(i), 'on');
    curr_x = curr_x + w_tc + gap;
end

% Plot GR
axes(ax_tc(1)); 
hold on;

% 1. Logika Shading 
for i = 1:length(DEPT)-1
    if isnan(GR(i)) || isnan(GR(i+1)); continue; end
    
    if GR(i) <= GR_cutoff
        col = [1 1 0]; % Kuning cerah (Zona Permeabel)
    else
        col = [0 0.6 0]; % Hijau (Zona Shale)
    end
    patch([GR(i) GR_cutoff GR_cutoff GR(i)], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], col, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
end

% 2. Plot Garis Log dan Garis Cut-off
plot(GR, DEPT, 'g', 'LineWidth', 1.5);
xline(GR_cutoff, 'r--', 'LineWidth', 1.3);

% Pengaturan Sumbu (TANPA MENGUBAH GRID)
xlim([0 150]); xticks([0 50 100 150]); 
ylabel('Depth (ft)', 'FontWeight', 'bold', 'FontSize', 11);

% -----------------------------------------------------------------
% Plot Resistivity
axes(ax_tc(2)); LLD_plot = LLD; LLS_plot = LLS; LLD_plot(LLD_plot<=0)=NaN; LLS_plot(LLS_plot<=0)=NaN; 
semilogx(LLD_plot, DEPT, 'r', 'LineWidth', 1.5); semilogx(LLS_plot, DEPT, 'm', 'LineWidth', 1.0);
set(gca,'XScale','log'); xlim([0.2 2000]); xticks([0.2 1 10 100 1000 2000]); set(gca,'YTickLabel',[]);

% Plot Density-Neutron
axes(ax_tc(3)); NPHI_scaled = 2.7 - (NPHI .* (1.0/0.6)); separation = NPHI_scaled - RHOB;
for i = 1:length(DEPT)-1
    if isnan(RHOB(i)) || isnan(NPHI_scaled(i)); continue; end
    if separation(i) > 0, patch([RHOB(i) NPHI_scaled(i) NPHI_scaled(i) RHOB(i)], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], [1 1 0], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    else, patch([RHOB(i) NPHI_scaled(i) NPHI_scaled(i) RHOB(i)], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], [0 0.6 0], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); end
end
plot(RHOB, DEPT, 'r', 'LineWidth', 1.5); plot(NPHI_scaled, DEPT, 'b--', 'LineWidth', 1.5); 
xlim([1.7 2.7]); xticks([1.7 1.95 2.2 2.45 2.7]); set(gca,'YTickLabel',[]);

% Link & Set Kujung Limit
linkaxes(ax_tc,'y'); set(ax_tc(1), 'YLim', [Top_Formasi Base_Formasi]);

% EXPORT TRIPLE COMBO
exportgraphics(fig_tc, 'Triple_Combo_Kujung_600DPI.png', 'Resolution', 600);

%% ========================================================================
%  FIGURE 2: PETROPHYSICAL LOG (VSH, PHI, SW) - CLEAN VERSION (NO CUT-OFF LINES)
% =========================================================================
fig_petro = figure('Name','Petrophysical Log','Color','w','Position',[100 100 1000 850], 'NumberTitle','off');

% Penyesuaian Lebar (Lithology Dihapus, Lebar dibagi 3 Rata)
margin_left = 0.08; margin_right = 0.05; margin_bottom = 0.05; 
h_header = 0.12; h_track = 0.82; gap = 0.005;
y_header = margin_bottom + h_track; y_track = margin_bottom;
total_width_petro = 1.0 - margin_left - margin_right;
w_pet = (total_width_petro - (2 * gap)) / 3;

lbl_pet = {
    {'VSHALE', '(%)', '0', '1'}, ...               
    {'EFFECTIVE POROSITY', '(%)', '0', '1'}, ...           
    {'WATER SATURATION', '(%)', '0', '1'}                
};

% Track Vsh, Phi, Sw
ax_pet = gobjects(1, 3); curr_x = margin_left;
for i = 1:3
    % 1. Menggambar Header (Ukuran Font Diperbesar)
    axes('Position',[curr_x, y_header, w_pet, h_header]); axis off; 
    rectangle('Position',[0 0 1 1],'LineWidth',1); 
    
    % Judul Track (Misal diperbesar menjadi 14)
    text(0.5, 0.70, lbl_pet{i}{1}, 'Horiz','center','FontWeight','bold','FontSize',14);
    
    % Satuan (%) -> INI YANG ANDA CARI (Diperbesar menjadi 15)
    text(0.5, 0.30, lbl_pet{i}{2}, 'Horiz','center','FontSize',15,'FontWeight','bold');
    
    % Skala 0 dan 1 (Diperbesar menjadi 12)
    text(0.02, 0.1, lbl_pet{i}{3}, 'Horiz','left','FontSize',12,'FontWeight','bold');
    text(0.98, 0.1, lbl_pet{i}{4}, 'Horiz','right','FontSize',12,'FontWeight','bold');
    
    % 2. Menggambar Area Plotting
    ax_pet(i) = axes('Position', [curr_x, y_track, w_pet, h_track]);
    set(ax_pet(i), 'FontName', font_name, 'FontSize', axis_size, 'Box', 'on', 'LineWidth', 1.0, ...
        'TickDir', 'out', 'XColor', 'k', 'YColor', 'k', 'YDir', 'reverse', 'XAxisLocation', 'bottom', ...
        'XTickLabel', [], 'XGrid', 'on', 'YGrid', 'on', 'GridColor', grid_color, 'GridAlpha', 0.8, ...
        'XMinorGrid', 'on', 'YMinorGrid', 'on', 'MinorGridColor', grid_color, 'MinorGridAlpha', 0.4);
    
    % Hanya tampilkan label angka kedalaman di track paling kiri (VSH)
    if i == 1
        ylabel(ax_pet(i), 'Depth (ft)', 'FontWeight', 'bold', 'FontSize', 11);
    else
        set(ax_pet(i), 'YTickLabel', []);
    end
    
    hold(ax_pet(i), 'on'); 
    curr_x = curr_x + w_pet + gap;
end

% ---------------------------------------------------------
% Plot Track 1: Vshale (Sesuai Gambar: Arsiran Abu-abu)
% ---------------------------------------------------------
axes(ax_pet(1)); 
% Arsir terlebih dahulu agar garis log berada di depan arsiran
for i = 1:length(DEPT)-1
    if isnan(Vsh_Final(i)) || isnan(Vsh_Final(i+1)); continue; end
    patch([Vsh_Final(i) 1 1 Vsh_Final(i)], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], [0.85 0.85 0.85], 'FaceAlpha', 1.0, 'EdgeColor', 'none'); 
end
plot(Vsh_Final, DEPT, 'Color', [0.4 0.3 0.1], 'LineWidth', 1.5); % Garis coklat
xlim([0 1]); xticks([0 0.25 0.5 0.75 1]);

% MEMUNCULKAN GRID: Paksa garis kotak-kotak tampil di lapisan paling atas
set(ax_pet(1), 'Layer', 'top');

% ---------------------------------------------------------
% Plot Track 2: Porosity (Sesuai Gambar: Tanpa Arsiran, Cyan)
% ---------------------------------------------------------
axes(ax_pet(2)); 
plot(phi_eff, DEPT, 'c', 'LineWidth', 1.5);
% xline(CO_Phi, 'r--', 'LineWidth', 1.5); % DIMATIKAN
xlim([0 1]); xticks([0 0.25 0.5 0.75 1]); set(gca, 'XDir', 'normal');

% ---------------------------------------------------------
% Plot Track 3: Water Saturation (Sesuai Gambar: Tanpa Arsiran, Kuning)
% ---------------------------------------------------------
axes(ax_pet(3)); 
plot(Sw_Final, DEPT, 'y', 'LineWidth', 1.5);
% xline(CO_Sw, 'r--', 'LineWidth', 1.5); % DIMATIKAN
xlim([0 1]); xticks([0 0.25 0.5 0.75 1]);

% ---------------------------------------------------------
% Link Sumbu Kedalaman & Set Skala
% ---------------------------------------------------------
linkaxes(ax_pet,'y'); 
set(ax_pet(1), 'YLim', [Top_Formasi Base_Formasi]);

% =========================================================================
% EXPORT PETROPHYSICS KE RESOLUSI TINGGI
% =========================================================================
exportgraphics(fig_petro, 'Petrophysics_Kujung_Clean_600DPI.png', 'Resolution', 600);
fprintf('=> Sukses! Gambar Petrophysics Log (Tanpa Litho & Tanpa Cut-off) High-Res telah disimpan.\n');

%% ========================================================================
%  FIGURE 3: NET RESERVOIR & NET PAY IDENTIFICATION (WITH CUT-OFF & FLAGS)
% =========================================================================
fig_pay = figure('Name','Net Pay & Reservoir Log','Color','w','Position',[200 150 1000 850], 'NumberTitle','off');

% Penyesuaian Lebar
margin_left = 0.08; margin_right = 0.05; margin_bottom = 0.05; 
h_header = 0.12; h_track = 0.82; gap = 0.005;
y_header = margin_bottom + h_track; y_track = margin_bottom;
total_width_petro = 1.0 - margin_left - margin_right;
w_pet = (total_width_petro - (2 * gap)) / 3;

lbl_pet = {
    {'VSHALE', '(%)', '0', '1'}, ...               
    {'EFFECTIVE POROSITY', '(%)', '0', '1'}, ...           
    {'WATER SATURATION', '(%)', '0', '1'}                
};

% Track Vsh, Phi, Sw
ax_pay = gobjects(1, 3); curr_x = margin_left;
for i = 1:3
    % 1. Menggambar Header (Ukuran Font Diperbesar Sesuai Permintaan)
    axes('Position',[curr_x, y_header, w_pet, h_header]); axis off; 
    rectangle('Position',[0 0 1 1],'LineWidth',1); 
    
    % Judul Track
    text(0.5, 0.70, lbl_pet{i}{1}, 'Horiz','center','FontWeight','bold','FontSize',14);
    % Satuan (%)
    text(0.5, 0.30, lbl_pet{i}{2}, 'Horiz','center','FontSize',15,'FontWeight','bold');
    % Skala 0 dan 1
    text(0.02, 0.1, lbl_pet{i}{3}, 'Horiz','left','FontSize',12,'FontWeight','bold');
    text(0.98, 0.1, lbl_pet{i}{4}, 'Horiz','right','FontSize',12,'FontWeight','bold');
    
    % 2. Menggambar Area Plotting
    ax_pay(i) = axes('Position', [curr_x, y_track, w_pet, h_track]);
    set(ax_pay(i), 'FontName', font_name, 'FontSize', axis_size, 'Box', 'on', 'LineWidth', 1.0, ...
        'TickDir', 'out', 'XColor', 'k', 'YColor', 'k', 'YDir', 'reverse', 'XAxisLocation', 'bottom', ...
        'XTickLabel', [], 'XGrid', 'on', 'YGrid', 'on', 'GridColor', grid_color, 'GridAlpha', 0.8, ...
        'XMinorGrid', 'on', 'YMinorGrid', 'on', 'MinorGridColor', grid_color, 'MinorGridAlpha', 0.4);
    
    if i == 1
        ylabel(ax_pay(i), 'Depth (ft)', 'FontWeight', 'bold', 'FontSize', 12);
    else
        set(ax_pay(i), 'YTickLabel', []);
    end
    
    hold(ax_pay(i), 'on'); 
    curr_x = curr_x + w_pet + gap;
end

% ---------------------------------------------------------
% Plot Track 1: Vshale
% ---------------------------------------------------------
axes(ax_pay(1)); 
for i = 1:length(DEPT)-1
    if isnan(Vsh_Final(i)) || isnan(Vsh_Final(i+1)); continue; end
    patch([Vsh_Final(i) 1 1 Vsh_Final(i)], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], [0.85 0.85 0.85], 'FaceAlpha', 1.0, 'EdgeColor', 'none'); 
end
plot(Vsh_Final, DEPT, 'Color', [0.4 0.3 0.1], 'LineWidth', 1.5); 
xline(CO_Vsh, 'r--', 'LineWidth', 1.8); % Cut-off Dinyalakan & Dipertebal
xlim([0 1]); xticks([0 0.25 0.5 0.75 1]);
set(ax_pay(1), 'Layer', 'top'); % Menjaga Grid tetap di depan arsiran abu-abu

% ---------------------------------------------------------
% Plot Track 2: Porosity
% ---------------------------------------------------------
axes(ax_pay(2)); 
plot(phi_eff, DEPT, 'c', 'LineWidth', 1.5);
xline(CO_Phi, 'r--', 'LineWidth', 1.8); % Cut-off Dinyalakan & Dipertebal
xlim([0 1]); xticks([0 0.25 0.5 0.75 1]); set(gca, 'XDir', 'normal');
set(ax_pay(2), 'Layer', 'top');

% ---------------------------------------------------------
% Plot Track 3: Water Saturation (Figure 3 - Dengan Cut-off)
% ---------------------------------------------------------
axes(ax_pay(3)); 
plot(Sw_Final, DEPT, 'y', 'LineWidth', 1.5);

% MENAMBAHKAN KEMBALI GARIS CUT-OFF SW
xline(CO_Sw, 'r--', 'LineWidth', 1.8); 

% Label nilai cut-off agar lebih informatif bagi reviewer
text(CO_Sw + 0.02, Top_Formasi + 50, ['CO Sw: ', num2str(CO_Sw, '%.2f')], ...
    'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

xlim([0 1]); xticks([0 0.25 0.5 0.75 1]);
set(ax_pay(3), 'Layer', 'top');
% ---------------------------------------------------------
% Shading Flag (Edge Banding): Net Reservoir & Net Pay
% ---------------------------------------------------------
for i = 1:3
    axes(ax_pay(i)); hold on; 
    xl = xlim(); 
    % Lebar pita diperbesar menjadi 7.5% agar lebih mencolok untuk jurnal
    strip_w = (xl(2) - xl(1)) * 0.075; 
    
    % Hijau (Net Reservoir) di sisi kiri, Merah (Net Pay) di sisi kanan
    func_draw_patches_fast(ax_pay(i), DEPT, netResFlag, [xl(1), xl(1)+strip_w], [0.0 0.8 0.0], 0.8);
    func_draw_patches_fast(ax_pay(i), DEPT, netPayFlag, [xl(2)-strip_w, xl(2)], [1.0 0.0 0.0], 0.8);
end

% ---------------------------------------------------------
% Link Sumbu Kedalaman & Set Skala
% ---------------------------------------------------------
linkaxes(ax_pay,'y'); 
set(ax_pay(1), 'YLim', [Top_Formasi Base_Formasi]);

% =========================================================================
% EXPORT FIGURE 3 KE RESOLUSI TINGGI
% =========================================================================
exportgraphics(fig_pay, 'Net_Pay_Identification_600DPI.png', 'Resolution', 600);
fprintf('=> Sukses! Gambar Figure 3 (Net Pay & Cut-off) High-Res telah disimpan.\n');
%% =========================================================================
% FUNCTION BLOCKS (WAJIB DI BARIS PALING BAWAH SCRIPT)
% =========================================================================

function update_view(src, ~)
    fig = get(src, 'Parent'); data = get(fig, 'UserData');
    val = get(data.menu_scale, 'Value');
    sl_raw = get(data.slider, 'Value'); 
    pct = (sl_raw - data.min_depth) / data.dept_range;
    
    switch val
        case 1, sp = data.dept_range; 
        case 2, sp = data.Base_Formasi - data.Top_Formasi;
        case 3, sp = 1000; 
        case 4, sp = 500; 
        case 5, sp = 200; 
        case 6, sp = 100; 
    end
    
    av_scroll = data.dept_range - sp; if av_scroll<0, av_scroll=0; end
    
    if val == 2 
        curr_top = data.Top_Formasi;
        if av_scroll > 0
            pct_new = 1 - ((curr_top - data.min_depth) / av_scroll);
            set(data.slider, 'Value', data.min_depth + pct_new * data.dept_range);
        end
    else
        curr_top = data.min_depth + (1 - pct) * av_scroll;
    end
    
    if curr_top < data.min_depth, curr_top = data.min_depth; end
    if curr_top + sp > data.max_depth
        if sp ~= data.dept_range, curr_top = data.max_depth - sp; else, curr_top = data.min_depth; end
    end
    
    new_lim = [curr_top, curr_top + sp];
    for i=1:length(data.axes_handles), set(data.axes_handles(i), 'YLim', new_lim); end
    set(data.lbl_depth, 'String', sprintf('%.0f', curr_top));
    
    if sp < data.dept_range
        step = (sp/4) / av_scroll; if step>1, step=1; end
        set(data.slider, 'SliderStep', [0.01, step], 'Enable', 'on');
    else
        set(data.slider, 'Enable', 'off'); 
    end
end

function func_draw_patches_fast(ax, depth, flag, x_bounds, color, alpha)
    dFlag = diff([0; flag(:); 0]);
    start_indices = find(dFlag == 1);
    end_indices = find(dFlag == -1) - 1; 
    
    if ~isempty(start_indices)
        num_intervals = length(start_indices);
        
        X_mat = zeros(4, num_intervals);
        Y_mat = zeros(4, num_intervals);
        
        for idx = 1:num_intervals
            z_start = depth(start_indices(idx));
            z_end = depth(end_indices(idx));
            
            X_mat(:, idx) = [x_bounds(1); x_bounds(2); x_bounds(2); x_bounds(1)];
            Y_mat(:, idx) = [z_start; z_start; z_end; z_end];
        end
        
        p = patch('XData', X_mat, 'YData', Y_mat, 'FaceColor', color, ...
                  'FaceAlpha', alpha, 'EdgeColor', 'none', ...
                  'HandleVisibility', 'off', 'Parent', ax);
                  
        uistack(p, 'bottom');
    end
end

%% ========================================================================
%                     10. EXPORT HASIL KE EXCEL
% =========================================================================
fprintf('\nSedang mengekspor data ke Excel, mohon tunggu...\n');

% 1. Memastikan semua variabel berbentuk kolom (Column Vectors)
Depth_ft = DEPT(:);
GR_api   = GR(:);
Res_Deep = LLD(:);
Rho_bulk = RHOB(:);
Nphi_v   = NPHI(:);
Dt_usft  = DT(:);

% Parameter Hasil Hitungan
Vshale_v     = Vsh_Final(:);
Porosity_Eff = phi_eff(:);
Water_Sat    = Sw_Final(:);

% Flags & Lithology (Diubah jadi angka 1 atau 0 agar mudah difilter di Excel)
Net_Res_Flag = double(netResFlag(:)); 
Net_Pay_Flag = double(netPayFlag(:));
Lith_Code    = LithCode(:); % 1=Shale, 2=Tight, 3=Water, 4=Pay

% KUNCI FILTERING EXCEL: Menambahkan status validitas data
Valid_Data   = double(valid_log(:)); % 1=Valid, 0=Rusak/Washout

% 2. Membuat Tabel MATLAB
Tabel_Petrofisika = table(Depth_ft, GR_api, Res_Deep, Rho_bulk, Nphi_v, Dt_usft, ...
                          Vshale_v, Porosity_Eff, Water_Sat, ...
                          Lith_Code, Net_Res_Flag, Net_Pay_Flag, Valid_Data);

% 3. Menyimpan Tabel ke File Excel
% (.xlsx)
nama_file_excel = 'Hasil_Petrofisika_Formasi_Kujung.xlsx';
writetable(Tabel_Petrofisika, nama_file_excel, 'Sheet', 'Data_Per_Depth');

fprintf('=> SUKSES! Seluruh data dari kedalaman %g ft - %g ft telah disimpan.\n', min(DEPT), max(DEPT));
fprintf('=> Buka Excel Anda, lalu filter kolom "Valid_Data" = 1 untuk menghitung rata-rata.\n');
fprintf('=> Silakan cek folder MATLAB Anda untuk file: %s\n', nama_file_excel);
fprintf('==========================================================================================================================\n\n');


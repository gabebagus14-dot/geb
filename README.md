%% ========================================================================
%               PETROPHYSICAL ANALYSIS - FINAL FIXED VERSION
%  ========================================================================
clc; clear; close all;

%% ========================================================================
%                     READ LAS FILE
% =========================================================================
filename = 'WCL0003714.LAS';   
fid = fopen(filename,'r');
if fid < 0
    error('LAS file not found');
end
line = '';
while ischar(line)
    line = fgetl(fid);
    if contains(upper(line),'~A')
        break
    end
end
ascii_data = [];
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    nums = str2num(line); 
    if ~isempty(nums)
        ascii_data = [ascii_data; nums];
    end
end
fclose(fid);
NULL = -999.25;
ascii_data(ascii_data == NULL) = NaN;
DEPT = ascii_data(:,1);
GR   = ascii_data(:,2);
CALI = ascii_data(:,3);
LLD  = ascii_data(:,4);
LLS  = ascii_data(:,5);
DT   = ascii_data(:,7);
RHOB = ascii_data(:,8);
NPHI = ascii_data(:,9);

%% ========================================================================
%              QUICKLOOK RESERVOIR IDENTIFICATION
% =========================================================================
GR_cutoff = mean(GR(~isnan(GR)));
flag_GR = GR <= GR_cutoff;
Rt_cutoff = prctile(LLD(~isnan(LLD)),60);
flag_Rt = LLD >= Rt_cutoff;
flag_cross = (NPHI <= 0.25) & (RHOB <= 2.35);
flag_reservoir = flag_GR & flag_Rt & flag_cross;
dz = median(diff(DEPT));
flag = flag_reservoir & ~isnan(DEPT);
idx  = find(diff([0; flag; 0]) ~= 0);
TopDepth  = [];
BaseDepth = [];
Thickness = [];
for i = 1:2:length(idx)
    top  = DEPT(idx(i));
    base = DEPT(idx(i+1)-1);
    TopDepth  = [TopDepth; top];
    BaseDepth = [BaseDepth; base];
    Thickness = [Thickness; base - top + dz];
end
ReservoirIntervalTable = table(TopDepth, BaseDepth, Thickness);

%% ========================================================================
%                PETROPHYSICAL CALCULATION 
% =========================================================================
Vsh     = NaN(size(DEPT));
phi_eff = NaN(size(DEPT));
Sw      = NaN(size(DEPT));
k       = NaN(size(DEPT)); % Inisialisasi K Permeability
zone = flag_reservoir & ~isnan(GR) & ~isnan(RHOB) & ...
       ~isnan(NPHI) & ~isnan(LLD); 

%% ========================================================================
%                  VOLUME SHALE CALCULATION
%                    (LARIONOV TERTIARY)
% =========================================================================
GR_valid = GR(~isnan(GR));
GR_clean = prctile(GR_valid, 5);   
GR_shale = prctile(GR_valid, 95);  
IGR = (GR - GR_clean) ./ (GR_shale - GR_clean);
Vsh = 0.083 .* (2.^(3.7 .* IGR) - 1);
Vsh(Vsh < 0) = 0;
Vsh(Vsh > 1) = 1;
if exist('flag_reservoir', 'var')
    Vsh(~flag_reservoir) = NaN; 
else
    warning('Flag Reservoir belum didefinisikan. Vsh ditampilkan full depth.');
end

%% ========================================================================
%                  POROSITY CALCULATION
% =========================================================================
rho_ma = 2.65;  
rho_f  = 1.00;  
idx_shale = Vsh > 0.8 & ~isnan(RHOB) & ~isnan(NPHI);
if sum(idx_shale) > 5
    rho_bsh = median(RHOB(idx_shale));
    phi_Nsh = median(NPHI(idx_shale));
else
    rho_bsh = 2.5; 
    phi_Nsh = 0.35;
end
phi_Dsh = (rho_ma - rho_bsh) / (rho_ma - rho_f);
phi_D = (rho_ma - RHOB) ./ (rho_ma - rho_f);
phi_Dc = phi_D - (Vsh .* phi_Dsh);
phi_Nc = NPHI - (Vsh .* phi_Nsh);
phi_eff_calc = sqrt((phi_Nc.^2 + phi_Dc.^2) ./ 2);
phi_eff(zone) = phi_eff_calc(zone);
phi_eff(phi_eff < 0) = 0;
phi_eff(phi_eff > 0.45) = 0.45; 

%% ========================================================================
%                  WATER SATURATION CALCULATION
%                         (INDONESIA)
% =========================================================================
a = 0.62; 
m = 2.15; 
n = 2;
Rw = 0.08; 
if sum(idx_shale) > 5
    Rsh = median(LLD(idx_shale));
else
    Rsh = 2.0;
end
Rt = LLD(zone);
Vsh_z = Vsh(zone);
Phi_z = phi_eff(zone);
term_shale = (Vsh_z .^ (1 - 0.5 * Vsh_z)) ./ sqrt(Rsh);
term_form = sqrt((Phi_z .^ m) ./ (a * Rw));
Sw_calc = ( (1 ./ sqrt(Rt)) ./ (term_shale + term_form) ) .^ (2 / n);
Sw(zone) = Sw_calc;
Sw(Sw < 0) = 0;
Sw(Sw > 1) = 1;

%% ========================================================================
%                  PERMEABILITY ESTIMATION (TIMUR)
% ========================================================================
k = NaN(size(DEPT)); % Pastikan ukuran sama dengan DEPT
idx_k = ~isnan(phi_eff) & ~isnan(Sw) & phi_eff > 0 & Sw > 0;
if any(idx_k)
    k(idx_k) = 8581 .* (phi_eff(idx_k).^4.4) ./ (Sw(idx_k).^2);
end
% Clamp agar aman di skala log
k(k < 0.01)   = 0.01;
k(k > 10000) = 10000;

%% ========================================================================
%                  SONIC–DENSITY CROSSPLOT
% =========================================================================
if ~exist('DT','var') || ~exist('DEPT','var'), error('Data log tidak ditemukan!'); end
idx = flag_reservoir & ~isnan(DT) & ~isnan(RHOB) & ~isnan(DEPT);
DTv   = DT(idx);        % us/ft
RHOBv = RHOB(idx);      % g/cc
DEPv  = DEPT(idx);      % ft (Depth untuk warna)
rho_f = 1.00;     
dt_f  = 189;      
QZ = [2.65 55.5];  % Quartz
LS = [2.71 47.6];  % Calcite
DL = [2.87 43.5];  % Dolomite
AN = [2.98 50.0];  % Anhydrite
HL = [2.03 67.0];  % Halite
% Porosity vectors for calculations
phi_fine = (0:0.1:45)/100; % Untuk menggambar garis kurva yang halus
phi_tick = (0:1:40)/100;   % Untuk tick mark setiap 1%

%% 3. FUNCTIONS
% A. DENSITY (Linear)
CalcRho = @(rho_ma, phi) (1 - phi) .* rho_ma + phi .* rho_f;
% B. WYLLIE TIME AVERAGE (Linear - Red Line)
CalcDt_Wyllie = @(dt_ma, phi) (1 - phi) .* dt_ma + phi .* dt_f;
% C. RAYMER-HUNT-GARDNER (Empirical - Black Curved Line)
CalcDt_RHG = @(dt_ma, phi) dt_ma ./ (1 - 1.6 * phi);

%% 4. PLOTTING SETUP
figure('Color','w','Position',[50 50 900 800]);
hold on; box on;
% --- DETAILED GRID SETUP (KOTAK-KOTAK MENDETAIL) ---
% Mengatur tick secara manual agar sangat rapat seperti kertas grafik
set(gca, 'XLim', [40 120], 'YLim', [1.9 3.0], 'YDir', 'reverse');
set(gca, 'XTick', 40:2:120);       % Garis vertikal setiap 2 us/ft
set(gca, 'YTick', 1.90:0.02:3.00); % Garis horizontal setiap 0.02 g/cc
grid on;
set(gca, 'GridLineStyle', '-', 'GridColor', [0.4 0.4 0.4], 'GridAlpha', 0.6);
set(gca, 'FontSize', 10);

%% 5. PLOT DATA (COLOR BY DEPTH)
% Menggunakan DEPv untuk warna
scatter(DTv, RHOBv, 30, DEPv, 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.75);
% Setup Colorbar untuk Depth
colormap(flipud(jet)); % Membalik warna (Merah=Dalam, Biru=Dangkal)
c = colorbar;
c.Label.String = 'Depth (ft)';
c.Label.FontSize = 11;
c.Label.FontWeight = 'bold';

%% 6. PLOT LITHOLOGY LINES & PRECISION RULERS
Matrices = {'Quartz', QZ; 'Calcite', LS; 'Dolomite', DL};
LabelRot = [45, 50, 52]; % Estimasi rotasi text untuk Qz, Ls, Dl
for i = 1:size(Matrices,1)
    name = Matrices{i,1}; rho_ma = Matrices{i,2}(1); dt_ma = Matrices{i,2}(2);
    
    % --- A. Draw Main Curves ---
    rho_c = CalcRho(rho_ma, phi_fine);
    dt_w  = CalcDt_Wyllie(dt_ma, phi_fine);
    dt_rhg= CalcDt_RHG(dt_ma, phi_fine);
    
    plot(dt_w, rho_c, 'r-', 'LineWidth', 1.5); % Wyllie (Red)
    plot(dt_rhg, rho_c, 'k-', 'LineWidth', 2.0); % RHG (Black)
    
    % --- B. Draw Precision Rulers (Ticks & Numbers) ---
    rho_t = CalcRho(rho_ma, phi_tick);
    dt_wt = CalcDt_Wyllie(dt_ma, phi_tick);
    dt_rt = CalcDt_RHG(dt_ma, phi_tick);
    
    % Loop untuk setiap 1% porositas
    for j = 1:length(phi_tick)
        p_val = phi_tick(j)*100; % Nilai porositas (0, 1, 2...)
        
        % 1. Gambar Tick Mark kecil ('|') di kedua garis
        plot(dt_wt(j), rho_t(j), 'r|', 'MarkerSize',6, 'LineWidth',1);
        plot(dt_rt(j), rho_t(j), 'k|', 'MarkerSize',6, 'LineWidth',1);
        
        % 2. Tambahkan Nomor hanya pada interval utama (0, 10, 20, 30, 40)
        if mod(p_val, 10) == 0
            % Label pada garis Merah (Wyllie)
            text(dt_wt(j)+0.5, rho_t(j)+0.01, num2str(p_val), ...
                'Color','r', 'FontSize',9, 'FontWeight','bold', ...
                'Rotation', LabelRot(i)-5, 'HorizontalAlignment','left');
            
            % Label pada garis Hitam (RHG) - sedikit offset agar tidak bertumpuk
             text(dt_rt(j)-0.5, rho_t(j)-0.01, num2str(p_val), ...
                'Color','k', 'FontSize',9, 'FontWeight','bold', ...
                'Rotation', LabelRot(i)+5, 'HorizontalAlignment','right');
        end
    end
    % --- C. Matrix Point Label ---
    plot(dt_ma, rho_ma, 'ko', 'MarkerFaceColor','w', 'MarkerSize', 8);
    text(dt_ma-1.5, rho_ma-0.03, name, 'FontWeight','bold', 'FontSize',11);
end

%% 7. FINAL TOUCHES
% Plot Mineral Lain
plot(AN(2), AN(1), 'ko', 'MarkerFaceColor','w', 'MarkerSize',8);
text(AN(2)+1.5, AN(1), 'Anhydrite', 'FontSize',10);
plot(HL(2), HL(1), 'ko', 'MarkerFaceColor','w', 'MarkerSize',8);
text(HL(2)+1.5, HL(1), 'Halite', 'FontSize',10);
% Judul dan Label Sumbu
xlabel('Interval Transit Time, \Delta t (\mus/ft)', 'FontSize',12, 'FontWeight','bold');
ylabel('Bulk Density, \rho_b (g/cc)', 'FontSize',12, 'FontWeight','bold');
title({'Sonic–Density Crossplot';'Precision Ruler & Depth-Colored Data'}, ...
      'FontSize',14, 'FontWeight','bold');
% Legend
h1 = plot(nan,nan, 'k-', 'LineWidth',2);
h2 = plot(nan,nan, 'r-', 'LineWidth',1.5);
legend([h1 h2], {'Empirical (Raymer-Hunt-Gardner)', 'Time Average (Wyllie)'}, ...
       'Location', 'northeast', 'FontSize',10);
hold off;

%% ==========================================================
%  8. AUTOMATIC LITHOLOGY INTERPRETATION (COMMAND WINDOW)
%     Based on Proximity to Raymer-Hunt-Gardner Curves
%     With Thickness Calculation
% ==========================================================
fprintf('\n=================================================================================\n');
fprintf('               LITHOLOGY INTERPRETATION (Sonic-Density RHG)                      \n');
fprintf('=================================================================================\n');
% --- 1. Persiapan Parameter ---
% [rho_ma, dt_ma]
mat_QZ = [2.65, 55.5];
mat_LS = [2.71, 47.6];
mat_DL = [2.87, 43.5];
rho_fl = 1.00; 
% Hitung dz (sampling rate) untuk perhitungan ketebalan yang akurat
dz = median(diff(DEPT));
if isnan(dz), dz = 0.5; end % fallback default
% --- 2. Fungsi Hitung Rho Model (Inverse RHG) ---
% Rumus: rho = (1 - phi)*rho_ma + phi*rho_fl
% Dimana phi dari RHG: phi = 0.625 * (1 - dt_ma/dt_log)
calc_rho_from_dt = @(dt_log, rho_m, dt_m) ...
    (1 - (0.625 * (1 - dt_m./dt_log))) * rho_m + ...
    (0.625 * (1 - dt_m./dt_log)) * rho_fl;
% --- 3. Loop Kalkulasi per Titik Data Reservoir ---
% Pastikan menggunakan data terfilter (DTv, RHOBv, DEPv) dari section sebelumnya
if ~exist('DTv','var') || ~exist('RHOBv','var') || ~exist('DEPv','var')
    warning('Data reservoir terfilter tidak ditemukan. Jalankan section plotting dulu.');
else
    Lithology_ID = zeros(length(DEPv), 1); 
    for i = 1:length(DEPv)
        dt_val = DTv(i);
        rho_val = RHOBv(i);
        
        % Hitung Rho ideal
        rho_model_QZ = calc_rho_from_dt(dt_val, mat_QZ(1), mat_QZ(2));
        rho_model_LS = calc_rho_from_dt(dt_val, mat_LS(1), mat_LS(2));
        rho_model_DL = calc_rho_from_dt(dt_val, mat_DL(1), mat_DL(2));
        
        % Cari selisih terkecil
        diff_QZ = abs(rho_val - rho_model_QZ);
        diff_LS = abs(rho_val - rho_model_LS);
        diff_DL = abs(rho_val - rho_model_DL);
        
        [~, min_idx] = min([diff_QZ, diff_LS, diff_DL]);
        Lithology_ID(i) = min_idx; % 1=Quartz, 2=Calcite, 3=Dolomite
    end
    % --- 4. Grouping per Zona (Continuous Depth) ---
    gap_threshold = 2.0; % ft (Toleransi jeda antar zona)
    jump_idx = find(diff(DEPv) > gap_threshold);
    zone_start_idx = [1; jump_idx + 1];
    zone_end_idx   = [jump_idx; length(DEPv)];
    lit_names = {'Quartz Sandstone', 'Limestone', 'Dolomite'};
    % --- 5. Print Hasil Tabel ---
    % Header Tabel
    fprintf('%-8s | %-23s | %-20s | %-10s\n', ...
            'ZONA', 'DEPTH INTERVAL (ft)', 'DOMINANT LITHOLOGY', 'THICKNESS');
    fprintf('---------------------------------------------------------------------------------\n');
    total_net_res = 0; % Untuk menghitung total
    
    % >>> PERBAIKAN: JANGAN PAKAI 'k' SEBAGAI ITERATOR (Ganti 'iz') <<<
    for iz = 1:length(zone_start_idx)
        s = zone_start_idx(iz);
        e = zone_end_idx(iz);
        
        % Ambil Depth Top & Bottom
        z_top = DEPv(s);
        z_bot = DEPv(e);
        
        % Hitung Ketebalan (Base - Top + 1 step sampling)
        thick = (z_bot - z_top) + dz;
        total_net_res = total_net_res + thick;
        
        % Cari litologi dominan
        zona_liths = Lithology_ID(s:e);
        mode_lith  = mode(zona_liths);
        
        % Print Baris
        fprintf('Zona %-3d : %8.1f - %8.1f    | %-20s | %7.1f ft\n', ...
            iz, z_top, z_bot, lit_names{mode_lith}, thick);
    end
    fprintf('---------------------------------------------------------------------------------\n');
    fprintf('TOTAL NET RESERVOIR THICKNESS DETECTED: %8.1f ft\n', total_net_res);
    fprintf('=================================================================================\n');
end

%% ==========================================================
%  10. MATRIX IDENTIFICATION PLOT (MID PLOT) - MILLIMETER PAPER STYLE
%      (High Density Grid, Exact Reference Match)
% ==========================================================
figure('Color','w','Position',[50 50 1000 950]); 
% --- 1. DATA PREPARATION ---
idx_mid = flag_reservoir & ~isnan(RHOB) & ~isnan(DT) & ~isnan(NPHI) & ~isnan(DEPT);
R_vec  = RHOB(idx_mid);
DT_vec = DT(idx_mid);
N_vec  = NPHI(idx_mid); 
D_vec  = DEPT(idx_mid);
% MID Calculation
rho_f_mid = 1.00; dt_f_mid = 189;
rho_maa = (R_vec - N_vec .* rho_f_mid) ./ (1 - N_vec);
dt_maa  = (DT_vec - N_vec .* dt_f_mid) ./ (1 - N_vec);
% --- 2. MAIN AXES SETUP (GRID RAPAT) ---
axMain = axes('Position', [0.1 0.1 0.85 0.85]); 
hold(axMain, 'on'); box(axMain, 'on');
% Limit Axis (Sesuai Gambar Referensi)
xlim(axMain, [35 75]); 
ylim(axMain, [2.4 3.1]); 
set(axMain, 'YDir', 'reverse');
% --- SETUP GRID "KOTAK KECIL" (MILLIMETER STYLE) ---
% 1. Major Ticks (Garis Tebal Utama)
set(axMain, 'XTick', 35:5:75);       % Tiap 5 unit
set(axMain, 'YTick', 2.4:0.1:3.1);   % Tiap 0.1 unit
% 2. Minor Ticks (Garis Tipis/Kotak Kecil) - INI KUNCINYA
axMain.XAxis.MinorTickValues = 35:1:75;      % Tiap 1 unit (Rapat)
axMain.YAxis.MinorTickValues = 2.4:0.01:3.1; % Tiap 0.01 unit (Sangat Rapat)
% 3. Aktifkan Tampilan Grid
grid(axMain, 'on');          % Nyalakan Major Grid
set(axMain, 'XMinorGrid', 'on', 'YMinorGrid', 'on'); % Nyalakan Minor Grid (Wajib!)
% 4. Styling Grid (Warna Slate Blue ala Chart Fisik)
% Grid Utama (Lebih Gelap)
set(axMain, 'GridColor', [0.2 0.3 0.4], 'GridAlpha', 0.6, 'LineWidth', 1.0);
% Grid Kecil (Lebih Terang tapi Terlihat)
set(axMain, 'MinorGridColor', [0.2 0.3 0.4], 'MinorGridAlpha', 0.3, 'MinorGridLineStyle', '-');
% --- 3. MINERAL DEFINITIONS & PLOT ---
QZ = [55.5, 2.65]; LS = [47.6, 2.71]; DL = [43.5, 2.87]; 
AN = [50.0, 2.98]; Musc = [46, 3.06]; Alb = [49.5, 2.60]; Lang = [52.5, 2.79];
% Draw Triangle
centroid = (QZ + LS + DL) / 3; 
draw_comb_line(axMain, QZ, LS, 4, centroid); 
draw_comb_line(axMain, LS, DL, 4, centroid); 
draw_comb_line(axMain, DL, QZ, 4, centroid); 
% Plot Minerals (Red & Black)
ms_main = 6;
plot(axMain, QZ(1), QZ(2), 'ro', 'MarkerFaceColor','w', 'LineWidth',1.5, 'MarkerSize',ms_main);
plot(axMain, LS(1), LS(2), 'ro', 'MarkerFaceColor','w', 'LineWidth',1.5, 'MarkerSize',ms_main);
plot(axMain, DL(1), DL(2), 'ro', 'MarkerFaceColor','w', 'LineWidth',1.5, 'MarkerSize',ms_main);
plot(axMain, [AN(1) Musc(1) Alb(1) Lang(1)], [AN(2) Musc(2) Alb(2) Lang(2)], ...
    'ko', 'MarkerFaceColor','w', 'LineWidth',1.2, 'MarkerSize',ms_main);
% Labels
text(axMain, QZ(1)+1, QZ(2), 'Quartz', 'Color','r', 'FontWeight','bold', 'FontSize',10);
text(axMain, LS(1)-4, LS(2), 'Calcite', 'Color','r', 'FontWeight','bold', 'FontSize',10);
text(axMain, DL(1)-4, DL(2), 'Dolomite', 'Color','r', 'FontWeight','bold', 'FontSize',10);
text(axMain, AN(1)+1, AN(2), 'Anhydrite', 'Color','k', 'FontWeight','bold', 'FontSize',9);
text(axMain, Musc(1)-4, Musc(2), 'Muscovite', 'Color','k', 'FontSize',8);
text(axMain, Alb(1)-3, Alb(2), 'Albite', 'Color','k', 'FontSize',8);
text(axMain, Lang(1)+1, Lang(2), 'Langbeinite', 'Color','k', 'FontSize',8);
% --- 4. PLOT DATA USER ---
scatter(axMain, dt_maa, rho_maa, 25, D_vec, 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.7);
% Colorbar
colormap(axMain, flipud(jet));
cb = colorbar(axMain, 'Location', 'east');
cb.Position = [0.88 0.1 0.03 0.55]; 
cb.Label.String = 'Depth (ft)';
cb.Label.FontWeight = 'bold';
% Labels
xlabel(axMain, '\Delta t_{maa}, Apparent Matrix Transit Time (\mus/ft)', 'FontWeight','bold', 'FontSize',11);
ylabel(axMain, '\rho_{maa}, Apparent Matrix Density (g/cc)', 'FontWeight','bold', 'FontSize',11);
title(axMain, {'Matrix Identification Plot (MID)';'Apparent Matrix Parameters'}, 'FontWeight','bold', 'FontSize',14);
% --- 5. FIXED ARROWS (CYAN) ---
% Menggunakan plot garis biasa agar nempel di koordinat data
plot(axMain, [35 46], [2.68 2.68], 'c-', 'LineWidth', 2); % Horizontal Line
plot(axMain, 46, 2.68, 'c>', 'MarkerFaceColor','c', 'MarkerSize',8); % Arrowhead
plot(axMain, [52.5 52.5], [3.1 2.8], 'c-', 'LineWidth', 2); % Vertical Line
plot(axMain, 52.5, 2.8, 'c^', 'MarkerFaceColor','c', 'MarkerSize',8); % Arrowhead
% --- 6. INSET BOX (WITH DENSE GRID) ---
axInset = axes('Position', [0.65 0.72 0.23 0.18]); 
hold(axInset, 'on'); 
% Background Putih Solid (Menutupi grid belakang)
set(axInset, 'Color', 'w', 'Box', 'on', 'YDir', 'reverse');
xlim(axInset, [60 75]); ylim(axInset, [2.0 2.15]);
% Grid Inset (Sama Rapatnya)
set(axInset, 'XTick', 60:5:75, 'YTick', 2.0:0.05:2.15);
axInset.XAxis.MinorTickValues = 60:1:75;       % Rapat 1 unit
axInset.YAxis.MinorTickValues = 2.0:0.01:2.15; % Rapat 0.01 unit
grid(axInset, 'on'); 
set(axInset, 'XMinorGrid', 'on', 'YMinorGrid', 'on'); % Aktifkan
set(axInset, 'GridColor',[0.2 0.3 0.4], 'GridAlpha',0.6);
set(axInset, 'MinorGridColor',[0.2 0.3 0.4], 'MinorGridAlpha',0.3, 'MinorGridLineStyle','-');
% Plot Inset Data
Halite = [67.0, 2.03]; Ortho = [69, 2.13]; 
plot(axInset, Halite(1), Halite(2), 'ko', 'MarkerFaceColor','w', 'LineWidth',1.2, 'MarkerSize',5);
plot(axInset, Ortho(1), Ortho(2), 'ko', 'MarkerFaceColor','w', 'LineWidth',1.2, 'MarkerSize',5);
text(axInset, Halite(1)-4, Halite(2), 'Halite', 'Color','k', 'FontSize',8, 'FontWeight','bold');
text(axInset, Ortho(1)-4, Ortho(2)+0.015, 'Orthoclase', 'Color','k', 'FontSize',8);
hold(axMain, 'off'); hold(axInset, 'off');

% --- HELPER FUNCTION ---
function draw_comb_line(ax, p1, p2, n_ticks, center_pt)
    plot(ax, [p1(1) p2(1)], [p1(2) p2(2)], 'k-', 'LineWidth', 1.5);
    v = p2 - p1;
    t_vals = linspace(0, 1, n_ticks+2); t_vals = t_vals(2:end-1);
    
    % Normal vector calculation
    n_cand1 = [-v(2), v(1)]; n_cand2 = [v(2), -v(1)]; 
    midpoint = (p1 + p2) / 2;
    vec_to_center = center_pt - midpoint;
    
    aspect_ratio = 60; 
    dot1 = n_cand1(1)*vec_to_center(1) + (n_cand1(2)*aspect_ratio)*vec_to_center(2);
    if dot1 > 0, n_use = n_cand1; else, n_use = n_cand2; end
    
    len_vis = sqrt(n_use(1)^2 + (n_use(2)*aspect_ratio)^2);
    n_use = n_use / len_vis;
    
    tick_len_x = 0.8; 
    for t = t_vals
        pt_start = p1 + t*v;
        dx = n_use(1) * tick_len_x;
        dy = n_use(2) * (tick_len_x / aspect_ratio) * 60;
        plot(ax, [pt_start(1) pt_start(1)+dx], [pt_start(2) pt_start(2)+dy], 'k-', 'LineWidth', 1);
    end
end

%% ==========================================================
%  14. PRINT LITHOLOGY TABLE (MID PLOT) - COMMAND WINDOW
%      (Table Format Matches Your Request)
% ==========================================================
% --- 1. DATA PREPARATION & CALCULATION ---
% Gunakan filter yang sama dengan Section 10 (MID Plot)
idx_print = flag_reservoir & ~isnan(RHOB) & ~isnan(DT) & ~isnan(NPHI) & ~isnan(DEPT);
if sum(idx_print) == 0
    % Jika tidak ada data, lewati
    fprintf('\nTidak ada data reservoir untuk interpretasi MID.\n');
else
    % Ambil Data
    Rv = RHOB(idx_print);
    Dv = DT(idx_print);
    Nv = NPHI(idx_print);
    D_print = DEPT(idx_print);
    
    % Parameter Fluida (Fresh Water)
    rho_f_mid = 1.00; 
    dt_f_mid  = 189;
    
    % Hitung Apparent Matrix (Rumus MID Plot)
    rho_maa_v = (Rv - Nv .* rho_f_mid) ./ (1 - Nv);
    dt_maa_v  = (Dv - Nv .* dt_f_mid) ./ (1 - Nv);
    
    % --- 2. REFERENSI MINERAL [dt_ma, rho_ma] ---
    Ref_QZ = [55.5, 2.65]; % Quartz
    Ref_LS = [47.6, 2.71]; % Calcite
    Ref_DL = [43.5, 2.87]; % Dolomite
    Ref_AN = [50.0, 2.98]; % Anhydrite
    Ref_HA = [67.0, 2.03]; % Halite
    
    % --- 3. IDENTIFIKASI LITOLOGI (WEIGHTED DISTANCE) ---
    % Faktor pembobot karena skala DT (0-100) != skala Rho (2-3)
    W_rho = 50; 
    
    Lith_ID = zeros(length(D_print), 1);
    
    for i = 1:length(D_print)
        dt_i  = dt_maa_v(i);
        rho_i = rho_maa_v(i);
        
        % Jarak ke Mineral
        d_QZ = sqrt((dt_i - Ref_QZ(1))^2 + (W_rho * (rho_i - Ref_QZ(2)))^2);
        d_LS = sqrt((dt_i - Ref_LS(1))^2 + (W_rho * (rho_i - Ref_LS(2)))^2);
        d_DL = sqrt((dt_i - Ref_DL(1))^2 + (W_rho * (rho_i - Ref_DL(2)))^2);
        d_AN = sqrt((dt_i - Ref_AN(1))^2 + (W_rho * (rho_i - Ref_AN(2)))^2);
        d_HA = sqrt((dt_i - Ref_HA(1))^2 + (W_rho * (rho_i - Ref_HA(2)))^2);
        
        % Cari Minimum
        [~, id_min] = min([d_QZ, d_LS, d_DL, d_AN, d_HA]);
        Lith_ID(i) = id_min; 
    end
    
    % Mapping Nama (Sesuaikan agar pas di kolom)
    LitNames = {'Quartz Sandstone', 'Limestone', 'Dolomite', 'Anhydrite', 'Halite'};
    
    % --- 4. GROUPING ZONES & PRINTING ---
    dz = median(diff(DEPT)); 
    if isnan(dz), dz = 0.5; end
    
    % Deteksi Zona (Jika ada gap > 2 ft atau litologi berubah drastis - opsional)
    % Di sini kita grouping berdasarkan depth continuity saja sesuai contoh Anda
    jump_loc = find(diff(D_print) > 2.0);
    z_starts = [1; jump_loc + 1];
    z_ends   = [jump_loc; length(D_print)];
    % --- CETAK HEADER (FORMAT PERSIS) ---
    fprintf('\n');
    fprintf('=================================================================================\n');
    fprintf('               LITHOLOGY INTERPRETATION (MID Plot)                               \n');
    fprintf('=================================================================================\n');
    fprintf('%-8s | %-23s | %-20s | %-10s \n', 'ZONA', 'DEPTH INTERVAL (ft)', 'DOMINANT LITHOLOGY', 'THICKNESS');
    fprintf('---------------------------------------------------------------------------------\n');
    Total_Thick = 0;
    
    % >>> PERBAIKAN: JANGAN PAKAI 'k' SEBAGAI ITERATOR (Ganti 'iz') <<<
    for iz = 1:length(z_starts)
        s = z_starts(iz);
        e = z_ends(iz);
        
        % Data Zona
        top_d  = D_print(s);
        base_d = D_print(e);
        thick  = (base_d - top_d) + dz;
        Total_Thick = Total_Thick + thick;
        
        % Litologi Dominan (Modus)
        mode_lit = mode(Lith_ID(s:e));
        lit_str  = LitNames{mode_lit};
        
        % Print Baris dengan Format Rapi
        % %-8s : Rata kiri
        % %8.1f: Rata kanan angka float
        fprintf('Zona %-3d : %8.1f - %8.1f    | %-20s | %7.1f ft\n', ...
            iz, top_d, base_d, lit_str, thick);
    end
    
    fprintf('---------------------------------------------------------------------------------\n');
    fprintf('TOTAL RESERVOIR THICKNESS : %36.1f ft\n', Total_Thick);
    fprintf('=================================================================================\n');
end

%% ==========================================================
%  NET RESERVOIR & NET PAY FLAG
% ==========================================================
Vsh_cut = 0.50;
phi_cut = 0.15;
Sw_cut  = 0.70;
k_cut   = 1;     % mD (opsional)
flag_net_res = (Vsh <= Vsh_cut) & (phi_eff >= phi_cut);
flag_net_pay = flag_net_res & (Sw <= Sw_cut) & (k >= k_cut); % k sudah aman

%% ---- 3.5 Lithology Classification 
LithCode = NaN(size(DEPT));
% 1 = Shale, 2 = Sand, 3 = Limestone, 4 = Siltstone
for i = 1:length(DEPT)
    if isnan(GR(i)) || isnan(RHOB(i)) || isnan(NPHI(i))
        continue
    end
    % SHALE
    if GR(i) > 75
        LithCode(i) = 1;
    % LIMESTONE
    elseif GR(i) < 60 && RHOB(i) >= 2.45 && NPHI(i) <= 0.15
        LithCode(i) = 3;
    % SANDSTONE
    elseif GR(i) < 60 && RHOB(i) < 2.45
        LithCode(i) = 2;
    % SILTSTONE 
    else
        LithCode(i) = 4;
    end
end

%% ==========================================================
% 4. VISUALIZATION – QUICKLOOK & PETROPHYSICS (FIXED)
% ==========================================================
figure('Color','w','Position',[50 50 2000 900]);
sgtitle('TUGAS AKHIR GABE','Color','k','FontSize',18);

%% ---- TRACK 1: LITHOLOGY ----
ax0 = subplot(1,9,1);
hold on
for i = 1:length(DEPT)-1
    if isnan(LithCode(i)), continue; end
    z1 = DEPT(i);
    z2 = DEPT(i+1);
    switch LithCode(i)
        case 1, col = [0 0.7 0];        % Shale
        case 2, col = [1 1 0];          % Sand
        case 3, col = [0 0.4 1];        % Limestone
        case 4, col = [0.6 0.4 0.2];    % Siltstone
    end
    patch([0 1 1 0],[z1 z1 z2 z2],col,'EdgeColor','none');
end
set(gca,'YDir','reverse','Color','w','XColor','k','YColor','k')
xlim([0 1])
xticks([])
ylabel('Depth (ft)','Color','k')
title('Lithology','Color','k')
grid on

%% ---- TRACK 2: GR ----
ax1 = subplot(1,9,2);
hold on
% Plot GR dulu
plot(GR, DEPT,'g','LineWidth',1.3);
% Shading mengikuti kurva GR
for i = 1:length(DEPT)-1
    if isnan(GR(i)) || isnan(GR(i+1))
        continue
    end
    z1 = DEPT(i);
    z2 = DEPT(i+1);
    if GR(i) <= GR_cutoff
        col = [1 1 0];      % KUNING → permeabel
    else
        col = [0 0.6 0];    % HIJAU → non-permeabel
    end
    patch([GR(i) GR_cutoff GR_cutoff GR(i)], ...
          [z1 z1 z2 z2], ...
          col, ...
          'FaceAlpha',0.35, ...
          'EdgeColor','none');
end
% Cutoff line
xline(GR_cutoff,'r--','GR cutoff','LineWidth',1.2);
set(gca,'YDir','reverse','Color','w','XColor','k','YColor','k')
xlabel('GR (API)','Color','k')
title('Gamma Ray','Color','k')
xlim([0 150])
grid on

%% ---- TRACK 3: RESISTIVITY ----
ax2 = subplot(1,9,3);
semilogx(LLD,DEPT,'r','LineWidth',1); hold on
semilogx(LLS,DEPT,'m','LineWidth',1);
set(gca,'YDir','reverse','Color','w','XColor','k','YColor','k')
xlabel('Rt (ohm.m)','Color','k')
title('Resistivity','Color','k')
xlim([0.2 2000]); grid on
legend({'LLD','LLS'},'TextColor','k','Location','best')

%% ---- TRACK 4: DENSITY–NEUTRON (IP-LIKE HC ZONE) ----
ax3 = subplot(1,9,4);
hold on
% --- Plot RHOB & NPHI (tetap sama) ---
plot(RHOB,DEPT,'r','LineWidth',1);
xlim([1.7 2.7])
% NPHI scaled to RHOB axis (jangan ubah formula ini)
NPHI_scaled = 1.7 + (0.6 - NPHI).*(1.0/0.6);
plot(NPHI_scaled,DEPT,'c','LineWidth',1);
% --------------------
% IP-like adaptive filtering (robust & defensible)
% --------------------
% 1) Zona permeabel (gate pertama) - gunakan Vsh atau GR sesuai data
flag_perm = false(size(DEPT));
if exist('Vsh','var')
    flag_perm = Vsh < 0.35;        % default IP-like cutoff, bisa disesuaikan
else
    % fallback: jika tidak ada Vsh, gunakan GR jika tersedia
    if exist('GR','var')
        flag_perm = GR < (nanmedian(GR));  % sederhana: lebih rendah dari median GR
    else
        flag_perm = true(size(DEPT)); % jika tidak ada informasi, jangan blokir semua
    end
end
% 2) Crossover magnitude
dphi = RHOB - NPHI_scaled;
% Build mask of indices with usable data and in permeable zones
valid_mask = flag_perm & ~isnan(dphi) & ~isnan(LLD) & ~isnan(RHOB);
% If not enough valid samples, set sensible defaults
if sum(valid_mask) >= 5
    dphi_cut_low  = prctile(dphi(valid_mask),60);  % threshold untuk oil-level crossover
    dphi_cut_high = prctile(dphi(valid_mask),85);  % threshold untuk gas-strong crossover
    rt_cut         = prctile(LLD(valid_mask),40);  % resistivity relatif rendah
    rhob_perm_med  = prctile(RHOB(flag_perm & ~isnan(RHOB)),40); % densitas relatif rendah di zone permeabel
else
    % fallback fixed but reasonable defaults if dataset kecil
    dphi_cut_low  = 0.05;
    dphi_cut_high = 0.15;
    rt_cut        = 20;        % ohm.m (empiris)
    rhob_perm_med = 2.25;
end
% 3) Flags final (mutually exclusive)
flag_crossover = dphi > dphi_cut_low;
flag_crossover_strong = dphi > dphi_cut_high;
flag_rt_low = LLD < rt_cut;
% Gas: permeable + crossover kuat + resistivitas relatif rendah + RHOB relatif rendah
flag_gas = flag_perm & flag_crossover_strong & flag_rt_low & (RHOB < rhob_perm_med);
% Oil: permeable + crossover (at least low) + not gas + resistivity not extremely high
flag_oil = flag_perm & flag_crossover & ~flag_gas & (LLD < prctile(LLD(valid_mask),80));
% Ensure mutual exclusivity (safety)
flag_oil(flag_gas) = false;
% --------------------
% Draw shading: oil first, gas last (gas on top so it is not visually hidden)
% --------------------
for i = 1:length(DEPT)-1
    if isnan(RHOB(i)) || isnan(NPHI_scaled(i)) || isnan(LLD(i))
        continue
    end
    z1 = DEPT(i);
    z2 = DEPT(i+1);
    % OIL (KUNING) — digambar dulu
    if flag_oil(i)
        patch([NPHI_scaled(i) RHOB(i) RHOB(i) NPHI_scaled(i)], ...
              [z1 z1 z2 z2], ...
              [1 1 0], 'FaceAlpha',0.30, 'EdgeColor','none');
    end
    % GAS (MERAH) — digambar terakhir (prioritas visual)
    if flag_gas(i)
        patch([NPHI_scaled(i) RHOB(i) RHOB(i) NPHI_scaled(i)], ...
              [z1 z1 z2 z2], ...
              [1 0 0], 'FaceAlpha',0.40, 'EdgeColor','none');
    end
end
% --------------------
% Axis style & legend (manual legend ensures correct colors)
% --------------------
set(gca,'YDir','reverse','Color','w','XColor','k','YColor','k','FontSize',10)
xlabel('RHOB / NPHI','Color','k')
title('Density–Neutron','Color','k')
xlim([1.7 2.7])
grid on
% Create dummy handles for legend so colors map correctly
h_rhob = plot(nan,nan,'r','LineWidth',1);
h_nphi = plot(nan,nan,'c','LineWidth',1);
h_gas  = patch(nan,nan,[1 0 0],'FaceAlpha',0.40,'EdgeColor','none');
h_oil  = patch(nan,nan,[1 1 0],'FaceAlpha',0.30,'EdgeColor','none');
legend([h_rhob h_nphi h_gas h_oil], {'RHOB','NPHI','Gas','Oil'}, ...
       'TextColor','k','Location','best');

%% ---- TRACK 5: Vshale (Reservoir Only) ----
ax4 = subplot(1,9,5);
hold on
plot(Vsh, DEPT, 'Color', [0.4 0.3 0.1], 'LineWidth', 1.2); % Warna Coklat Gelap
% Shading (Arsiran) agar terlihat volume-nya
% Area diarsir dari Kurva ke Kanan (Shale volume) atau Kiri (Sand volume)
% Disini kita arsir area shale (Kanan kurva) untuk visualisasi kotoran
for i = 1:length(DEPT)-1
    if isnan(Vsh(i)) || isnan(Vsh(i+1))
        continue; 
    end
    % Gambar area abu-abu/kuning
    patch([Vsh(i) 1 1 Vsh(i)], ...
          [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], ...
          [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % Abu-abu (Shale)
end
% Garis Cutoff (Merah Putus-putus)
xline(0.35, 'r--', 'LineWidth', 1.5);
% Kosmetik Axis
set(gca, 'YDir','reverse', 'Color','w', 'XColor','k', 'YColor','k');
xlabel('Vsh (v/v)', 'Color','k', 'FontSize',8, 'FontWeight','bold');
title('Vshale', 'Color','k', 'FontSize',10, 'FontWeight','bold');
xlim([0 1]);
grid on;
set(gca, 'XTick', 0:0.2:1); % Grid per 0.2

%% ---- TRACK 6: EFFECTIVE POROSITY ----
ax5 = subplot(1,9,6);
plot(phi_eff,DEPT,'c','LineWidth',1);
set(gca,'YDir','reverse','Color','w','XColor','k','YColor','k')
xlabel('\phi_{eff}','Color','k')
title('Effective Porosity','Color','k')
xlim([0 0.4])
grid on

%% ---- TRACK 7: WATER SATURATION ----
ax6 = subplot(1,9,7);
plot(Sw,DEPT,'y','LineWidth',1);
set(gca,'YDir','reverse','Color','w','XColor','k','YColor','k')
xlabel('Sw','Color','k')
title('Water Saturation','Color','k')
xlim([0 1])
grid on

%% ---- TRACK 8: PERMEABILITY (FINAL – SUDAH AMAN) ----
ax7 = subplot(1,9,8);
hold on
% Karena variabel 'k' sudah tidak dirusak oleh loop, kode ini aman
idx_plot = ~isnan(k) & ~isnan(DEPT);
semilogx(k(idx_plot), DEPT(idx_plot), 'b-', 'LineWidth', 1.2);
if exist('flag_reservoir','var')
    idx_res_plot = idx_plot & flag_reservoir;
    semilogx(k(idx_res_plot), DEPT(idx_res_plot), 'r.', 'MarkerSize', 4);
end
set(gca,'YDir','reverse','Color','w','XColor','k','YColor','k')
set(gca,'XScale','log')
xlim([0.01 10000])
set(gca,'XTick',[0.01 0.1 1 10 100 1000 10000])
set(gca,'XTickLabel',{'0.01','0.1','1','10','100','1K','10K'})
xlabel('Permeability (mD)')
title('Permeability')
grid on
grid minor
hold off

%% ---- LINK DEPTH AXIS ----
linkaxes([ax0 ax1 ax2 ax3 ax4 ax5 ax6 ax7],'y');

%% ========================================================================
%               SUMMARY TABLE (AVERAGE PROPERTIES PER ZONE)
%      (Letakkan bagian ini di paling akhir script Anda)
% =========================================================================
% 1. SETUP HEADER TABEL
fprintf('\n\n');
fprintf('=================================================================================\n');
fprintf('               PETROPHYSICS CALCULATION SUMMARY (AVERAGE VALUES)                 \n');
fprintf('=================================================================================\n');
fprintf('%-8s | %-23s | %-5s | %-9s | %-5s | %-8s\n', ...
        'ZONA', 'DEPTH INTERVAL (ft)', 'VSH', 'POROSITAS', 'SW', 'PERM(mD)');
fprintf('---------------------------------------------------------------------------------\n');
% 2. DETEKSI INTERVAL ZONA
% Kita gunakan flag_reservoir. Jika flag_net_pay ada, Anda bisa menggantinya ke situ.
% Pastikan flag bersih dari NaN agar tidak error
if exist('flag_reservoir','var')
    flag_summary = flag_reservoir;
else
    % Fallback jika flag belum didefinisikan
    flag_summary = Vsh < 0.5 & phi_eff > 0.05; 
end
% Filter tambahan: pastikan data perhitungan ada (tidak NaN)
flag_summary = flag_summary & ~isnan(DEPT) & ~isnan(Vsh) & ~isnan(phi_eff) & ~isnan(Sw);
% Logika mencari start & end index dari zona yang kontinyu
padded_flag = [0; flag_summary; 0]; % Tambah 0 di ujung agar deteksi tepi akurat
diff_flag   = diff(padded_flag);
start_idx   = find(diff_flag == 1);      % Indeks awal zona
end_idx     = find(diff_flag == -1) - 1; % Indeks akhir zona
% 3. LOOP KALKULASI & PRINT
if isempty(start_idx)
    fprintf('TIDAK ADA ZONA POTENSI YANG TERDETEKSI.\n');
else
    % >>> PERBAIKAN: JANGAN PAKAI 'k' SEBAGAI ITERATOR (Ganti 'iz') <<<
    for iz = 1:length(start_idx)
        % Ambil rentang indeks baris untuk zona ke-i
        idx_range = start_idx(iz):end_idx(iz);
        
        % Ambil Depth Top & Bottom
        z_top = DEPT(start_idx(iz));
        z_bot = DEPT(end_idx(iz));
        
        % Hitung Rata-rata (Mean) parameter di zona tersebut
        % Menggunakan 'mean' biasa. Jika MATLAB versi lama error, pakai nanmean
        avg_vsh  = mean(Vsh(idx_range));
        avg_phi  = mean(phi_eff(idx_range));
        avg_sw   = mean(Sw(idx_range));
        
        % Khusus Permeability: gunakan geometric mean atau arithmetic mean (di sini arithmetic)
        % Handle jika k belum dihitung atau NaN
        if exist('k','var')
            k_zona = k(idx_range);
            avg_k = mean(k_zona(~isnan(k_zona))); 
            if isempty(avg_k) || isnan(avg_k), avg_k = 0; end
        else
            avg_k = 0;
        end
        
        % PRINT BARIS (Format presisi angka desimal disesuaikan agar rapi)
        fprintf('Zona %-3d : %8.1f - %8.1f           | %.3f |   %.3f   | %.3f | %8.2f\n', ...
                iz, z_top, z_bot, avg_vsh, avg_phi, avg_sw, avg_k);
    end
end
fprintf('---------------------------------------------------------------------------------\n');

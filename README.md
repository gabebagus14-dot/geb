%% ========================================================================
%                          PETROPHYSICAL ANALYSIS 
%                    TUGAS AKHIR GABE BAGUS ARTHANTHA
% =========================================================================
clc; clear; close all;

%% ========================================================================
%   1. INPUT DATA LOG & QUALITY CONTROL (QC)
% =========================================================================
[file, path] = uigetfile('*.las', 'Pilih File Data LAS');
if isequal(file,0), error('ANALISIS BERHENTI: Tidak ada file yang dipilih.'); end
filename = fullfile(path, file);
fprintf('\n=>  Membaca file: %s\n', file);

fid = fopen(filename,'r'); curve_names = {}; in_curve = false;
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if startsWith(upper(line), '~C'), in_curve = true; continue;
    elseif startsWith(upper(line), '~') && in_curve, in_curve = false; end
    if startsWith(upper(line), '~A'), break; end
    if in_curve && ~isempty(line) && ~startsWith(line, '#')
        tokens = regexp(line, '^([^.\s]+)', 'tokens');
        if ~isempty(tokens), curve_names{end+1} = upper(tokens{1}{1}); end
    end
end

ascii_data = [];
while true
    line = fgetl(fid); if ~ischar(line), break; end
    nums = str2num(line); if ~isempty(nums), ascii_data = [ascii_data; nums]; end
end
fclose(fid);
NULL = -999.25; ascii_data(ascii_data <= NULL) = NaN; 

idx_dept = []; for a={'DEPT', 'DEPTH', 'DEP'}, idx=find(strcmpi(curve_names, a{1}),1); if ~isempty(idx), idx_dept=idx; break; end; end
if isempty(idx_dept), error('ANALISIS BERHENTI: Kolom Depth tidak ditemukan!'); end
DEPT = ascii_data(:, idx_dept);

idx_gr   = []; for a={'GR', 'GAMMA', 'CGR', 'SGR'}, idx=find(strcmpi(curve_names, a{1}),1); if ~isempty(idx), idx_gr=idx; break; end; end
idx_lld  = []; for a={'LLD', 'RT', 'ILD', 'RD', 'AHT90', 'RDEP'}, idx=find(strcmpi(curve_names, a{1}),1); if ~isempty(idx), idx_lld=idx; break; end; end
idx_lls  = []; for a={'LLS', 'RXO', 'ILM', 'RS', 'AHT10', 'RMIC'}, idx=find(strcmpi(curve_names, a{1}),1); if ~isempty(idx), idx_lls=idx; break; end; end
idx_rhob = []; for a={'RHOB', 'ZDEN', 'DEN', 'RHOZ', 'DENB'}, idx=find(strcmpi(curve_names, a{1}),1); if ~isempty(idx), idx_rhob=idx; break; end; end
idx_nphi = []; for a={'NPHI', 'CNC', 'PHIN', 'NPOR', 'HNPO'}, idx=find(strcmpi(curve_names, a{1}),1); if ~isempty(idx), idx_nphi=idx; break; end; end
idx_dt   = []; for a={'DT', 'AC', 'DTC', 'DTCO', 'ACCO'}, idx=find(strcmpi(curve_names, a{1}),1); if ~isempty(idx), idx_dt=idx; break; end; end
idx_cali = []; for a={'CALI', 'CAL', 'HCAL', 'CALS'}, idx=find(strcmpi(curve_names, a{1}),1); if ~isempty(idx), idx_cali=idx; break; end; end

n_rows = length(DEPT);
if ~isempty(idx_gr),   GR = ascii_data(:, idx_gr);     else, GR = NaN(n_rows, 1);   end
if ~isempty(idx_lld),  LLD = ascii_data(:, idx_lld);   else, LLD = NaN(n_rows, 1);  end
if ~isempty(idx_lls),  LLS = ascii_data(:, idx_lls);   else, LLS = NaN(n_rows, 1);  end
if ~isempty(idx_rhob), RHOB = ascii_data(:, idx_rhob); else, RHOB = NaN(n_rows, 1); end
if ~isempty(idx_nphi), NPHI = ascii_data(:, idx_nphi); else, NPHI = NaN(n_rows, 1); end
if ~isempty(idx_dt),   DT = ascii_data(:, idx_dt);     else, DT = NaN(n_rows, 1);   end
if ~isempty(idx_cali), CALI = ascii_data(:, idx_cali); else, CALI = NaN(n_rows, 1); end

dz = median(diff(DEPT), 'omitnan'); 
if isnan(dz) || dz == 0, error('ANALISIS BERHENTI: Resolusi kedalaman (dz) tidak valid.'); end

%% ========================================================================
%   2. FULL WELL TRIPLE COMBO LOG
% =========================================================================

margin_left = 0.08; margin_right = 0.05; margin_bottom = 0.05; 
h_header = 0.10; h_track = 0.82; gap = 0.005;
y_header = margin_bottom + h_track; y_track = margin_bottom;
font_name = 'Arial'; grid_color = [0.75 0.75 0.75]; axis_size = 9;
total_width_tc = 1.0 - margin_left - margin_right; w_tc = (total_width_tc - (2 * gap)) / 3; 
lbl_tc = {{'GR (Green) / CALI (Black)', '(API) / (INCH)', '0 / 6', '150 / 16'}, {'RESISTIVITY', '(OHM.M)', '0.2', '2000'}, {'RHOB (Red) / NPHI (Blue)', '(G/CC) / (V/V)', '1.7 / 0.6', '2.7 / 0'}};

GR_valid_fw = GR(~isnan(GR) & GR > 0);
if ~isempty(GR_valid_fw)
    GR_clean_fw = prctile(GR_valid_fw, 5); GR_shale_fw = prctile(GR_valid_fw, 95);
    GR_cutoff_fw = GR_clean_fw + 0.5 * (GR_shale_fw - GR_clean_fw); 
    fprintf('=> Nilai GR Cut-off (Shale Baseline 50%%) Full Well : %.2f API\n', GR_cutoff_fw);
end

fig_full = figure('Name','Full Well Triple Combo Log','Color','w','Position',[10 50 1000 800], 'NumberTitle','off');
ax_full = gobjects(1, 3); curr_x = margin_left;
for i = 1:3
    axes('Position',[curr_x, y_header, w_tc, h_header]); axis off; rectangle('Position',[0 0 1 1],'LineWidth',1); 
    text(0.5, 0.75, lbl_tc{i}{1}, 'Horiz','center','FontWeight','bold','FontSize',10); text(0.5, 0.45, lbl_tc{i}{2}, 'Horiz','center','FontSize',9,'FontAngle','italic');
    text(0.02, 0.1, lbl_tc{i}{3}, 'Horiz','left','FontSize',9,'FontWeight','bold'); text(0.98, 0.1, lbl_tc{i}{4}, 'Horiz','right','FontSize',9,'FontWeight','bold');
    
    ax_full(i) = axes('Position', [curr_x, y_track, w_tc, h_track]);
    set(ax_full(i), 'FontName', font_name, 'FontSize', axis_size, 'Box', 'on', 'LineWidth', 1.0, 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k', 'YDir', 'reverse', 'XAxisLocation', 'bottom', 'XTickLabel', [], 'XGrid', 'on', 'YGrid', 'on', 'GridColor', grid_color, 'GridAlpha', 0.8, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'MinorGridColor', grid_color, 'MinorGridAlpha', 0.4);
    hold(ax_full(i), 'on'); curr_x = curr_x + w_tc + gap;
end

% Track 1: GR & CALI (FULL LOG)
axes(ax_full(1));
if ~isempty(GR_valid_fw)
    for i = 1:length(DEPT)-1
        if isnan(GR(i)) || isnan(GR(i+1)); continue; end
        if GR(i) <= GR_cutoff_fw, col = [1 1 0]; else, col = [0 0.6 0]; end
        patch([GR(i) GR_cutoff_fw GR_cutoff_fw GR(i)], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], col, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
    end
end
plot(GR, DEPT, 'g', 'LineWidth', 1.2); if ~isempty(GR_valid_fw), xline(GR_cutoff_fw, 'r--', 'LineWidth', 1.5); end
% --- PLOT CALIPER ---
CALI_scaled = (CALI - 4) .* (150 / (14 - 4));
BitSize_scaled = (6.0 - 4) * (150 / (14 - 4)); % Bit Size aktual = 6.0 inch
plot(CALI_scaled, DEPT, 'k-', 'LineWidth', 1.5); % Garis Caliper Hitam
xline(BitSize_scaled, 'k--', 'LineWidth', 1.5); % Garis Putus-putus Bit Size
% ----------------------------------------------------------
xlim([0 150]); ylabel('Depth (ft)', 'FontWeight', 'bold', 'FontSize', 11);
set(ax_full(1), 'Layer', 'top');
% Track 2: Resistivity (FULL LOG)
axes(ax_full(2)); LLD_plot = LLD; LLS_plot = LLS; LLD_plot(LLD_plot<=0)=NaN; LLS_plot(LLS_plot<=0)=NaN; 
semilogx(LLD_plot, DEPT, 'r', 'LineWidth', 1.2); semilogx(LLS_plot, DEPT, 'm', 'LineWidth', 0.8); 
set(gca,'XScale','log'); xlim([0.2 2000]); set(gca,'YTickLabel',[]);
% Track 3: Density-Neutron (FULL LOG)
axes(ax_full(3)); NPHI_scaled = 2.7 - (NPHI .* (1.0/0.6)); separation = NPHI_scaled - RHOB;
for i = 1:length(DEPT)-1
    if isnan(RHOB(i)) || isnan(NPHI_scaled(i)); continue; end
    if separation(i) > 0, patch([RHOB(i) NPHI_scaled(i) NPHI_scaled(i) RHOB(i)], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], [1 1 0], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    else, patch([RHOB(i) NPHI_scaled(i) NPHI_scaled(i) RHOB(i)], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], [0 0.6 0], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); end
end
plot(RHOB, DEPT, 'r', 'LineWidth', 1.2); plot(NPHI_scaled, DEPT, 'b--', 'LineWidth', 1.2); 
xlim([1.7 2.7]); set(gca,'YTickLabel',[]); set(ax_full(3), 'Layer', 'top'); 
linkaxes(ax_full,'y'); set(ax_full(1), 'YLim', [min(DEPT) max(DEPT)]); 
exportgraphics(fig_full, '0_Full_Well_Triple_Combo_600DPI.png', 'Resolution', 600);

%% ========================================================================
%   3. INPUT ZONA TARGET & VISUALISASI TRIPLE COMBO ZONA PROSPEK
% =========================================================================
pause(2); 

prompt_zona = {
    'Top Formasi / Zona Target (ft):', ...
    'Base Formasi / Zona Target (ft):', ...
    'Top Zona Data Rusak/Washout (ft) [Isi 0 jika aman]:', ...
    'Base Zona Data Rusak/Washout (ft) [Isi 0 jika aman]:'
};
dlgtitle_zona = 'Identifikasi Zona Prospek & QC Washout'; dims = [1 65];
definput_zona = {'5500', '7400', '6600', '6800'};
ans_zona = inputdlg(prompt_zona, dlgtitle_zona, dims, definput_zona);
if isempty(ans_zona), error('ANALISIS BERHENTI: Input dibatalkan.'); end

Top_Formasi  = str2double(ans_zona{1}); Base_Formasi = str2double(ans_zona{2});
Washout_Top  = str2double(ans_zona{3}); Washout_Bot  = str2double(ans_zona{4});

flag_reservoir = (DEPT >= Top_Formasi) & (DEPT <= Base_Formasi);
if Washout_Top > 0 && Washout_Bot > 0
    valid_log = ~(DEPT >= Washout_Top & DEPT <= Washout_Bot);
else
    valid_log = true(size(DEPT));
end
zone = flag_reservoir & valid_log & ~isnan(GR) & ~isnan(LLD) & ~isnan(RHOB) & ~isnan(NPHI);

fig_tc = figure('Name','Triple Combo Log Target','Color','w','Position',[50 50 1000 800], 'NumberTitle','off');
ax_tc = gobjects(1, 3); curr_x = margin_left;
for i = 1:3
    axes('Position',[curr_x, y_header, w_tc, h_header]); axis off; rectangle('Position',[0 0 1 1],'LineWidth',1); 
    text(0.5, 0.75, lbl_tc{i}{1}, 'Horiz','center','FontWeight','bold','FontSize',10); text(0.5, 0.45, lbl_tc{i}{2}, 'Horiz','center','FontSize',9,'FontAngle','italic');
    text(0.02, 0.1, lbl_tc{i}{3}, 'Horiz','left','FontSize',9,'FontWeight','bold'); text(0.98, 0.1, lbl_tc{i}{4}, 'Horiz','right','FontSize',9,'FontWeight','bold');
    
    ax_tc(i) = axes('Position', [curr_x, y_track, w_tc, h_track]);
    set(ax_tc(i), 'FontName', font_name, 'FontSize', axis_size, 'Box', 'on', 'LineWidth', 1.0, 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k', 'YDir', 'reverse', 'XAxisLocation', 'bottom', 'XTickLabel', [], 'XGrid', 'on', 'YGrid', 'on', 'GridColor', grid_color, 'GridAlpha', 0.8, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'MinorGridColor', grid_color, 'MinorGridAlpha', 0.4);
    hold(ax_tc(i), 'on'); curr_x = curr_x + w_tc + gap;
end

% Track 1: GR & CALI (ZONA PROSPEK)
axes(ax_tc(1));
if ~isempty(GR_valid_fw)
    for i = 1:length(DEPT)-1
        if isnan(GR(i)) || isnan(GR(i+1)); continue; end
        if GR(i) <= GR_cutoff_fw, col = [1 1 0]; else, col = [0 0.6 0]; end
        patch([GR(i) GR_cutoff_fw GR_cutoff_fw GR(i)], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], col, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
    end
end
plot(GR, DEPT, 'g', 'LineWidth', 1.5); if ~isempty(GR_valid_fw), xline(GR_cutoff_fw, 'r--', 'LineWidth', 1.3); end

% --- PLOT CALIPER & SHADING WASHOUT ---
CALI_scaled = (CALI - 4) .* (150 / (14 - 4));
BitSize_scaled = (6.0 - 4) * (150 / (14 - 4)); % Bit Size = 6.0 inch

for i = 1:length(DEPT)-1
    if isnan(CALI_scaled(i)) || isnan(CALI_scaled(i+1)); continue; end
    if CALI_scaled(i) > BitSize_scaled
        patch([BitSize_scaled CALI_scaled(i) CALI_scaled(i) BitSize_scaled], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], [1 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
end
plot(CALI_scaled, DEPT, 'k-', 'LineWidth', 2.0); % Garis Caliper Hitam Tebal
xline(BitSize_scaled, 'k--', 'LineWidth', 2.0); % Garis Putus-putus Bit Size 6 inch
% --------------------------------------

xlim([0 150]); xticks([0 50 100 150]);
ylabel('Depth (ft)', 'FontWeight', 'bold', 'FontSize', 11); set(ax_tc(1), 'Layer', 'top');

% Track 2: Resistivity (ZONA PROSPEK)
axes(ax_tc(2)); LLD_plot = LLD; LLS_plot = LLS; LLD_plot(LLD_plot<=0)=NaN; LLS_plot(LLS_plot<=0)=NaN; 
semilogx(LLD_plot, DEPT, 'r', 'LineWidth', 1.5); semilogx(LLS_plot, DEPT, 'm', 'LineWidth', 1.0); 
set(gca,'XScale','log'); xlim([0.2 2000]); xticks([0.2 1 10 100 1000 2000]); set(gca,'YTickLabel',[]);

% Track 3: Density-Neutron (ZONA PROSPEK)
axes(ax_tc(3)); NPHI_scaled = 2.7 - (NPHI .* (1.0/0.6)); separation = NPHI_scaled - RHOB;
for i = 1:length(DEPT)-1
    if isnan(RHOB(i)) || isnan(NPHI_scaled(i)); continue; end
    if separation(i) > 0, patch([RHOB(i) NPHI_scaled(i) NPHI_scaled(i) RHOB(i)], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], [1 1 0], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    else, patch([RHOB(i) NPHI_scaled(i) NPHI_scaled(i) RHOB(i)], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], [0 0.6 0], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); end
end
plot(RHOB, DEPT, 'r', 'LineWidth', 1.5); plot(NPHI_scaled, DEPT, 'b--', 'LineWidth', 1.5); 
xlim([1.7 2.7]); xticks([1.7 1.95 2.2 2.45 2.7]); set(gca,'YTickLabel',[]); set(ax_tc(3), 'Layer', 'top'); 

linkaxes(ax_tc,'y'); set(ax_tc(1), 'YLim', [Top_Formasi Base_Formasi]); 
exportgraphics(fig_tc, '1_Triple_Combo_Target_600DPI.png', 'Resolution', 600);

%% ========================================================================
%   4. IDENTIFIKASI LITOLOGI (INPUT 3 ZONA & VISUALISASI CROSSPLOT)
% =========================================================================
prompt_lit = {
    'Zona 1 Top (ft) [Wajib Diisi]:', ...
    'Zona 1 Base (ft) [Wajib Diisi]:', ...
    'Zona 2 Top (ft) [Ketik 0 jika reservoir tipis & tidak dibagi]:', ...
    'Zona 2 Base (ft) [Ketik 0 jika reservoir tipis & tidak dibagi]:', ...
    'Zona 3 Top (ft) [Ketik 0 jika reservoir tipis & tidak dibagi]:', ...
    'Zona 3 Base (ft) [Ketik 0 jika reservoir tipis & tidak dibagi]:'
};
dlgtitle_lit = 'Input Pembagian Zona Litologi (Crossplot)'; dims = [1 75]; 
definput_lit = {'5500', '6000', '6000', '7000', '7000', '7400'};
ans_lit = inputdlg(prompt_lit, dlgtitle_lit, dims, definput_lit);
if isempty(ans_lit), error('ANALISIS BERHENTI: Input Zona Litologi dibatalkan.'); end

Zona_Kualitatif = [
    str2double(ans_lit{1}), str2double(ans_lit{2});
    str2double(ans_lit{3}), str2double(ans_lit{4});
    str2double(ans_lit{5}), str2double(ans_lit{6})
];

ZoneColors = [
    1.0 0.0 0.0; 
    0.0 0.0 1.0; 
    0.0 0.8 0.0; 
];

idx_formasi = (DEPT >= Top_Formasi & DEPT <= Base_Formasi & valid_log);
Title_FS = 22; Label_FS = 18; Tick_FS  = 14; Text_FS  = 14;   

% -------------------------------------------------------------------------
% 4.1 DENSITY-NEUTRON CROSSPLOT
% -------------------------------------------------------------------------
fig_dn = figure('Name','Density-Neutron Crossplot','Color','w','Position',[100 50 1000 950]);
hold on; box on; xlim([-0.05 0.45]); ylim([1.9 3.0]);
set(gca, 'XDir', 'normal', 'YDir', 'reverse', 'FontSize', Tick_FS, 'LineWidth', 2.0, 'FontWeight', 'bold'); 
grid on; set(gca, 'GridLineStyle', ':', 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.6);
xlabel('Neutron Porosity (v/v) - Limestone Matrix', 'FontWeight','bold', 'FontSize',Label_FS);
ylabel('Bulk Density, \rho_b (g/cc)', 'FontWeight','bold', 'FontSize',Label_FS);
title('Density-Neutron Crossplot', 'FontWeight','bold', 'FontSize',Title_FS);

phi_ref = 0:0.01:0.50; 
pts_n_ss=[-0.015,0.04,0.10,0.16,0.22,0.28,0.34]; pts_d_ss=[2.65,2.56,2.46,2.36,2.26,2.16,2.06]; 
n_line_ss=interp1(pts_n_ss,pts_n_ss,phi_ref,'linear','extrap'); d_line_ss=interp1(pts_n_ss,pts_d_ss,phi_ref,'pchip','extrap');
n_line_ls=phi_ref; d_line_ls=(1-phi_ref)*2.71+phi_ref*1.0; 
pts_n_dl=[0.01,0.07,0.13,0.19,0.25,0.31,0.38]; pts_d_dl=[2.87,2.76,2.65,2.54,2.43,2.32,2.20];
n_line_dl=interp1(pts_n_dl,pts_n_dl,phi_ref,'linear','extrap'); d_line_dl=interp1(pts_n_dl,pts_d_dl,phi_ref,'pchip','extrap');

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
legend_names = cell(1, size(Zona_Kualitatif, 1));
for i = 1:size(Zona_Kualitatif, 1)
    idx_z = (DEPT >= Zona_Kualitatif(i,1) & DEPT <= Zona_Kualitatif(i,2)) & valid_log & ~isnan(NPHI) & ~isnan(RHOB);
    h_zones(i) = scatter(NPHI(idx_z), RHOB(idx_z), 55, ZoneColors(i,:), 'filled', 'MarkerEdgeColor','k', 'LineWidth',1.0, 'MarkerFaceAlpha',0.9);
    legend_names{i} = sprintf('Z%d (%d-%d)', i, Zona_Kualitatif(i,1), Zona_Kualitatif(i,2));
end
legend(h_zones, legend_names, 'Location','southwest', 'FontSize', Text_FS, 'FontWeight','bold');
plot(-0.01, 2.98, 'ko', 'MarkerFaceColor','w', 'LineWidth',2.0, 'MarkerSize',10); text(0.01, 2.98, 'Anhydrite', 'FontSize',Text_FS, 'FontWeight','bold');
plot(1.0, 1.0, 'bo', 'MarkerFaceColor','b', 'MarkerSize',10); text(0.40, 1.95, '\rho_f = 1.0 g/cc', 'FontSize',Text_FS, 'FontWeight','bold', 'Color','b'); hold off;
exportgraphics(fig_dn, '1_Density_Neutron_Crossplot_600DPI.png', 'Resolution', 600);

% -------------------------------------------------------------------------
% 4.2 MATRIX IDENTIFICATION PLOT (MID PLOT)
% -------------------------------------------------------------------------
fig_mid = figure('Name','MID Plot','Color','w','Position',[50 50 1100 1000]); 
axMain = axes('Position', [0.12 0.12 0.82 0.78]); hold(axMain, 'on'); box(axMain, 'on');
xlim(axMain, [35 75]); ylim(axMain, [2.4 3.1]); 
set(axMain, 'YDir', 'reverse', 'FontSize', Tick_FS, 'LineWidth', 2.0, 'FontWeight', 'bold');
set(axMain, 'XTick', 35:5:75, 'YTick', 2.4:0.1:3.1);  
grid(axMain, 'on'); set(axMain, 'XMinorGrid', 'on', 'YMinorGrid', 'on'); 
set(axMain, 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.6, 'GridLineStyle', ':');

rho_f_mid=1.00; dt_f_mid=189;
rho_maa = (RHOB - NPHI.*rho_f_mid) ./ (1 - NPHI); dt_maa  = (DT - NPHI.*dt_f_mid) ./ (1 - NPHI);
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
    idx_z_mid = (DEPT >= Zona_Kualitatif(i,1) & DEPT <= Zona_Kualitatif(i,2)) & valid_log & ~isnan(rho_maa) & ~isnan(dt_maa);
    scatter(axMain, dt_maa(idx_z_mid), rho_maa(idx_z_mid), 55, ZoneColors(i,:), 'filled', 'MarkerEdgeColor','k', 'LineWidth',1.0, 'MarkerFaceAlpha',0.9);
end
xlabel(axMain, '\Delta t_{maa}, Apparent Matrix Transit Time (\mus/ft)', 'FontWeight','bold', 'FontSize',Label_FS);
ylabel(axMain, '\rho_{maa}, Apparent Matrix Density (g/cc)', 'FontWeight','bold', 'FontSize',Label_FS);
title(axMain, 'Matrix Identification Plot (MID)', 'FontWeight','bold', 'FontSize',Title_FS);
legend(axMain, legend_names, 'Location','southwest', 'FontSize', Text_FS, 'FontWeight','bold'); hold(axMain, 'off');
exportgraphics(fig_mid, '2_MID_Plot_600DPI.png', 'Resolution', 600);

% -------------------------------------------------------------------------
% 4.3 SONIC–DENSITY CROSSPLOT 
% -------------------------------------------------------------------------
fig_sd = figure('Name','Sonic-Density Crossplot','Color','w','Position',[50 50 1100 1000]); hold on; box on;
set(gca, 'XLim', [40 120], 'YLim', [1.9 3.0], 'YDir', 'reverse', 'FontSize', Tick_FS, 'LineWidth', 2.0, 'FontWeight', 'bold');
set(gca, 'XTick', 40:10:120); set(gca, 'YTick', 1.90:0.1:3.00); set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; set(gca, 'GridLineStyle', '-', 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.6);
xlabel('Interval Transit Time, \Delta t (\mus/ft)', 'FontSize',Label_FS, 'FontWeight','bold');
ylabel('Bulk Density, \rho_b (g/cc)', 'FontSize',Label_FS, 'FontWeight','bold');
title('Sonic–Density Crossplot', 'FontSize',Title_FS, 'FontWeight','bold');

rho_f=1.00; dt_f=189; QZ=[2.65 55.5]; LS=[2.71 47.6]; DL=[2.87 43.5]; AN=[2.98 50.0]; HL=[2.03 67.0];
phi_fine=(0:0.1:45)/100; phi_tick=(0:1:40)/100;  
CalcRho = @(rho_ma, phi) (1 - phi) .* rho_ma + phi .* rho_f;
CalcDt_Wyllie = @(dt_ma, phi) (1 - phi) .* dt_ma + phi .* dt_f; CalcDt_RHG = @(dt_ma, phi) dt_ma ./ (1 - 1.6 * phi);
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
    plot(dt_ma, rho_ma, 'ko', 'MarkerFaceColor','w', 'MarkerSize', 12, 'LineWidth', 1.5); text(dt_ma-2.0, rho_ma-0.03, name, 'FontWeight','bold', 'FontSize',Text_FS);
end
plot(AN(2), AN(1), 'ko', 'MarkerFaceColor','w', 'MarkerSize',12, 'LineWidth', 1.5); text(AN(2)+2.0, AN(1), 'Anhydrite', 'FontSize',Text_FS, 'FontWeight','bold');
plot(HL(2), HL(1), 'ko', 'MarkerFaceColor','w', 'MarkerSize',12, 'LineWidth', 1.5); text(HL(2)+2.0, HL(1), 'Halite', 'FontSize',Text_FS, 'FontWeight','bold');

idx_bg_sd = idx_formasi & ~isnan(DT) & ~isnan(RHOB);
scatter(DT(idx_bg_sd), RHOB(idx_bg_sd), 25, [0.6 0.6 0.6], 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5);
h_zones_sd = zeros(1, size(Zona_Kualitatif, 1));
for i = 1:size(Zona_Kualitatif, 1)
    idx_z = (DEPT >= Zona_Kualitatif(i,1) & DEPT <= Zona_Kualitatif(i,2)) & valid_log & ~isnan(DT) & ~isnan(RHOB);
    h_zones_sd(i) = scatter(DT(idx_z), RHOB(idx_z), 55, ZoneColors(i,:), 'filled', 'MarkerEdgeColor','k', 'LineWidth',1.0, 'MarkerFaceAlpha',0.9);
end
h1 = plot(nan,nan, 'k-', 'LineWidth',3.0); h2 = plot(nan,nan, 'r-', 'LineWidth',2.5);
legend([h1 h2 h_zones_sd(1)], {'Empirical (Raymer-5.Hunt-Gardner)', 'Time Average (Wyllie)', 'Zona Potensi'}, 'Location', 'northeast', 'FontSize',Text_FS, 'FontWeight','bold'); hold off;
exportgraphics(fig_sd, '3_Sonic_Density_Crossplot_600DPI.png', 'Resolution', 600);

% -------------------------------------------------------------------------
% 4.4 NEUTRON-SONIC CROSSPLOT 
% -------------------------------------------------------------------------
fig_ns = figure('Name','Neutron-Sonic Crossplot','Color','w','Position',[50 20 1100 1100]); 
ax = axes('Position',[0.12 0.1 0.8 0.8]); hold on; box on; set(ax, 'XTick', [], 'YTick', []); 
title('Neutron-Sonic Crossplot', 'FontWeight','bold', 'FontSize',Title_FS);
xlim([-5 45]); ylim([40 120]); 
col_major = [0.0 0.3 0.4]; col_minor = [0.0 0.6 0.7]; 
for x = -5:1:45, if mod(x,5)==0, plot([x x], [40 120], 'Color', col_major, 'LineWidth', 2.0); text(x, 38.0, num2str(x), 'Horiz','center', 'FontSize',Tick_FS, 'FontWeight','bold'); else, plot([x x], [40 120], 'Color', col_minor, 'LineWidth', 1.0, 'LineStyle', '-'); end; end
for y = 40:2:120, if mod(y,10)==0, plot([-5 45], [y y], 'Color', col_major, 'LineWidth', 2.0); text(-5.5, y, num2str(y), 'Horiz','right', 'FontSize',Tick_FS, 'FontWeight','bold'); else, plot([-5 45], [y y], 'Color', col_minor, 'LineWidth', 1.0, 'LineStyle', '-'); end; end
for ym = 140:20:380, y_ft = ym / 3.28084; if y_ft >= 40 && y_ft <= 120, plot([45 45.8], [y_ft y_ft], 'Color','k', 'LineWidth',2.0); text(46.5, y_ft, num2str(ym), 'Horiz','left', 'FontSize',12, 'Color','k', 'FontWeight','bold'); end; end

dt_qz=55.5; dt_ls=47.6; dt_dl=43.5; dt_fl=189; phi=0:0.1:46; phi_pct=phi;
w_qz=(1-phi/100)*dt_qz+(phi/100)*dt_fl; w_ls=(1-phi/100)*dt_ls+(phi/100)*dt_fl; w_dl=(1-phi/100)*dt_dl+(phi/100)*dt_fl;
pts_n=[0,10,20,30,40,45]; pts_d_qz=[55.5,67.5,84.0,106.0,138.0,158.0]; emp_qz=spline(pts_n,pts_d_qz,phi_pct);
pts_d_ls=[47.6,58.5,73.0,93.0,120.0,138.0]; emp_ls=spline(pts_n,pts_d_ls,phi_pct); pts_d_dl=[43.5,53.5,66.5,84.0,108.0,124.0]; emp_dl=spline(pts_n,pts_d_dl,phi_pct);
plot(phi_pct, w_qz, 'r-', 'LineWidth', 2.5); plot(phi_pct, w_ls, 'r-', 'LineWidth', 2.5); plot(phi_pct, w_dl, 'r-', 'LineWidth', 2.5); 
plot(phi_pct, emp_qz, 'k-', 'LineWidth', 3.0); plot(phi_pct, emp_ls, 'k-', 'LineWidth', 3.0); plot(phi_pct, emp_dl, 'k-', 'LineWidth', 3.0);
text(11, 55.5 + 11*1.335, '  Quartz  ', 'Rotation',45, 'FontSize',Text_FS, 'FontWeight','bold', 'BackgroundColor','w', 'EdgeColor','none', 'HorizontalAlignment', 'center'); 
text(16, 47.6 + 16*1.414, '  Calcite  ', 'Rotation',47, 'FontSize',Text_FS, 'FontWeight','bold', 'BackgroundColor','w', 'EdgeColor','none', 'HorizontalAlignment', 'center'); 
text(21, 43.5 + 21*1.455, '  Dolomite  ', 'Rotation',48, 'FontSize',Text_FS, 'FontWeight','bold', 'BackgroundColor','w', 'EdgeColor','none', 'HorizontalAlignment', 'center');
plot(0, dt_qz, 'ko','MarkerFaceColor','w', 'MarkerSize', 10, 'LineWidth', 1.5); plot(0, dt_ls, 'ko','MarkerFaceColor','w', 'MarkerSize', 10, 'LineWidth', 1.5); plot(0, dt_dl, 'ko','MarkerFaceColor','w', 'MarkerSize', 10, 'LineWidth', 1.5);
rectangle('Position', [1 103 26 15], 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 2.0); 
plot([3 7], [113 113], 'r-', 'LineWidth', 3.0); text(9, 113, 'Time Average', 'FontSize',Text_FS, 'FontWeight','bold'); plot([3 7], [107 107], 'k-', 'LineWidth', 3.0); text(9, 107, 'Empirical', 'FontSize',Text_FS, 'FontWeight','bold');
rectangle('Position', [22 43 23 7], 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5); text(33.5, 46.5, ['\Delta t_f = ' num2str(dt_fl) ' \mu s/ft'], 'Horiz','center', 'FontSize',Text_FS, 'FontWeight','bold');

NPHI_pct = NPHI * 100; idx_bg_ns = idx_formasi & ~isnan(DT) & ~isnan(NPHI);
scatter(NPHI_pct(idx_bg_ns), DT(idx_bg_ns), 25, [0.6 0.6 0.6], 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5);
for i = 1:size(Zona_Kualitatif, 1)
    idx_z = (DEPT >= Zona_Kualitatif(i,1) & DEPT <= Zona_Kualitatif(i,2)) & valid_log & ~isnan(DT) & ~isnan(NPHI);
    scatter(NPHI_pct(idx_z), DT(idx_z), 55, ZoneColors(i,:), 'filled', 'MarkerEdgeColor','k', 'LineWidth',1.0, 'MarkerFaceAlpha',0.9);
end
hold off; exportgraphics(fig_ns, '4_Neutron_Sonic_Crossplot_600DPI.png', 'Resolution', 600);

% -------------------------------------------------------------------------
% 4.5 TABEL KESIMPULAN LITOLOGI
% -------------------------------------------------------------------------
calc_rho_from_dt = @(dt_log, rho_m, dt_m) (1 - (0.625 * (1 - dt_m./dt_log))) * rho_m + (0.625 * (1 - dt_m./dt_log)) * 1.0;
fprintf('\n==========================================================================================================================\n');
fprintf('                              COMPARISON OF LITHOLOGY ESTIMATION FROM VARIOUS CROSSPLOTS\n');
fprintf('==========================================================================================================================\n');
fprintf('        |                 |                                       Lithology                                      |\n');
fprintf('        |                 |--------------------------------------------------------------------------------------|\n');
fprintf('  Zone  |   Zone Range    |   Neutron-Density   |    Sonic-Density    |    Neutron-Sonic    |      MID Plot      |\n');
fprintf('        |      (ft)       |      Crossplot      |      Crossplot      |      Crossplot      |                    |\n');
fprintf('--------------------------------------------------------------------------------------------------------------------------\n');

for i = 1:size(Zona_Kualitatif, 1)
    z_top = Zona_Kualitatif(i,1); z_bot = Zona_Kualitatif(i,2);
    idx = (DEPT >= z_top & DEPT <= z_bot) & valid_log & ~isnan(RHOB) & ~isnan(NPHI) & ~isnan(DT);
    
    if sum(idx) == 0
        fprintf('   %d    |  %4d - %4d   |       NO DATA       |       NO DATA       |       NO DATA       |      NO DATA       |\n', i, z_top, z_bot); continue;
    end
    
    N_m = mean(NPHI(idx)); R_m = mean(RHOB(idx)); D_m = mean(DT(idx));
    
    % Litologi NPHI-RHOB
    rho_maa_val = (R_m - N_m*1.0) / (1 - N_m);
    if rho_maa_val < 2.68, lit_ND = 'Gas Effect'; elseif rho_maa_val < 2.80, lit_ND = 'Limestone'; else, lit_ND = 'Dolomite'; end
    
    % Litologi DT-RHOB
    rho_QZ = calc_rho_from_dt(D_m, 2.65, 55.5); rho_LS = calc_rho_from_dt(D_m, 2.71, 47.6); rho_DL = calc_rho_from_dt(D_m, 2.87, 43.5);
    [~, min_id_SD] = min(abs(R_m - [rho_QZ, rho_LS, rho_DL])); nama_lit_SD = {'Gas Effect', 'Limestone', 'Dolomite'}; 
    if R_m < (rho_QZ - 0.05), lit_SD = 'Gas Effect'; else, lit_SD = nama_lit_SD{min_id_SD}; end
    
    % Litologi NPHI-DT
    dt_maa_val = (D_m - N_m*189) / (1 - N_m);
    if dt_maa_val > 51, lit_NS = 'Gas Effect'; elseif dt_maa_val > 44, lit_NS = 'Limestone'; else, lit_NS = 'Dolomite'; end
    
    % Litologi MID
    d1 = sqrt((dt_maa_val - 55.5)^2 + (50*(rho_maa_val - 2.65))^2); d2 = sqrt((dt_maa_val - 47.6)^2 + (50*(rho_maa_val - 2.71))^2); d3 = sqrt((dt_maa_val - 43.5)^2 + (50*(rho_maa_val - 2.87))^2);
    [~, min_id_MID] = min([d1, d2, d3]);
    if rho_maa_val < 2.60 || dt_maa_val > 58, lit_MID = 'Gas Effect'; else, lit_MID = nama_lit_SD{min_id_MID}; end
    
    fprintf('   %d    |  %4d - %4d   | %-19s | %-19s | %-19s | %-18s |\n', i, z_top, z_bot, lit_ND, lit_SD, lit_NS, lit_MID);
end
fprintf('==========================================================================================================================\n');

%% ========================================================================
%   5. PERHITUNGAN PETROFISIKA INTI, PICKETT PLOT & STATISTIK CUT-OFF
% =========================================================================
pause(2); 
prompt_petro = {
    'Densitas Matriks (rho_ma) [Sandstone = 2.65, Limestone = 2.71, Dolomite = 2.87]:', ...
    'Densitas Fluida (rho_f) [Fresh Water = 1.0, Salt Water = 1.1, Oil = 0.85, ]:', ...
    'Parameter Archie - a (Tortuosity Factor) [Sandstone = 0.62 , Karbonat = 1.0]:', ...
    'Parameter Archie - m (Cementation Factor) [Sandstone = 2.15, Karbonat = 2.0]:', ...
    'Parameter Archie - n (Saturation Exponent) [Umumnya 2.0]:', ...
    'Top Zona Air (ft) [Untuk menghitung Rw murni dari log, misal: 7200]:', ...
    'Base Zona Air (ft) [misal: 7250]:'
};
dlgtitle_petro = 'Input Parameter Petrofisika & Zona Air'; dims = [1 85];
definput_petro = {'2.71', '1.0', '1.0', '2.0', '2.0', '7200', '7250'};
ans_petro = inputdlg(prompt_petro, dlgtitle_petro, dims, definput_petro);
if isempty(ans_petro), error('ANALISIS BERHENTI: Input Petrofisika dibatalkan oleh pengguna.'); end
rho_ma = str2double(ans_petro{1});
rho_f  = str2double(ans_petro{2});
a = str2double(ans_petro{3});
m = str2double(ans_petro{4});
n = str2double(ans_petro{5});
z_top_water = str2double(ans_petro{6});
z_bot_water = str2double(ans_petro{7});

% --- PERHITUNGAN BASELINE & VSHALE ---
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
fprintf('\n=> Rho_sh = %.4f g/cc, NPHI_sh = %.4f v/v, Rsh = %.4f ohm.m\n', rho_sh, phi_Nsh, Rsh);
GR_valid = GR(~isnan(GR) & GR > 0);
GR_clean = prctile(GR_valid, 5); 
GR_shale = prctile(GR_valid, 95); 
Vsh = (GR - GR_clean) ./ (GR_shale - GR_clean); 
% SAFETY CLAMP 
Vsh(Vsh < 0) = 0; Vsh(Vsh > 1) = 1; 

% --- PERHITUNGAN POROSITAS EFEKTIF ---
phi_D = (rho_ma - RHOB) ./ (rho_ma - rho_f); 
phi_Dsh = (rho_ma - rho_sh) ./ (rho_ma - rho_f);
fprintf('\n========================================================\n');
fprintf('                     PARAMETER\n');
fprintf('==========================================================\n');
fprintf(' - Densitas Shale (pb_sh / rho_sh)   : %.4f g/cc\n', rho_sh);
fprintf(' - Porositas Neutron Shale (Por_Nsh) : %.4f v/v\n', phi_Nsh);
fprintf(' - Porositas Density Shale (phi_Dsh) : %.4f v/v\n', phi_Dsh);
fprintf(' - Clean Sand Baseline (GR Min)      : %.2f API\n', GR_clean);
fprintf(' - Shale Baseline (GR Max)           : %.2f API\n', GR_shale);
fprintf('========================================================\n\n');
phi_Dc = phi_D - (Vsh .* phi_Dsh); 
phi_Nc = NPHI - (Vsh .* phi_Nsh); 
phi_eff = sqrt((phi_Nc.^2 + phi_Dc.^2) ./ 2);
% SAFETY CLAMP
phi_eff(phi_eff < 0) = 0; phi_eff(phi_eff > 1) = 1; 

% --- PERHITUNGAN WATER SATURATION & Rw ---
idx_water = (DEPT >= z_top_water) & (DEPT <= z_bot_water) & ~isnan(LLD) & ~isnan(phi_eff);
Rwa_water_zone = (LLD(idx_water) .* (phi_eff(idx_water).^m)) ./ a;
Rw = median(Rwa_water_zone, 'omitnan');
if isnan(Rw)
    error('ANALISIS BERHENTI: Gagal menghitung Rw. Tidak ada data log valid di zona air!');
end
fprintf('=> Nilai Rw : %.6f ohm.m\n', Rw);
Sw = ( (a .* Rw) ./ ((phi_eff.^m) .* LLD) ) .^ (1/n);
% SAFETY CLAMP
Sw(Sw < 0) = 0; Sw(Sw > 1) = 1; 

% --- VISUALISASI PICKETT PLOT ---
fig_pickett = figure('Name', 'Pickett Plot - Final Validation', 'Color', 'w', 'Position', [200 200 900 800]);
h_all = loglog(LLD(zone), phi_eff(zone), 'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none'); 
hold on; box on; 
phi_line_range = logspace(log10(0.01), log10(1.0), 100); 
sw_lines = [1.0, 0.5, 0.25]; colors = {'b', 'g', 'r'};   
for i = 1:length(sw_lines)
    sw_val = sw_lines(i);
    rt_calc = (a * Rw) ./ ((phi_line_range.^m) .* (sw_val^n));
    plot(rt_calc, phi_line_range, 'Color', colors{i}, 'LineWidth', 3.0); 
    text_phi = 0.05; text_rt = (a * Rw) / (text_phi^m * sw_val^n);
    if text_rt < 2000 
        text(text_rt, text_phi, [' Sw ' num2str(sw_val*100) '%'], 'Color', colors{i}, 'FontWeight', 'bold', ...
            'FontSize', 14, 'BackgroundColor', 'w', 'EdgeColor', colors{i}, 'Margin', 2);
    end
end
if sum(idx_water) > 0
    h_water = loglog(LLD(idx_water), phi_eff(idx_water), 'bo', 'MarkerSize', 9, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2); 
end
plot(Rw, 1, 'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'r', 'LineWidth', 1.5); 
title(['Pickett Plot | Rw = ' num2str(Rw, '%.6f') ' ohm.m'], 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Deep Resistivity, R_t (ohm.m)', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Effective Porosity, \phi_e (%)', 'FontWeight', 'bold', 'FontSize', 16);
set(gca, 'XLim', [0.1 2000], 'YLim', [0.01 1.0], 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2.5, 'TickDir', 'in'); 
grid on; 
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridLineStyle', '--', 'GridAlpha', 0.6, 'MinorGridLineStyle', ':', 'MinorGridAlpha', 0.4);
if sum(idx_water) > 0
      legend([h_all, h_water], {'All Data', 'Water Zone (Verified)'}, 'Location', 'southwest', 'FontSize', 13, 'LineWidth', 1.5);
else
    legend(h_all, {'All Data'}, 'Location', 'southwest', 'FontSize', 13, 'LineWidth', 1.5);
end
hold off;
exportgraphics(fig_pickett, '1_Pickett_Plot_600DPI.png', 'Resolution', 600);

% --- INPUT CUT-OFF & FLAG IDENTIFICATION ---
pause(2); 
prompt_co = {
    'Cut-off Volume Shale (CO_Vsh) [%, misal: 0.37]:', ...
    'Cut-off Porositas Efektif (CO_Phi) [%, misal: 0.065]:', ...
    'Cut-off Water Saturation (CO_Sw) [%, misal: 0.55]:'
};
dlgtitle_co = 'Input Nilai Cut-Off Formasi'; dims = [1 65];
definput_co = {'0.37', '0.065', '0.55'};
ans_co = inputdlg(prompt_co, dlgtitle_co, dims, definput_co);
if isempty(ans_co), error('ANALISIS BERHENTI: Input Cut-Off dibatalkan.'); end
CO_Vsh = str2double(ans_co{1});
CO_Phi = str2double(ans_co{2});
CO_Sw  = str2double(ans_co{3});

% LOGICAL SHADING FLAGS
netResFlag = (Vsh <= CO_Vsh) & (phi_eff >= CO_Phi) & valid_log;
netPayFlag = netResFlag & (Sw <= CO_Sw) & valid_log;
netResFlag(isnan(Vsh) | isnan(phi_eff)) = 0;
netPayFlag(isnan(Vsh) | isnan(phi_eff) | isnan(Sw)) = 0;

vsh_form = Vsh(zone); 
phi_form = phi_eff(zone); 
sw_form  = Sw(zone);

% --- VISUALISASI ANALISIS CUT-OFF STATISTIK ---
Title_FS = 20; Label_FS = 16; Tick_FS  = 14; Text_FS  = 15;
fig_cutoff = figure('Name', 'Analisis Cut-Off Statistik', 'Color', 'w', 'Position', [50 100 1600 550], 'NumberTitle', 'off');

% Subplot 1: Histogram Porositas
subplot(1,3,1);
histogram(phi_form(vsh_form <= CO_Vsh), 30, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k', 'LineWidth', 1.2); hold on; grid on;
xline(CO_Phi, 'k--', 'LineWidth', 5.0); 
title('Histogram PHIE', 'FontWeight', 'bold', 'FontSize', Title_FS);
xlabel('Effective Porosity (PHIE)', 'FontWeight', 'bold', 'FontSize', Label_FS); ylabel('Frekuensi', 'FontWeight', 'bold', 'FontSize', Label_FS);
xlim([0 max(0.4, max(phi_form)+0.05)]);
text(CO_Phi + 0.015, max(ylim)*0.50, sprintf('PHIE \\geq %.3f', CO_Phi), 'Color', 'k', 'FontWeight', 'bold', 'FontSize', Text_FS, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, 'Margin', 3);
set(gca, 'FontSize', Tick_FS, 'LineWidth', 1.5, 'FontWeight', 'bold', 'GridLineStyle', ':', 'GridAlpha', 0.6);

% Subplot 2: Crossplot VSH vs PHIE
subplot(1,3,2);
scatter(vsh_form, phi_form, 45, DEPT(zone), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.8); colormap(jet); hold on; grid on;
xline(CO_Vsh, 'k--', 'LineWidth', 5.0); yline(CO_Phi, 'k--', 'LineWidth', 5.0); 
title('Crossplot VSH vs PHIE', 'FontWeight', 'bold', 'FontSize', Title_FS);
xlabel('Volume Shale (VSH)', 'FontWeight', 'bold', 'FontSize', Label_FS); ylabel('Effective Porosity (PHIE)', 'FontWeight', 'bold', 'FontSize', Label_FS);
xlim([0 1.0]); ylim([0 max(0.4, max(phi_form)+0.05)]);
text(CO_Vsh + 0.03, max(ylim)*0.50, sprintf('VSH \\leq %.2f', CO_Vsh), 'FontWeight', 'bold', 'FontSize', Text_FS, 'Color', 'k', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, 'Margin', 3);
text(max(xlim)*0.70, CO_Phi + 0.02, sprintf('PHIE \\geq %.3f', CO_Phi), 'FontWeight', 'bold', 'FontSize', Text_FS, 'Color', 'k', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, 'Margin', 3, 'HorizontalAlignment', 'center');
set(gca, 'FontSize', Tick_FS, 'LineWidth', 1.5, 'FontWeight', 'bold', 'GridLineStyle', ':', 'GridAlpha', 0.6);

% Subplot 3: Crossplot PHIE vs SW
subplot(1,3,3);
scatter(phi_form, sw_form, 45, DEPT(zone), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.8); hold on; grid on;
xline(CO_Phi, 'k--', 'LineWidth', 5.0); yline(CO_Sw, 'k--', 'LineWidth', 5.0);  
title('Crossplot PHIE vs SW', 'FontWeight', 'bold', 'FontSize', Title_FS);
xlabel('Effective Porosity (PHIE)', 'FontWeight', 'bold', 'FontSize', Label_FS); ylabel('Water Saturation (SW)', 'FontWeight', 'bold', 'FontSize', Label_FS);
xlim([0 max(0.4, max(phi_form)+0.05)]); ylim([0 1.0]);
text(CO_Phi + 0.02, 0.85, sprintf('PHIE \\geq %.3f', CO_Phi), 'FontWeight', 'bold', 'FontSize', Text_FS, 'Color', 'k', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, 'Margin', 3);
text(max(xlim)*0.70, CO_Sw + 0.04, sprintf('SW \\leq %.2f', CO_Sw), 'FontWeight', 'bold', 'FontSize', Text_FS, 'Color', 'k', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, 'Margin', 3, 'HorizontalAlignment', 'center');
set(gca, 'FontSize', Tick_FS, 'LineWidth', 1.5, 'FontWeight', 'bold', 'GridLineStyle', ':', 'GridAlpha', 0.6);
exportgraphics(fig_cutoff, '2_Analisis_Cutoff_Statistik_600DPI.png', 'Resolution', 600);

%% ========================================================================
%   6. PLOT DISTRIBUSI PETROFISIKA (TANPA FLAG WARNA)
% =========================================================================
fig_petro = figure('Name','Distribusi Petrofisika (Mentah)','Color','w','Position',[100 100 1000 850], 'NumberTitle','off');
margin_left = 0.08; margin_right = 0.05; margin_bottom = 0.05; 
h_header = 0.12; h_track = 0.82; gap = 0.005; y_header = margin_bottom + h_track; y_track = margin_bottom;
font_name = 'Arial'; grid_color = [0.75 0.75 0.75]; axis_size = 9;
total_width_petro = 1.0 - margin_left - margin_right; w_pet = (total_width_petro - (2 * gap)) / 3;
lbl_pet = {{'VSHALE', '(%)', '0', '1'}, {'EFFECTIVE POROSITY', '(%)', '0', '1'}, {'WATER SATURATION', '(%)', '0', '1'}};

ax_pet = gobjects(1, 3); curr_x = margin_left;
for i = 1:3
    axes('Position',[curr_x, y_header, w_pet, h_header]); axis off; rectangle('Position',[0 0 1 1],'LineWidth',1); 
    text(0.5, 0.70, lbl_pet{i}{1}, 'Horiz','center','FontWeight','bold','FontSize',14);
    text(0.5, 0.30, lbl_pet{i}{2}, 'Horiz','center','FontSize',15,'FontWeight','bold');
    text(0.02, 0.1, lbl_pet{i}{3}, 'Horiz','left','FontSize',12,'FontWeight','bold'); text(0.98, 0.1, lbl_pet{i}{4}, 'Horiz','right','FontSize',12,'FontWeight','bold');
    
    ax_pet(i) = axes('Position', [curr_x, y_track, w_pet, h_track]);
    set(ax_pet(i), 'FontName', font_name, 'FontSize', axis_size, 'Box', 'on', 'LineWidth', 1.0, 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k', 'YDir', 'reverse', 'XAxisLocation', 'bottom', 'XTickLabel', [], 'XGrid', 'on', 'YGrid', 'on', 'GridColor', grid_color, 'GridAlpha', 0.8, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'MinorGridColor', grid_color, 'MinorGridAlpha', 0.4);
    if i == 1, ylabel(ax_pet(i), 'Depth (ft)', 'FontWeight', 'bold', 'FontSize', 11); else, set(ax_pet(i), 'YTickLabel', []); end
    hold(ax_pet(i), 'on'); curr_x = curr_x + w_pet + gap;
end

% Track Vshale
axes(ax_pet(1)); 
for i = 1:length(DEPT)-1
    if isnan(Vsh(i)) || isnan(Vsh(i+1)); continue; end
    patch([Vsh(i) 1 1 Vsh(i)], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], [0.85 0.85 0.85], 'FaceAlpha', 1.0, 'EdgeColor', 'none'); 
end
plot(Vsh, DEPT, 'Color', [0.4 0.3 0.1], 'LineWidth', 1.5); 
xlim([0 1]); xticks([0 0.25 0.5 0.75 1]); set(ax_pet(1), 'Layer', 'top');

% Track Porosity
axes(ax_pet(2)); plot(phi_eff, DEPT, 'c', 'LineWidth', 1.5); xlim([0 1]); xticks([0 0.25 0.5 0.75 1]); set(gca, 'XDir', 'normal');

% Track Sw
axes(ax_pet(3)); plot(Sw, DEPT, 'y', 'LineWidth', 1.5); xlim([0 1]); xticks([0 0.25 0.5 0.75 1]);

linkaxes(ax_pet,'y'); set(ax_pet(1), 'YLim', [Top_Formasi Base_Formasi]);
exportgraphics(fig_petro, '3_Distribusi_Petrofisika_Mentah_600DPI.png', 'Resolution', 600);

%% ========================================================================
%   7. PLOT FINAL (DENGAN FLAG WARNA), PAY SUMMARY & EXPORT EXCEL
% =========================================================================
fig_final = figure('Name','Final Petrophysical Distribution','Color','w','Position',[50 50 1000 850], 'NumberTitle','off');
total_width_fin = 1.0 - margin_left - margin_right; w_fin = (total_width_fin - (2 * gap)) / 3; % Kembali ke 3 track
lbl_fin = {{'VSHALE', '(%)', '0', '1'}, {'EFFECTIVE POROSITY', '(%)', '0', '1'}, {'WATER SATURATION', '(%)', '0', '1'}};

ax_fin = gobjects(1, 3); curr_x = margin_left;
for i = 1:3
    axes('Position',[curr_x, y_header, w_fin, h_header]); axis off; rectangle('Position',[0 0 1 1],'LineWidth',1); 
    text(0.5, 0.70, lbl_fin{i}{1}, 'Horiz','center','FontWeight','bold','FontSize',14);
    text(0.5, 0.30, lbl_fin{i}{2}, 'Horiz','center','FontSize',14,'FontWeight','bold');
    text(0.02, 0.1, lbl_fin{i}{3}, 'Horiz','left','FontSize',12,'FontWeight','bold'); text(0.98, 0.1, lbl_fin{i}{4}, 'Horiz','right','FontSize',12,'FontWeight','bold');
    
    ax_fin(i) = axes('Position', [curr_x, y_track, w_fin, h_track]);
    set(ax_fin(i), 'FontName', font_name, 'FontSize', axis_size, 'Box', 'on', 'LineWidth', 1.0, 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k', 'YDir', 'reverse', 'XAxisLocation', 'bottom', 'XTickLabel', [], 'XGrid', 'on', 'YGrid', 'on', 'GridColor', grid_color, 'GridAlpha', 0.8, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'MinorGridColor', grid_color, 'MinorGridAlpha', 0.4);
    if i == 1, ylabel(ax_fin(i), 'Depth (ft)', 'FontWeight', 'bold', 'FontSize', 11); else, set(ax_fin(i), 'YTickLabel', []); end
    hold(ax_fin(i), 'on'); curr_x = curr_x + w_fin + gap;
end

% Lebar pita flag di sisi track (0.05 = 5% dari lebar track)
flag_w = 0.05; 

% Track 1: Vshale
axes(ax_fin(1)); 
func_draw_patches_fast(ax_fin(1), DEPT, netResFlag, [0, flag_w], [0 0.8 0], 0.65);         % PITA HIJAU (Kiri)
func_draw_patches_fast(ax_fin(1), DEPT, netPayFlag, [1-flag_w, 1.0], [1 0 0], 0.75);       % PITA MERAH (Kanan)
for i = 1:length(DEPT)-1
    if isnan(Vsh(i)) || isnan(Vsh(i+1)); continue; end
    patch([Vsh(i) 1 1 Vsh(i)], [DEPT(i) DEPT(i) DEPT(i+1) DEPT(i+1)], [0.85 0.85 0.85], 'FaceAlpha', 1.0, 'EdgeColor', 'none'); 
end
plot(Vsh, DEPT, 'Color', [0.4 0.3 0.1], 'LineWidth', 1.5); 
xlim([0 1]); xticks([0 0.25 0.5 0.75 1]); set(ax_fin(1), 'Layer', 'top');

% Track 2: Porosity
axes(ax_fin(2)); 
func_draw_patches_fast(ax_fin(2), DEPT, netResFlag, [0, flag_w], [0 0.8 0], 0.65);         % PITA HIJAU (Kiri)
func_draw_patches_fast(ax_fin(2), DEPT, netPayFlag, [1-flag_w, 1.0], [1 0 0], 0.75);       % PITA MERAH (Kanan)
plot(phi_eff, DEPT, 'c', 'LineWidth', 1.5); 
xlim([0 1]); xticks([0 0.25 0.5 0.75 1]); set(ax_fin(2), 'Layer', 'top');

% Track 3: Sw
axes(ax_fin(3)); 
func_draw_patches_fast(ax_fin(3), DEPT, netResFlag, [0, flag_w], [0 0.8 0], 0.65);         % PITA HIJAU (Kiri)
func_draw_patches_fast(ax_fin(3), DEPT, netPayFlag, [1-flag_w, 1.0], [1 0 0], 0.75);       % PITA MERAH (Kanan)
plot(Sw, DEPT, 'y', 'LineWidth', 1.5); 
xlim([0 1]); xticks([0 0.25 0.5 0.75 1]); set(ax_fin(3), 'Layer', 'top');

linkaxes(ax_fin,'y'); set(ax_fin(1), 'YLim', [Top_Formasi Base_Formasi]);
exportgraphics(fig_final, '4_Distribusi_Petrofisika_Final_Flag_600DPI.png', 'Resolution', 600);

% --- EXPORT EXCEL & SUMMARY ---
idx_form = find(DEPT >= Top_Formasi & DEPT <= Base_Formasi & valid_log); 
gross_thick = length(idx_form) * dz;
idx_net_res = find(DEPT >= Top_Formasi & DEPT <= Base_Formasi & Vsh <= CO_Vsh & phi_eff >= CO_Phi & valid_log); 
net_res_thick = length(idx_net_res) * dz;
idx_net_pay = find(DEPT >= Top_Formasi & DEPT <= Base_Formasi & Vsh <= CO_Vsh & phi_eff >= CO_Phi & Sw <= CO_Sw & valid_log); 
net_pay_thick = length(idx_net_pay) * dz;

if net_pay_thick > 0
    avg_vsh_pay = mean(Vsh(idx_net_pay), 'omitnan'); 
    avg_phi_pay = mean(phi_eff(idx_net_pay), 'omitnan'); 
    avg_sw_pay  = mean(Sw(idx_net_pay), 'omitnan');
else
    avg_vsh_pay = 0; avg_phi_pay = 0; avg_sw_pay = 0;
end

fprintf('\n=========================================================================\n');
fprintf('                 RESERVOIR PAY SUMMARY \n');
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

ResStatus_Code = NaN(size(DEPT));
ResStatus_Code(Vsh > CO_Vsh) = 1;                                      
ResStatus_Code(Vsh <= CO_Vsh & phi_eff < CO_Phi) = 2;                  
ResStatus_Code(Vsh <= CO_Vsh & phi_eff >= CO_Phi & Sw > CO_Sw) = 3;    
ResStatus_Code(Vsh <= CO_Vsh & phi_eff >= CO_Phi & Sw <= CO_Sw) = 4;   
Keterangan_Status = repmat({'Unknown'}, size(ResStatus_Code));
Keterangan_Status(ResStatus_Code == 1) = {'Shale'};
Keterangan_Status(ResStatus_Code == 2) = {'Tight Reservoir'};
Keterangan_Status(ResStatus_Code == 3) = {'Water Reservoir'};
Keterangan_Status(ResStatus_Code == 4) = {'Net Pay'};

rho_maa_calc = (RHOB - NPHI .* rho_f) ./ (1 - NPHI);
Lithology_Name = repmat({'Unknown'}, size(DEPT));
Lithology_Name(Vsh > CO_Vsh) = {'Shale'}; 
idx_clean = (Vsh <= CO_Vsh);
Lithology_Name(idx_clean & rho_maa_calc < 2.60) = {'Gas Effect / Coal'};
Lithology_Name(idx_clean & rho_maa_calc >= 2.60 & rho_maa_calc < 2.68) = {'Sandstone'};
Lithology_Name(idx_clean & rho_maa_calc >= 2.68 & rho_maa_calc < 2.78) = {'Limestone'};
Lithology_Name(idx_clean & rho_maa_calc >= 2.78) = {'Dolomite'};

Depth_ft = DEPT(zone); GR_api = GR(zone); Caliper_in = CALI(zone);
Res_Deep = LLD(zone); Rho_bulk = RHOB(zone); Nphi_v = NPHI(zone); Dt_usft = DT(zone);
Vshale_v = Vsh(zone); Porosity_Eff = phi_eff(zone); Water_Sat = Sw(zone);
Net_Res_Flag = double(netResFlag(zone)); Net_Pay_Flag = double(netPayFlag(zone));
Valid_Data = double(valid_log(zone)); 
Status_Reservoir = Keterangan_Status(zone); Litologi_Utama = Lithology_Name(zone);

Tabel_Petrofisika = table(Depth_ft, GR_api, Caliper_in, Res_Deep, Rho_bulk, Nphi_v, Dt_usft, Vshale_v, Porosity_Eff, Water_Sat, Litologi_Utama, Status_Reservoir, Net_Res_Flag, Net_Pay_Flag, Valid_Data);
Tabel_Raw_LAS = table(DEPT(:), GR(:), CALI(:), LLD(:), LLS(:), RHOB(:), NPHI(:), DT(:), 'VariableNames', {'DEPTH', 'GR', 'CALI', 'LLD', 'LLS', 'RHOB', 'NPHI', 'DT'});
nama_file_excel = 'Hasil_Petrofisika_Dinamis.xlsx';
writetable(Tabel_Petrofisika, nama_file_excel, 'Sheet', 'Data_Petrofisika_Target');
writetable(Tabel_Raw_LAS, nama_file_excel, 'Sheet', 'Raw_Data_LAS_Full');

fprintf('=> SUKSES! File Excel telah dibuat: %s\n', nama_file_excel);
fprintf('   - Sheet 1: Data Petrofisika Target (%g ft - %g ft)\n', Top_Formasi, Base_Formasi);
fprintf('   - Sheet 2: Data Konversi Raw LAS (Full Depth)\n');
fprintf('=> ANALISIS SELESAI. SILAKAN CEK FOLDER MATLAB ANDA!\n');
fprintf('==========================================================================================================================\n\n');

%% ========================================================================
%                     FUNCTION BLOCKS 
% =========================================================================
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

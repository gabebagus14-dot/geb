# geb
Petrophysics
%% ==========================================================
%               PETROPHYSICAL ANALYSIS
%  ==========================================================

clc; clear; close all;

%% ==========================================================
% 1. READ LAS FILE
% ==========================================================

filename = 'WCL0003714.LAS';   
fid = fopen(filename,'r');
if fid < 0
    error('LAS file not found');
end

% Cari section ~A
line = '';
while ischar(line)
    line = fgetl(fid);
    if contains(upper(line),'~A')
        break
    end
end

% Baca ASCII data
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

% Ganti NULL LAS dengan NaN
NULL = -999.25;
ascii_data(ascii_data == NULL) = NaN;

% Assign log sesuai urutan LAS
DEPT = ascii_data(:,1);
GR   = ascii_data(:,2);
CALI = ascii_data(:,3);
LLD  = ascii_data(:,4);
LLS  = ascii_data(:,5);
DT   = ascii_data(:,7);
RHOB = ascii_data(:,8);
NPHI = ascii_data(:,9);

%% ==========================================================
% 2. QUICKLOOK RESERVOIR IDENTIFICATION
% ==========================================================

% --- (1) GR rendah
GR_min = min(GR(~isnan(GR)));
GR_max = max(GR(~isnan(GR)));
GR_cutoff = GR_min + (GR_max - GR_min)/2;
flag_GR = GR <= GR_cutoff;

% --- (2) Resistivity tinggi
Rt_cutoff = prctile(LLD(~isnan(LLD)),60);
flag_Rt = LLD >= Rt_cutoff;

% --- (3) NPHI–RHOB crossover
flag_cross = (NPHI <= 0.25) & (RHOB <= 2.35);

% --- QUICKLOOK RESERVOIR FLAG
flag_reservoir = flag_GR & flag_Rt & flag_cross;

%% ===== TAMBAHKAN BLOK INI DI SINI =====
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


%% ==========================================================
% 3. PETROPHYSICAL CALCULATION (ONLY QUICKLOOK ZONE)
% ==========================================================

% Pre-allocate
Vsh     = NaN(size(DEPT));
phi_eff = NaN(size(DEPT));
Sw      = NaN(size(DEPT));
k       = NaN(size(DEPT));

zone = flag_reservoir & ~isnan(GR) & ~isnan(RHOB) & ...
       ~isnan(DT) & ~isnan(LLD);

%% ---- 3.1 Volume Shale (Larionov Tertiary)
IGR = (GR(zone) - GR_min) ./ (GR_max - GR_min);
Vsh(zone) = 0.083 .* (2.^(3.7 .* IGR) - 1);
Vsh(Vsh < 0) = 0;
Vsh(Vsh > 1) = 1;

%% ---- 3.2 Effective Porosity
rho_ma = 2.65; 
rho_fl = 1.0;
phi_den = (rho_ma - RHOB(zone)) ./ (rho_ma - rho_fl);

dt_ma = 55; 
dt_fl = 180;
phi_sonic = (DT(zone) - dt_ma) ./ (dt_fl - dt_ma);

phi_total = (phi_den + phi_sonic) / 2;
phi_eff(zone) = phi_total .* (1 - Vsh(zone));
phi_eff(phi_eff < 0) = 0;
phi_eff(phi_eff > 0.4) = 0.4;

%% ---- 3.3 Water Saturation (Indonesia Equation)
a = 1; m = 2.15; n = 2;
Rw  = 0.08;
Rsh = 2.0;

Rt = LLD(zone);
A = (Vsh(zone).^(1 - Vsh(zone)/2)) ./ sqrt(Rsh);
B = sqrt((phi_eff(zone).^m) ./ (a*Rw));

Sw(zone) = ((1./sqrt(Rt) - A) ./ B).^(2/n);
Sw(Sw < 0) = 0;
Sw(Sw > 1) = 1;

%% ---- 3.4 Permeability (Timur Equation)

k = NaN(size(DEPT));

phi_cut_k = 0.15;
Sw_cut_k  = 0.70;

idx_k = zone & phi_eff >= phi_cut_k & Sw <= Sw_cut_k;

k(idx_k) = 8581 .* (phi_eff(idx_k).^4.4) ./ (Sw(idx_k).^2);

k(k > 1e4) = 1e4;   % <<< WAJIB


%% ==========================================================
% NEUTRON–DENSITY CROSSPLOT (RESERVOIR ZONE ONLY)
% ==========================================================

idxND = flag_reservoir & ~isnan(NPHI) & ~isnan(RHOB);

if sum(idxND)==0
    warning('Tidak ada data reservoir untuk N–D crossplot');
else

    phiN = NPHI(idxND) * 100;     % φN,LS (%)
    rhoB = RHOB(idxND);

    phi = linspace(0,45,200);     % porosity range (%)

    rho_fl = 1.00;               % fresh water
    rho_qz = 2.65;               % quartz
    rho_ls = 2.71;               % calcite
    rho_dl = 2.87;               % dolomite

    rho_qz_line = (1 - phi/100).*rho_qz + (phi/100).*rho_fl;
    rho_ls_line = (1 - phi/100).*rho_ls + (phi/100).*rho_fl;
    rho_dl_line = (1 - phi/100).*rho_dl + (phi/100).*rho_fl;

    figure('Color','k','Position',[120 120 750 650]); hold on; grid on

    scatter(phiN, rhoB, 28, DEPT(idxND), 'filled', ...
        'MarkerEdgeColor','none','MarkerFaceAlpha',0.8);

    plot(phi, rho_qz_line,'r','LineWidth',2)
    plot(phi, rho_ls_line,'b','LineWidth',2)
    plot(phi, rho_dl_line,'g','LineWidth',2)

    xlabel('\phi_{N,LS} Neutron Limestone Porosity (%)','Color','w')
    ylabel('Bulk Density RHOB (g/cc)','Color','w')
    title('Neutron–Density Crossplot (Reservoir Zone Only)','Color','w')

    xlim([0 45])
    ylim([1.9 3.0])
    set(gca,'YDir','reverse','Color','k','XColor','w','YColor','w')

    colormap(jet)
    cb = colorbar;
    cb.Label.String = 'Depth (ft)';
    cb.Color = 'w';

    legend({'Reservoir Data','Quartz SS','Calcite','Dolomite'}, ...
           'TextColor','w','Location','southwest')
end



%% ==========================================================
%  NEUTRON–SONIC CROSSPLOT (RESERVOIR ZONE ONLY) - FIXED
% ==========================================================

idxNS = flag_reservoir & ~isnan(NPHI) & ~isnan(DT);

fprintf('\nJumlah data Neutron–Sonic reservoir: %d\n', sum(idxNS));

if sum(idxNS) < 10
    warning('Data reservoir terlalu sedikit untuk Neutron–Sonic crossplot');
else

    % --- DATA ---
    phiN = NPHI(idxNS) * 100;   % %
    DTv  = DT(idxNS);           % us/ft
    DepthNS = DEPT(idxNS);

    % --- PARAMETERS (Figure 4.12) ---
    dt_fl = 189;        % us/ft (fresh water)
    dt_qz = 55.5;       % quartz
    dt_ls = 47.5;       % limestone
    dt_dl = 43.5;       % dolomite

    phi = linspace(0,45,200);   % %

    % --- Wyllie Time-Average Curves ---
    DT_qz = dt_qz + phi/100 .* (dt_fl - dt_qz);
    DT_ls = dt_ls + phi/100 .* (dt_fl - dt_ls);
    DT_dl = dt_dl + phi/100 .* (dt_fl - dt_dl);

    %% ---- PLOT ----
    figure('Color','k','Position',[100 100 800 650]); hold on;

    scatter(phiN, DTv, 25, DepthNS, 'filled', 'MarkerFaceAlpha',0.75);
    colormap(jet); cb = colorbar;
    cb.Label.String = 'Depth (ft)';
    cb.Color = 'w';

    plot(phi, DT_qz,'r','LineWidth',2);
    plot(phi, DT_ls,'b','LineWidth',2);
    plot(phi, DT_dl,'g','LineWidth',2);

    xlabel('\phi_{N,LS} (%)','Color','w');
    ylabel('\Delta t (us/ft)','Color','w');
    title('Neutron–Sonic Crossplot (Reservoir Zone Only)','Color','w');

    xlim([0 45]);
    ylim([40 120]);      % <<< PENTING: TIDAK DIBALIK
    set(gca,'Color','k','XColor','w','YColor','w');
    grid on;

    legend({'Reservoir Data','Quartz','Limestone','Dolomite'}, ...
           'TextColor','w','Location','best');

end



%% ==========================================================
%  LITHOLOGY DESCRIPTION FROM NEUTRON–SONIC (COMMAND WINDOW)
% ==========================================================

fprintf('\n====================================================\n');
fprintf(' LITHOLOGY INTERPRETATION (Neutron–Sonic Crossplot)\n');
fprintf(' Reservoir Zones Only\n');
fprintf('====================================================\n');

% --- Pastikan pakai data yg sama dgn crossplot ---
idxNS = flag_reservoir & ~isnan(NPHI) & ~isnan(DT);

phiN = NPHI(idxNS) * 100;   % % Limestone scale
DTv  = DT(idxNS);           % us/ft
DepthNS = DEPT(idxNS);

% --- Lithology reference lines ---
phi_ref = phiN;   % pakai phi yg sama utk jarak

DT_qz_ref = dt_qz + phi_ref/100 .* (dt_fl - dt_qz);
DT_ls_ref = dt_ls + phi_ref/100 .* (dt_fl - dt_ls);
DT_dl_ref = dt_dl + phi_ref/100 .* (dt_fl - dt_dl);

% --- Hitung jarak ke tiap garis ---
d_qz = abs(DTv - DT_qz_ref);
d_ls = abs(DTv - DT_ls_ref);
d_dl = abs(DTv - DT_dl_ref);

% --- Tentukan litologi dominan ---
Lith_NS = strings(length(DTv),1);

for i = 1:length(DTv)
    [~, idx_min] = min([d_qz(i), d_ls(i), d_dl(i)]);

    switch idx_min
        case 1
            Lith_NS(i) = "Quartz Sandstone";
        case 2
            Lith_NS(i) = "Limestone (Calcite)";
        case 3
            Lith_NS(i) = "Dolomite";
    end
end

%% ==========================================================
%  GROUP INTO ZONES (TOP–BASE)
% ==========================================================

dz = median(diff(DEPT));
flag = idxNS;

idx_zone = find(diff([0; flag; 0]) ~= 0);

zone_id = 1;

for i = 1:2:length(idx_zone)

    z_top  = DEPT(idx_zone(i));
    z_base = DEPT(idx_zone(i+1)-1);

    idx_z = DepthNS >= z_top & DepthNS <= z_base;

    lith_unique = unique(Lith_NS(idx_z));
    lith_unique = lith_unique(~ismissing(lith_unique));

    fprintf('\nZone %d\n', zone_id);
    fprintf('Depth Interval : %.1f – %.1f ft\n', z_top, z_base);

    if numel(lith_unique) == 1
        fprintf('Dominant Lithology : %s\n', lith_unique(1));
    else
        fprintf('Dominant Lithology : Mixed (');
        fprintf('%s', lith_unique(1));
        for j = 2:numel(lith_unique)
            fprintf(', %s', lith_unique(j));
        end
        fprintf(')\n');
    end

    fprintf('Interpretation   : Based on Neutron–Sonic proximity\n');

    zone_id = zone_id + 1;
end

fprintf('\n====================================================\n');
fprintf(' END OF LITHOLOGY DESCRIPTION\n');
fprintf('====================================================\n');

%% ==========================================================
%  LITHOLOGY INTERPRETATION FROM N–D CROSSPLOT
%  (Reservoir Zone Only – Command Window Output)
% ==========================================================

idxND = flag_reservoir & ~isnan(NPHI) & ~isnan(RHOB);

if sum(idxND)==0
    disp('Tidak ada data reservoir untuk interpretasi litologi N–D');
else

    phiN = NPHI(idxND) * 100;   % φN,LS (%)
    rhoB = RHOB(idxND);
    depthND = DEPT(idxND);

    % Mineral properties
    rho_fl = 1.00;
    rho_qz = 2.65;
    rho_ls = 2.71;
    rho_dl = 2.87;

    LithND = strings(size(phiN));

    for i = 1:length(phiN)

        % Hitung RHOB mineral pada porositas titik tsb
        rho_qz_p = (1 - phiN(i)/100)*rho_qz + (phiN(i)/100)*rho_fl;
        rho_ls_p = (1 - phiN(i)/100)*rho_ls + (phiN(i)/100)*rho_fl;
        rho_dl_p = (1 - phiN(i)/100)*rho_dl + (phiN(i)/100)*rho_fl;

        % Jarak ke tiap garis mineral
        d_qz = abs(rhoB(i) - rho_qz_p);
        d_ls = abs(rhoB(i) - rho_ls_p);
        d_dl = abs(rhoB(i) - rho_dl_p);

        [~,idxMin] = min([d_qz d_ls d_dl]);

        switch idxMin
            case 1
                LithND(i) = "Sandstone (Quartz)";
            case 2
                LithND(i) = "Limestone (Calcite)";
            case 3
                LithND(i) = "Dolomite";
        end
    end

    %% ---- RINGKASAN PER ZONA RESERVOIR ----
    fprintf('\n===== LITHOLOGY INTERPRETATION (N–D CROSSPLOT) =====\n');

    % Kelompokkan per interval reservoir
    flag = flag_reservoir & ~isnan(DEPT);
    idx = find(diff([0; flag; 0])~=0);

    for k = 1:2:length(idx)

        top = DEPT(idx(k));
        base = DEPT(idx(k+1)-1);

        zoneIdx = depthND >= top & depthND <= base;

        if sum(zoneIdx)==0, continue; end

        lithZone = LithND(zoneIdx);
        lithMode = mode(categorical(lithZone));

        fprintf('Depth %.1f – %.1f ft  →  Dominant Lithology: %s\n', ...
                top, base, string(lithMode));
    end

    fprintf('===================================================\n');
end

%% ==========================================================
%  CALCULATE M–N PARAMETERS (RESERVOIR ZONE ONLY)
% ==========================================================

idxMN = flag_reservoir & ~isnan(DT) & ~isnan(RHOB) & ~isnan(NPHI);

DepthMN = DEPT(idxMN);

% Konstanta fluida & matriks (Quartz)
DT_fl   = 189;     
DT_ma   = 55.5;    
RHOB_fl = 1.00;    
RHOB_ma = 2.65;    

M = (DT(idxMN) - DT_fl) ./ (DT_ma - DT_fl);
N = (RHOB(idxMN) - RHOB_fl) ./ (RHOB_ma - RHOB_fl);

%% ==========================================================
%  DEFINE RESERVOIR ZONES FOR M–N
% ==========================================================

flag = flag_reservoir & ~isnan(DEPT);
idxZ = find(diff([0; flag; 0]) ~= 0);

nZone = length(idxZ)/2;

ZoneTop  = zeros(nZone,1);
ZoneBase = zeros(nZone,1);

for i = 1:nZone
    ZoneTop(i)  = DEPT(idxZ(2*i-1));
    ZoneBase(i) = DEPT(idxZ(2*i)-1);
end

zoneColor = lines(nZone);

%% ==========================================================
%  PLOT M–N CROSSPLOT (COLORED BY ZONE)
% ==========================================================

figure('Color','w','Position',[100 100 750 650]);
hold on; grid on; box on;

legendText = {};

for z = 1:nZone
    idxZone = DepthMN >= ZoneTop(z) & DepthMN <= ZoneBase(z);
    scatter(N(idxZone), M(idxZone), 35, ...
            zoneColor(z,:), 'filled');
    legendText{end+1} = sprintf('Zone %d (%.0f–%.0f ft)', ...
                                 z, ZoneTop(z), ZoneBase(z));
end

xlabel('N','FontSize',12);
ylabel('M','FontSize',12);
title('M–N Lithology Crossplot (Reservoir Zones Only)','FontSize',14);

xlim([0.3 0.8])
ylim([0.5 1.1])

legend(legendText,'Location','bestoutside');

%% ==========================================================
%  MINERAL REFERENCE POINTS
% ==========================================================

plot(0.72,0.83,'ks','MarkerFaceColor','y'); text(0.725,0.83,' Quartz SS');
plot(0.63,0.84,'ks','MarkerFaceColor','c'); text(0.635,0.84,' Limestone');
plot(0.52,0.78,'ks','MarkerFaceColor','g'); text(0.525,0.78,' Dolomite');
plot(0.55,0.70,'ks','MarkerFaceColor','m'); text(0.555,0.70,' Anhydrite');
plot(0.50,0.62,'ks','MarkerFaceColor','r'); text(0.505,0.62,' Shale');

%% ==========================================================
%  INTERPRETATION GUIDELINES (OPTIONAL)
% ==========================================================

text(0.60,0.92,'Secondary Porosity','FontSize',10,'Color',[0 0 0.6]);
text(0.70,0.90,'Gas / Salt Effect','FontSize',10,'Color',[0.6 0 0]);

%% ==========================================================
%  LITHOLOGY INTERPRETATION PER RESERVOIR ZONE (M–N)
% ==========================================================

% Reference minerals (N, M)
MineralName = {'Quartz Sandstone','Limestone','Dolomite','Anhydrite','Shale'};
MineralRef = [ ...
    0.72 0.83;   % Quartz
    0.63 0.84;   % Limestone
    0.52 0.78;   % Dolomite
    0.55 0.70;   % Anhydrite
    0.50 0.62];  % Shale

fprintf('\n===== M–N LITHOLOGY INTERPRETATION (RESERVOIR ZONES) =====\n');

for z = 1:nZone

    idxZone = DepthMN >= ZoneTop(z) & DepthMN <= ZoneBase(z);

    Mz = median(M(idxZone));
    Nz = median(N(idxZone));

    % Distance to reference minerals
    dist = sqrt((MineralRef(:,1)-Nz).^2 + (MineralRef(:,2)-Mz).^2);
    [~,id] = min(dist);

    fprintf('Zone %d (%.0f–%.0f ft): Dominant Lithology = %s\n', ...
        z, ZoneTop(z), ZoneBase(z), MineralName{id});
end


%% ==========================================================
%  MATRIX IDENTIFICATION PLOT
%  Using Neutron, Density, and Sonic
%  Reservoir Zone Only
% ==========================================================

% ---- Filter zona reservoir ----
idxMat = flag_reservoir & ~isnan(RHOB) & ~isnan(NPHI) & ~isnan(DT);

if sum(idxMat) == 0
    warning('Tidak ada data reservoir untuk Matrix Identification Plot');
else

    % ---- Log data ----
    RHOBv = RHOB(idxMat);        % g/cc
    NPHIv = NPHI(idxMat);        % v/v (decimal)
    DTv   = DT(idxMat);          % us/ft
    DEPv  = DEPT(idxMat);

    % ---- Fluid properties (sesuai gambar) ----
    rho_f = 1.00;       % g/cc (fresh water)
    dt_f  = 189;        % us/ft (freshwater mud)

    % ---- Apparent matrix properties ----
    rho_ma = (RHOBv - NPHIv .* rho_f) ./ (1 - NPHIv);
    dt_ma  = (DTv   - NPHIv .* dt_f ) ./ (1 - NPHIv);

    %% ---- Reference matrix points (Fig 4.18) ----
    %            Quartz   Calcite   Dolomite
    rho_ref = [2.65       2.71       2.87];
    dt_ref  = [55.5       47.6       43.5];

    %% ---- Plot ----
    figure('Color','k','Position',[100 100 750 700]);
    hold on; grid on;

    scatter(dt_ma, rho_ma, 25, DEPv, 'filled', ...
        'MarkerFaceAlpha',0.75);

    colormap(jet)
    cb = colorbar;
    cb.Label.String = 'Depth (ft)';
    cb.Color = 'w';

    % ---- Plot matrix vertices ----
    plot(dt_ref(1), rho_ref(1),'ro','MarkerSize',9,'LineWidth',2)
    plot(dt_ref(2), rho_ref(2),'bo','MarkerSize',9,'LineWidth',2)
    plot(dt_ref(3), rho_ref(3),'go','MarkerSize',9,'LineWidth',2)

    text(dt_ref(1)+0.6, rho_ref(1),'Quartz','Color','r','FontWeight','bold')
    text(dt_ref(2)+0.6, rho_ref(2),'Calcite','Color','b','FontWeight','bold')
    text(dt_ref(3)+0.6, rho_ref(3),'Dolomite','Color','g','FontWeight','bold')

    % ---- Triangle (matrix domain) ----
    plot([dt_ref dt_ref(1)], [rho_ref rho_ref(1)], 'w--', 'LineWidth',1.2)

    % ---- Axis & style ----
    set(gca,'YDir','reverse','Color','k','XColor','w','YColor','w')
    xlabel('\Delta t_{ma}  (µs/ft)','Color','w')
    ylabel('\rho_{ma}  (g/cc)','Color','w')
    title('Matrix Identification Plot (Reservoir Zone Only)','Color','w')

    xlim([35 75])
    ylim([2.4 3.05])

end


%% ==========================================================
%  ZONA LITHOLOGY INTERPRETATION
%  Based on Matrix Identification Plot (N-D-S)
% ==========================================================

fprintf('\n====================================================\n');
fprintf(' ZONAL LITHOLOGY INTERPRETATION (Reservoir Only)\n');
fprintf(' Based on Neutron–Density–Sonic Matrix Plot\n');
fprintf('====================================================\n');

ZoneID = 1;

for z = 1:height(ReservoirIntervalTable)

    zTop  = ReservoirIntervalTable.TopDepth(z);
    zBase = ReservoirIntervalTable.BaseDepth(z);

    % data di dalam zona ini
    idxZ = (DEPv >= zTop) & (DEPv <= zBase);

    if sum(idxZ) == 0
        continue
    end

    % hitung jarak ke masing-masing mineral
    dQ = abs(rho_ma(idxZ) - rho_ref(1)) + abs(dt_ma(idxZ) - dt_ref(1));
    dC = abs(rho_ma(idxZ) - rho_ref(2)) + abs(dt_ma(idxZ) - dt_ref(2));
    dD = abs(rho_ma(idxZ) - rho_ref(3)) + abs(dt_ma(idxZ) - dt_ref(3));

    % klasifikasi tiap titik
    lith_idx = zeros(sum(idxZ),1);

    for i = 1:length(dQ)
        [~, ii] = min([dQ(i), dC(i), dD(i)]);
        lith_idx(i) = ii;
    end

    % litologi dominan
    dom = mode(lith_idx);

    switch dom
        case 1, domLith = 'Quartz Sandstone';
        case 2, domLith = 'Limestone';
        case 3, domLith = 'Dolomite';
    end

    fprintf('Zone %d (%.1f–%.1f ft): Dominant Lithology = %s\n', ...
        ZoneID, zTop, zBase, domLith);

    ZoneID = ZoneID + 1;
end


%% ==========================================================
%  NET RESERVOIR & NET PAY FLAG
% ==========================================================

Vsh_cut = 0.50;
phi_cut = 0.15;
Sw_cut  = 0.70;
k_cut   = 1;     % mD (opsional)

flag_net_res = (Vsh <= Vsh_cut) & (phi_eff >= phi_cut);
flag_net_pay = flag_net_res & (Sw <= Sw_cut) & (k >= k_cut);


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

    % SILTSTONE (TRANSITION)
    else
        LithCode(i) = 4;
    end
end


%% ==========================================================
% 4. VISUALIZATION – QUICKLOOK & PETROPHYSICS (FIXED)
% ==========================================================

figure('Color','k','Position',[50 50 2000 900]);
sgtitle('TUGAS AKHIR GABE','Color','w','FontSize',18);

%% ---- TRACK 1: LITHOLOGY ----
ax0 = subplot(1,9,1); hold on

for i = 1:length(DEPT)-1
    if isnan(LithCode(i)), continue; end

    z1 = DEPT(i); z2 = DEPT(i+1);

    switch LithCode(i)
        case 1, col = [0 0.7 0];        % Shale
        case 2, col = [1 1 0];          % Sand
        case 3, col = [0 0.4 1];        % Limestone
        case 4, col = [0.6 0.4 0.2];    % Siltstone
    end

    patch([0 1 1 0],[z1 z1 z2 z2],col,'EdgeColor','none');
end

set(gca,'YDir','reverse','Color','k','XColor','w','YColor','w')
xlim([0 1]); xticks([])
ylabel('Depth (ft)','Color','w')
title('Lithology','Color','w'); grid on


%% ---- TRACK 2: GR ----
ax1 = subplot(1,9,2);
hold on

% ===== DEFINISI IDX (WAJIB ADA) =====
flag = flag_GR & ~isnan(GR);
idx  = find(diff([0; flag; 0]) ~= 0);

% ===== SHADING ZONA PERMEABEL =====
for i = 1:2:length(idx)
    z_top  = DEPT(idx(i));
    z_base = DEPT(idx(i+1)-1);

    patch([0 GR_cutoff GR_cutoff 0], ...
          [z_top z_top z_base z_base], ...
          [1 1 0], ...
          'FaceAlpha',0.25, ...
          'EdgeColor','none');
end

% ===== PLOT GR DI ATAS SHADING =====
plot(GR, DEPT,'g','LineWidth',1.3)

% ===== GR CUTOFF =====
xline(GR_cutoff,'r--','GR cutoff','LineWidth',1.2);

set(gca,'YDir','reverse','Color','k','XColor','w','YColor','w')
xlabel('GR (API)','Color','w')
title('Gamma Ray','Color','w')
xlim([0 150])
grid on


%% ---- TRACK 3: RESISTIVITY ----
ax2 = subplot(1,9,3);
semilogx(LLD,DEPT,'r','LineWidth',1); hold on
semilogx(LLS,DEPT,'m','LineWidth',1);
set(gca,'YDir','reverse','Color','k','XColor','w','YColor','w')
xlabel('Rt (ohm.m)','Color','w')
title('Resistivity','Color','w')
xlim([0.2 2000]); grid on
legend({'LLD','LLS'},'TextColor','w','Location','best')

%% ---- TRACK 4: DENSITY–NEUTRON ----
ax3 = subplot(1,9,4);
plot(RHOB,DEPT,'r','LineWidth',1); hold on
xlim([1.7 2.7])

% NPHI (dibalik secara visual)
NPHI_scaled = 1.7 + (0.6 - NPHI)*(1.0/0.6);
plot(NPHI_scaled,DEPT,'c','LineWidth',1);
set(gca,'YDir','reverse','Color','k','XColor','w','YColor','w')
xlabel('RHOB / NPHI','Color','w')
title('Density–Neutron','Color','w')
xlim([1.7 2.7]); grid on
legend({'RHOB','NPHI'},'TextColor','w','Location','best')

%% ---- TRACK 5: Vshale ----
ax4 = subplot(1,9,5);
plot(Vsh,DEPT,'y','LineWidth',1);
set(gca,'YDir','reverse','Color','k','XColor','w','YColor','w')
xlabel('Vsh','Color','w')
title('Vshale','Color','w')
xlim([0 1]); grid on

% optional cutoff
xline(0.35,'r--','Vsh cutoff');
%% ---- TRACK 5: EFFECTIVE POROSITY ----
ax5 = subplot(1,9,6);
plot(phi_eff,DEPT,'c','LineWidth',1);
set(gca,'YDir','reverse','Color','k','XColor','w','YColor','w')
xlabel('\phi_{eff}','Color','w')
title('Effective Porosity','Color','w')
xlim([0 0.4]); grid on
%% ---- TRACK 6: WATER SATURATION ----
ax6 = subplot(1,9,7);
plot(Sw,DEPT,'y','LineWidth',1);
set(gca,'YDir','reverse','Color','k','XColor','w','YColor','w')
xlabel('Sw','Color','w')
title('Water Saturation','Color','w')
xlim([0 1]); grid on
%% ---- TRACK 7: PERMEABILITY ----
ax7 = subplot(1,9,8);

k_plot = k;
k_plot(k_plot <= 0) = NaN;

semilogx(k_plot, DEPT,'w','LineWidth',1); hold on
xline([0.1 1 10 100 1000],'w:');

set(gca,'YDir','reverse','Color','k','XColor','w','YColor','w')
xlabel('k (mD)','Color','w')
title('Permeability','Color','w')
xlim([1e-3 1e4])
grid on

%% ---- TRACK 8: NET RESERVOIR & NET PAY ----
axNR = subplot(1,9,9);   % sesuaikan jika layout Anda 1x9
hold on

for i = 1:length(DEPT)-1

    z1 = DEPT(i);
    z2 = DEPT(i+1);

    % NET RESERVOIR (biru muda)
    if flag_net_res(i)
        patch([0 1 1 0], [z1 z1 z2 z2], ...
              [0.3 0.6 1], 'FaceAlpha',0.30, 'EdgeColor','none');
    end

    % NET PAY (merah)
    if flag_net_pay(i)
        patch([0 1 1 0], [z1 z1 z2 z2], ...
              [1 0 0], 'FaceAlpha',0.45, 'EdgeColor','none');
    end
end

set(gca,'YDir','reverse','Color','k','XColor','w','YColor','w')
xlim([0 1])
xticks([])
title('Net Res / Net Pay','Color','w')
grid on


%% ---- LINK DEPTH AXIS ----
linkaxes([ax0 ax1 ax2 ax3 ax4 ax5 ax6 ax7 axNR],'y');


%% ==========================================================
%  NET THICKNESS & NTG (TARGET INTERVAL ONLY)
% ==========================================================

% Vertical sampling
dz = median(diff(DEPT));   % ft

% Target interval thickness (quicklook reservoir zone)
target_thickness = sum(flag_reservoir) * dz;

% Net Reservoir Thickness (di dalam target)
net_res_thickness = sum(flag_net_res & flag_reservoir) * dz;

% Net Pay Thickness (di dalam target)
net_pay_thickness = sum(flag_net_pay & flag_reservoir) * dz;

% NTG (target-based)
NTG_target = (net_res_thickness / target_thickness) ;

%% =========================================================
%                       DISPLAY RESULT
% ==========================================================

fprintf('\n===== NET THICKNESS & NTG (TARGET INTERVAL) =====\n');
fprintf('Target Thickness     : %.2f ft\n', target_thickness);
fprintf('Net Reservoir        : %.2f ft\n', net_res_thickness);
fprintf('Net Pay              : %.2f ft\n', net_pay_thickness);
fprintf('NTG (Target-based)   : %.3f\n', NTG_target);


%% ==========================================================
%  RESERVOIR TARGET INTERVAL (TOP–BASE)
% ==========================================================

dz = median(diff(DEPT));   % ft

flag = flag_reservoir & ~isnan(DEPT);

% Cari perubahan zona (0→1 atau 1→0)
idx = find(diff([0; flag; 0]) ~= 0);

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

% Buat tabel hasil
ReservoirIntervalTable = table(TopDepth, BaseDepth, Thickness);

%% =========================
% DISPLAY RESULT
% =========================
disp('===== RESERVOIR TARGET INTERVAL =====');
disp(ReservoirIntervalTable);



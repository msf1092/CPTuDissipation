clear; close all; format compact; clc
%% inputs

wat_table = input('Enter the water table level: \n'); 
chartlinewidth = 1.5;

%% for FINNINMÄKI:

load("data.mat")
load("data_cptu.mat")

a = fieldnames(data);

depthCPT = data_cptu.u(:,1);
qc = data_cptu.qc(:,2);
fs = data_cptu.fs(:,2);
u2 = data_cptu.u(:,2);
Ic = data_cptu.Ic(:,2);

% Plot dissipation tests and CPTu results sections in one plot
figure('Name','Individual dissipation tests', 'units','centimeters', 'Position', [0.1 0.6 18 numel(a)*4]) 
set(gca,'defaulttextinterpreter','latex');
tiledlayout(3*(numel(a)-1),7, 'Padding', 'compact', 'TileSpacing', 'compact'); 
for i = 1 : numel(a)-1
    b = char(a(i));
    leg{i} = b;
    b = data.(b);
    u = b(:,2);
    t = b(:,1);

    drem = leg{i};
    drem = strrep(drem, 'd','');
    ddiss = vpa(drem)/100; % depth where dissipation is performed!
    SectRange = 5; % the indices range for sections plots
    LocalParams{i} = fncPlotCPTSectDiss(depthCPT,ddiss,qc,fs,u2,Ic,SectRange, u,t, leg{i}, i);
end

%% Compressive, single plot
testNo = 5; % Specify the test no here

fig_name = "Compressive";
f = figure('Name',fig_name, 'Units','centimeters', 'Position',[2 2 9 6.5]);
set(f,'defaulttextinterpreter','latex');

linecolors = ["#0072BD"; "#0072BD"; "#0072BD"; "#0072BD"; "#77AC30"; "#77AC30"; ...
    "#77AC30"; "#77AC30"; "#00FF00"; "#00FF00"; "#00FF00"; "#00FF00"; "#000000"; "#000000"];
linetypes = ["-"; "--"; ":"; "-."; "-"; "--"; ":"; "-."];

b = char(a(testNo));
leg = b;
b = data.(b);
u = b(:,2);
t = b(:,1);

drem = strrep(leg, 'd','');
ddiss = vpa(drem)/100;
y0 = (double(ddiss) - wat_table) * 9.81; % the water pressure at rest

[m, ind_m] = max(u);
n = min(u(ind_m:end))-0.1;
norm_u = (u-n) / (m-n);

tt = (min(t):(max(t)-min(t))/999:max(t))';
p = 1.00000;
uu = csaps(t,u,p,tt);

% Computation of t50 , u50 (t from the maxi of u2)
u50 = n + (m-n)/2;
[~, idu] = min(abs(uu - u50));
t50 = tt(idu);

str = num2str(linecolors(1));
color = sscanf(str(2:end), '%2x%2x%2x', [1 3]) / 255;
semilogx(t,u, "LineStyle", linetypes(1), "LineWidth",chartlinewidth,"Color",'k','DisplayName', sprintf('$d= %0.2f \\hspace{0.1cm} m$', ddiss));
hold on
semilogx([1 max(t)],[y0 y0], 'LineStyle', '-.', 'Color','red', 'DisplayName', 'in-situ $u$', 'HandleVisibility', 'off')
scatter(t50,u50,50,"Marker","o","MarkerEdgeColor","blue", 'DisplayName','t50','HandleVisibility','off')
hold off
xLab = "Time [sec]";
yLab = "Pore water pressure, $u_2$ [kPa]";
xlabel(xLab,'FontSize',10,'Color','k','Interpreter','latex')
ylabel(yLab,'FontSize',10,'Color','k','Interpreter','latex')
xlim([0 10^3])
ylim([0 100])
q = char(leg);
legend('FontSize',8, 'Location','northeast', 'Interpreter','latex', "color","white");
ax = gca;
set(ax,'TickLabelInterpreter','latex')
text(t50,u50,sprintf('\\hspace{0.2cm} $t_{50}$ = %-5.1f s',t50), 'FontSize',8,'VerticalAlignment','middle','HorizontalAlignment','left','Interpreter','latex')
grid on
grid minor
box on

%% Sully et al., 1999:
testNo = 3; % Specify the test no here

b = char(a(testNo));
leg = b;
b = data.(b);
u = b(:,2);
t = b(:,1);
t_root = t .^ 0.5;

x = t_root;
y = u;

% Initial data
test_depth = vpa(strrep(leg,'d',''))/100;
y0 = (test_depth - wat_table) * 9.81; 

% Initial line ranges
x1_min = 0; x1_max = max(x); 
y1_min = y0; y1_max = y0;
x2_min = 0; x2_max = (1 / 2) * max(x); 
y2_min = 2 * max(y); y2_max = y0;

[u0_sqrt, u_max, t_50, u_50] = fnc_gui_with_sliders_Sully(x,y,y0,x1_min,x1_max,y1_min,y1_max,x2_min,x2_max,y2_min,y2_max,test_depth);

%% Chai et al., 2012:
% Method Three
b = char(a(testNo)); % specify the test 
leg = b;
b = data.(b);
u = b(:,2);
t = b(:,1);
lege = vpa(strrep(leg,'d',''))/100;

qt = data_cptu.qt(:,2);
Bq = data_cptu.Bq(:,2);

% Find the peak of the curve
[peakValue, peakIdx] = findpeaks(u, "NPeaks", 1);
y0 = (double(lege) - wat_table) * 9.81; % the water pressure at rest

t_cut = t(peakIdx:end); % Time values starting from the peak
t_cut = t_cut - t_cut(1); % Time starting from zero, from the peak
u_cut = u(peakIdx:end); % u values starting from the peak
u_norm = 2 * (u_cut - y0) / (max(u_cut) - y0);

% Remove duplicate u_norm values
[u_norm_unique, uniqueIdx] = unique(u_norm);
t_cut_unique = t_cut(uniqueIdx);

% Find t at u = 1 using interpolation
if min(u_norm_unique) <= 1 && max(u_norm_unique) >= 1
    t_at_u_1 = interp1(u_norm_unique, t_cut_unique, 1, 'linear'); % Linear interpolation
    fprintf('The value of t at u = 1 is: %.4f\n', t_at_u_1);
else
    t_at_u_1 = NaN; % If u=1 is outside the range
    fprintf('u = 1 is outside the range of u_norm. Interpolation not possible.\n');
end

Nkt_Fu = exp(6.41424 - 1.489 * Ic(depthCPT(:, 1) == double(lege), 1));
gamma = 11.46 + 0.33 * log10(depthCPT) + 3.10 * log10(fs*1000) + 0.70 * log10(qt*1000); % After Mayne 2010
[,ind_depth] = find(depthCPT == double(lege));
sig_v = (depthCPT(1,1)*17 + sum(diff(depthCPT(1:ind_depth,1)) .* gamma(1:ind_depth-1,1))) / 1000; % in MPa
sig_prime_v = (depthCPT(1,1)*17 + sum(diff(depthCPT(1:ind_depth,1)) .* gamma(1:ind_depth-1,1))  -  (double(lege) - wat_table) * 9.81) / 1000; % in MPa
Su_kt_Fu = (qt(depthCPT(:, 1) == double(lege), 1) - sig_v) / Nkt_Fu; % in MPa
% eq. 11:
Nke_Fu = exp(6.31207 - 1.435 * Ic(depthCPT(:, 1) == double(lege), 1));
Su_ke_Fu = (qt(depthCPT(:, 1) == double(lege), 1) - u2(depthCPT(:, 1) == double(lege), 1)) / Nke_Fu; % in MPa
% Su after Robertson & Cabal, 7th ed.
N_kt_avg = 14;
Su_RobCab = (qt(depthCPT(:, 1) == double(lege), 1) - sig_v) / N_kt_avg;


% SCPTu
depths_G = [1.5; 1.74; 2; 2.52; 3; 3.5; 3.84; 4.52; 5; 5.5];
G = [0; 228.610524823260; 119.406308277414; 78.17968445; 102.3092638; 85.30360452; 64.17454772; 92.88935609; 47.75957637; 46.77538535];

% plot G
figure('Name','Seismic tests-SCPTu', 'units','centimeters', 'Position', [2 2 5 10]) 
set(gca,'defaulttextinterpreter','latex');
hold on;
grid on; grid minor

% Loop through each segment and draw a horizontal line
for i = 1:(length(depths_G) - 1)
    plot([G(i+1), G(i+1)], [depths_G(i), depths_G(i+1)], 'LineWidth', chartlinewidth);
end
xlabel('$G$ (MPa)', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('Depth (m)', 'Interpreter', 'latex', 'FontSize', 10);
xlim([0 350])
hold off;
set(gca,'TickLabelInterpreter','latex')
set (gca, 'YDir','reverse')
box on


for i = 2:length(depths_G) - 1
    if double(lege) < depths_G(i) && double(lege) > depths_G(i-1)
        G_test = G(i,1);
        Ir_kt_Fu = G_test / Su_kt_Fu;
        Ir_ke_Fu = G_test / Su_ke_Fu;
        Ir_kt_RobCab = G_test / Su_RobCab;
    end
end
t_50c_Fu_kt = t_at_u_1 / (1 + 0.001*18.5 * (t(peakIdx)/t_at_u_1)^0.57 * (Ir_kt_Fu/200)^0.3 ); 
t_50c_Fu_ke = t_at_u_1 / (1 + 0.001*18.5 * (t(peakIdx)/t_at_u_1)^0.57 * (Ir_ke_Fu/200)^0.3 );
t_50c_RobCab = t_at_u_1 / (1 + 0.001*18.5 * (t(peakIdx)/t_at_u_1)^0.57 * (Ir_kt_RobCab/200)^0.3 );


% Plot the original and cut/normalized chart
figure('Name', 'Second method', 'units', 'centimeters', 'Position', [3 2 9 7]);
set(gca, 'defaulttextinterpreter', 'latex');
semilogx(t_cut, u_norm, 'Color','k', 'LineWidth', chartlinewidth, 'DisplayName', sprintf('$d= %0.2f \\hspace{0.1cm} m$', lege));
xlabel('Time (sec)', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('Normalized pore pressure ($-$)', 'Interpreter', 'latex', 'FontSize', 10);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 10, 'Location','best');
hold on;
if ~isnan(t_at_u_1)
    plot(t_at_u_1, 1, 'ro', 'MarkerSize', 8, 'DisplayName', '$t_{50}$', 'HandleVisibility', 'on');
    plot(t_50c_Fu_kt, 1.2, 'r^', 'MarkerSize', 8, 'DisplayName', '$s_u$, Fu, $k_t$', 'HandleVisibility', 'on');
    plot(t_50c_Fu_ke, 1.4, 'rx', 'MarkerSize', 8, 'DisplayName', '$s_u$, Fu, $k_e$', 'HandleVisibility', 'on');
    plot(t_50c_RobCab, 1.6, 'r+', 'MarkerSize', 8, 'DisplayName', '$s_u$, R\&C, $k_{avg}$', 'HandleVisibility', 'on');
end
hold off;
text(t_at_u_1, 1, sprintf('\\hspace{0.3cm}$t_{50}$=%0.2f sec',t_at_u_1), 'Interpreter', 'latex', 'FontSize', 9, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Color', 'blue');
text(t_50c_Fu_kt, 1.2, sprintf('\\hspace{0.0cm}$t_{50c}$=%0.2f sec\\hspace{0.2cm}',t_50c_Fu_kt), 'Interpreter', 'latex', 'FontSize', 9, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Color', 'blue');
text(t_50c_Fu_ke, 1.4, sprintf('\\hspace{0.0cm}$t_{50c}$=%0.2f sec\\hspace{0.2cm}',t_50c_Fu_ke), 'Interpreter', 'latex', 'FontSize', 9, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Color', 'blue');
text(t_50c_RobCab, 1.6, sprintf('\\hspace{0.0cm}$t_{50c}$=%0.2f sec\\hspace{0.2cm}',t_50c_RobCab), 'Interpreter', 'latex', 'FontSize', 9, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Color', 'blue');
set(gca, 'TickLabelInterpreter', 'latex');
% Grid and box
grid on; grid minor;
box on;

%% Method Two:
b = char(a(testNo)); % specify the test 
leg = b;
b = data.(b);
u = b(:,2);
t = b(:,1);
lege = vpa(strrep(leg,'d',''))/100;

% Find the peak of the curve
[peakValue, peakIdx] = findpeaks(u, "NPeaks", 1);
y0 = (double(lege) - wat_table) * 9.81;

t_cut = t(peakIdx:end); % Time values starting from the peak
u_cut = u(peakIdx:end); % u values starting from the peak
t_cut = t_cut - t_cut(1);
u_norm = 2 * (u_cut - y0) / (max(u_cut) - y0);

% Remove duplicate u_norm values
[u_norm_unique, uniqueIdx] = unique(u_norm);
t_cut_unique = t_cut(uniqueIdx);

% Find t at u = 1 using interpolation
if min(u_norm_unique) <= 1 && max(u_norm_unique) >= 1
    t_at_u_1 = interp1(u_norm_unique, t_cut_unique, 1, 'linear'); % Linear interpolation
    fprintf('The value of t at u = 1 is: %.4f\n', t_at_u_1);
else
    t_at_u_1 = NaN; % If u=1 is outside the range
    fprintf('u = 1 is outside the range of u_norm. Interpolation not possible.\n');
end

% Plot the original and cut/normalized chart
figure('Name', 'Third method', 'units', 'centimeters', 'Position', [3 2 9 6.5]);
set(gca, 'defaulttextinterpreter', 'latex');
semilogx(t_cut, u_norm, 'LineWidth', chartlinewidth, 'Color','k', 'DisplayName', sprintf('$d= %0.2f \\hspace{0.1cm} m$', lege));
xlabel('Time (sec)', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('Normalized pore pressure ($-$)', 'Interpreter', 'latex', 'FontSize', 10);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 10);
hold on;
if ~isnan(t_at_u_1)
    plot(t_at_u_1, 1, 'ro', 'MarkerSize', 8, 'DisplayName', '$u = 1$', 'HandleVisibility', 'off');
end
hold off;
text(t_at_u_1, 1, sprintf('\\hspace{0.3cm}$t_{50}$=%0.2f sec',t_at_u_1), 'Interpreter', 'latex', 'FontSize', 9, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Color', 'blue');
ylim([0 2])
set(gca,'TickLabelInterpreter','latex')
grid on; grid minor;
box on;


%% k estimation for each of the three methods:
% For dilative behavior
depth = [2.13; 2.70; 3.38; 4.39; 5.38; 6.28];
depth_dil = [2.13; 2.70; 3.38; 6.28];
t50_1_dil = [17.28; 6.14; 31.1; 22.43];
t50_2_dil = [29.51; 9.49; 34.03; 32.92];
t50_3_dil = [1.79; 0.49; 1.65; 1.6]; % last row is a fake one
% For contractive behavior
depth_comp = [4.39; 5.38];
t50_comp = [5.30; 3.50];
% k values from lab
depth_lab = [2.45; 3.35; 4.35];
k_lab = [9.16e-9; 4.70e-10; 2.50e-7];


figure('Name','Seismic tests-SCPTu', 'units','centimeters', 'Position', [20 2 18 9]) %'WindowState','maximized'
set(gca,'defaulttextinterpreter','latex');
% For compressive behavior
k_PF_comp = ((1 ./ (25 * t50_comp)) .^ 1.25) .* 0.01; % Sully et al. (1999)
k_M_etal_comp = ((1 ./ (720 .* t50_comp)) .^ 1.05) .* 0.01; % Sully et al. (1999)
% Parez-Fauriel k (1988) for three methods:
k_PF_1 = ((1 ./ (25 * t50_1_dil)) .^ 1.25) .* 0.01; % Sully et al. (1999)
k_PF_2 = ((1 ./ (25 * t50_2_dil)) .^ 1.25) .* 0.01; % Chai et al. (2012)
k_PF_3 = ((1 ./ (25 * t50_3_dil)) .^ 1.25) .* 0.01; % Chai et al. (2012)
% Ziaie-Moayed et al. (2009) k for three methods:
k_M_etal_1 = ((1 ./ (720 .* t50_1_dil)) .^ 1.05) .* 0.01; % Sully et al. (1999)
k_M_etal_2 = ((1 ./ (720 .* t50_2_dil)) .^ 1.05) .* 0.01; % Chai et al. (2012)
k_M_etal_3 = ((1 ./ (720 .* t50_3_dil)) .^ 1.05) .* 0.01; % Chai et al. (2012)

scatter(k_PF_comp,depth_comp,40,"Marker","^","MarkerEdgeColor","red","HandleVisibility","on","DisplayName","$k_h$ after Parez and Fauriel (1988), for compressive tests")
hold on
scatter(k_M_etal_comp,depth_comp,40,"Marker","v","MarkerEdgeColor","red","HandleVisibility","on","DisplayName","$k_h$ after Ziaie-Moayed et al. (2009), for compressive tests")
scatter(k_PF_1,depth_dil,50,"Marker","x","MarkerEdgeColor","blue","HandleVisibility","on","DisplayName","$k_h$ after Parez and Fauriel (1988), $t_{50}$ from Sully et al. (1999)");
scatter(k_PF_2,depth_dil,50,"Marker","x","MarkerEdgeColor","red","HandleVisibility","on","DisplayName","$k_h$ after Parez and Fauriel (1988), $t_{50}$ from Chai et al. (2012)");
scatter(k_PF_3(1:end-1),depth_dil(1:end-1),50,"Marker","x","MarkerEdgeColor","magenta","LineWidth",1.5,"HandleVisibility","on","DisplayName","$k_h$ after Parez and Fauriel (1988), $t_{50c}$ from Chai et al. (2012)");
scatter(k_M_etal_1,depth_dil,50,"Marker","o","MarkerEdgeColor","blue","HandleVisibility","on","DisplayName","$k_h$ after Ziaie-Moayed et al. (2009), $t_{50}$ from Sully et al. (1999)");
scatter(k_M_etal_2,depth_dil,50,"Marker","o","MarkerEdgeColor","red","HandleVisibility","on","DisplayName","$k_h$ after Ziaie-Moayed et al. (2009), $t_{50}$ from Chai et al. (2012)");
scatter(k_M_etal_3(1:end-1),depth_dil(1:end-1),50,"Marker","o","MarkerEdgeColor","magenta","LineWidth",1.5,"HandleVisibility","on","DisplayName","$k_h$ after Ziaie-Moayed et al. (2009), $t_{50c}$ from Chai et al. (2012)");
scatter(k_lab,depth_lab,200,"Marker",".","MarkerEdgeColor","blue","MarkerFaceColor","blue","HandleVisibility","on","DisplayName","$k_v$ from lab tests")
%plot the k values based on I_c
I_c = Ic;
d = depthCPT;
for i = 1 : length(depth)-1 % the last test is excluded, because the fs and u2 values are not trustable at the last depth (for estimation of k), which is the ending depth of the test, just above the bedrock.
    dd = depth(i);
    ind = find(d == dd);
    int = 0.2;
    d_low = dd - int/2;
    d_high = dd + int/2;
    if ind+3 <= length(I_c)
        Ic_kDepths(i,1) = mean(I_c(ind-3:ind+3,1));
    else
        Ic_kDepths(i,1) = mean(I_c(ind-4:ind,1));
    end
    if Ic_kDepths(i,1) > 3.60 % (2) Organic soils-clay
        k_Rob_low(i,1) = 1e-10;
        k_Rob_high(i,1) = 1e-8;

    elseif Ic_kDepths(i,1) >= 2.95 && Ic_kDepths(i,1) < 3.60 % (3) Clay
        k_Rob_low(i,1) = 1e-10;
        k_Rob_high(i,1) = 1e-9;

    elseif Ic_kDepths(i,1) >= 2.60 && Ic_kDepths(i,1) < 2.95 % (4) Silt mixture
        k_Rob_low(i,1) = 3e-9;
        k_Rob_high(i,1) = 1e-7;

    elseif Ic_kDepths(i,1) >= 2.05 && Ic_kDepths(i,1) < 2.60 % (5) Sand mixture
        k_Rob_low(i,1) = 1e-7;
        k_Rob_high(i,1) = 1e-5;

    elseif Ic_kDepths(i,1) >= 1.31 && Ic_kDepths(i,1) < 2.05 % (6) Sand
        k_Rob_low(i,1) = 1e-5;
        k_Rob_high(i,1) = 1e-3;

    elseif Ic_kDepths(i,1) < 1.31 % (7) Dense sand to gravelly sand
        k_Rob_low(i,1) = 1e-3;
        k_Rob_high(i,1) = 1;
    end
    x = [k_Rob_low(i,1) k_Rob_high(i,1) k_Rob_high(i,1) k_Rob_low(i,1)];
    y = [d_low d_low d_high d_high];
    if i == 1
        fill(x,y,'r',"FaceAlpha",0.1,"LineStyle","none","HandleVisibility","on","DisplayName","Robertson (2010)")
    else
        fill(x,y,'r',"FaceAlpha",0.1,"LineStyle","none","HandleVisibility","off","DisplayName","Robertson (2010)")
    end
    hold on

    % Robertson and Cabal, 2010
    if Ic_kDepths(i,1) > 1.0 && Ic_kDepths(i,1) <= 3.27
        k_RobCab(i,1) = 10 ^ (0.952 - 3.04 * Ic_kDepths(i,1));
    elseif Ic_kDepths(i,1) > 3.27 && Ic_kDepths(i,1) <= 4.0
        k_RobCab(i,1) = 10 ^ (-4.52 - 1.37 * Ic_kDepths(i,1));
    end
end
% k_h after Jamiolkowski et al., 1985:
for i = 1 : length(k_lab)
    dd = depth_lab(i,1);
    int = 0.2;
    d_low = dd - int/2;
    d_high = dd + int/2;
    % As proposed by Jamiolkowski (1985), 3:15 x kv is plotted for kv:
    x = [3*k_lab(i,1) 15*k_lab(i,1) 15*k_lab(i,1) 3*k_lab(i,1)];
    y = [d_low d_low d_high d_high];
    if i == 1
        fill(x,y,'y',"FaceAlpha",0.1,"LineStyle","-","HandleVisibility","on","DisplayName","$k_h$ range from lab $k_v$ after Jamiolkowski et al. (1985)")
    else
        fill(x,y,'y',"FaceAlpha",0.1,"LineStyle","-","HandleVisibility","off","DisplayName","$k_h$ range from lab $k_v$ after Jamiolkowski et al. (1985)")
    end
end
% Robertson and Cabal, 2010
scatter(k_RobCab,depth(1:end-1,1),50,"Marker","square","MarkerEdgeColor","blue","HandleVisibility","on","DisplayName","Robertson and Cabal (2010)");

xLab = "Permeability, $k$ [m/s]";
yLab = "Depth [m]";
xlabel(xLab,'FontSize',10,'Color','k','Interpreter','latex')
ylabel(yLab,'FontSize',10,'Color','k','Interpreter','latex')
xlim([1e-10 1e-2])
ylim([1 7])
set(gca,'XScale','log')
legend('FontSize',8, 'Location','eastoutside', 'Interpreter','latex', "color","white");
grid on
grid minor
ax = gca;
set(ax,'TickLabelInterpreter','latex')
set (gca, 'YDir','reverse')
set(gca,'box','on')


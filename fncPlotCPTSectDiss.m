function LocalParams = fncPlotCPTSectDiss(depth, ddiss, qc, fs, u2, Ic,SectRange, u,t,leg,i)
%it plots the qc,fs,u2,Ic versus depth in local sections where dissipation
%is performed.

% Step 1: Find the index of the row matching 'd' in the first column
rowIndex = find(depth(:, 1) == ddiss);

% Step 2: Define the range for SectRange rows before and after
startIndex = max(1, rowIndex - SectRange); % Ensure we don't go below the first row
endIndex = min(size(depth, 1), rowIndex + SectRange); % Ensure we don't go beyond the last row

% Step 3: Extract the submatrix
depthSect = depth(startIndex:endIndex); % Extract rows in the range
qcSect = qc(startIndex:endIndex); % in MPa
fsSect = 1000*fs(startIndex:endIndex); % in kPa
u2Sect = 1000*u2(startIndex:endIndex); % in kPa
IcSect = Ic(startIndex:endIndex);

LocalParams = [depthSect, qcSect, fsSect, u2Sect, IcSect];

chartlinewidth = 1.5;

ax = nexttile(1+(i-1)*21,[1,1]);
y = depthSect;
x = qcSect;
pd = fitdist(x, 'Normal'); % Fit a normal distribution
mu = pd.mu; % Mean
sigma = pd.sigma; % Standard deviation
x_gauss = linspace(min(qc), max(qc), 50); % Fine grid for Gaussian curve
f_gauss = (1 / (sqrt(2 * pi) * sigma)) * exp(-((x_gauss - mu).^2) / (2 * sigma^2)); % Gaussian PDF
plot(x_gauss, f_gauss, 'r-', 'LineWidth', chartlinewidth);
title(sprintf('$\\mu$ = %.2f', mu), 'FontSize', 8, 'Interpreter', 'latex');
xlim([0 3.5])
set(ax, 'XTickLabel', []); % Remove labels only
set(ax, 'YTickLabel', []); % Remove labels only
grid on, grid minor
box on


ax = nexttile(8+(i-1)*21,[2,1]);
y = depthSect;
x = qcSect;
plot(x, y, 'k-', 'LineWidth', chartlinewidth);
ylabel('Depth (m)','FontSize',10,'Interpreter','latex');
ax.XLim = [0 3.5];
hold on
text(qc((startIndex+endIndex)/2), (ax.YLim(1)+ax.YLim(2))/2, num2str(round(qc((startIndex+endIndex)/2),2)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'latex');
hold off
ylim([depth(startIndex) depth(endIndex)])
set(gca, "YDir","reverse")
set(ax,'TickLabelInterpreter','latex')
grid on, grid minor
box on


nexttile(2+(i-1)*21,[1,1])
y = depthSect;
x = fsSect;
pd = fitdist(x, 'Normal'); % Fit a normal distribution
mu = pd.mu; % Mean
sigma = pd.sigma; % Standard deviation
x_gauss = linspace(1000*min(fs), 1000*max(fs), 50); % Fine grid for Gaussian curve
f_gauss = (1 / (sqrt(2 * pi) * sigma)) * exp(-((x_gauss - mu).^2) / (2 * sigma^2)); % Gaussian PDF
plot(x_gauss, f_gauss, 'r-', 'LineWidth', chartlinewidth);
title(sprintf('$\\mu$ = %.2f', mu), 'FontSize', 8, 'Interpreter', 'latex');
xlim([0 80])
set(gca, 'XTickLabel', []); % Remove labels only
set(gca, 'YTickLabel', []); % Remove labels only
grid on, grid minor
box on

ax = nexttile(9+(i-1)*21,[2,1]);
y = depthSect;
x = fsSect;
plot(x, y, 'k-', 'LineWidth', chartlinewidth);
upperX = 80;
ax.XLim = [0 upperX];
hold on
text(1000*fs((startIndex+endIndex)/2), (ax.YLim(1)+ax.YLim(2))/2, num2str(round(1000*fs((startIndex+endIndex)/2),2)), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'FontWeight', 'bold','Interpreter','latex');
hold off
ylim([depth(startIndex) depth(endIndex)])
set(gca, "YDir","reverse")
set(gca, 'YTickLabel', []); % Remove labels only
set(ax,'TickLabelInterpreter','latex')
grid on, grid minor
box on


nexttile(3+(i-1)*21,[1,1])
y = depthSect;
x = u2Sect;
pd = fitdist(x, 'Normal'); % Fit a normal distribution
mu = pd.mu; % Mean
sigma = pd.sigma; % Standard deviation
x_gauss = linspace(min(x), max(x), 50); % Fine grid for Gaussian curve
f_gauss = (1 / (sqrt(2 * pi) * sigma)) * exp(-((x_gauss - mu).^2) / (2 * sigma^2)); % Gaussian PDF
plot(x_gauss, f_gauss, 'r-', 'LineWidth', chartlinewidth);
title(sprintf('$\\mu$ = %.2f', mu), 'FontSize', 8, 'Interpreter', 'latex');
xlim([-80 200])
set(gca, 'XTickLabel', []); % Remove labels only
set(gca, 'YTickLabel', []); % Remove labels only
grid on, grid minor
box on

ax = nexttile(10+(i-1)*21,[2,1]);
y = depthSect;
x = u2Sect;
plot(x, y, 'k-', 'LineWidth', chartlinewidth);
ax.XLim = [-80 200];
hold on
text(u2((startIndex+endIndex)/2), (ax.YLim(1)+ax.YLim(2))/2, num2str(round(1000*u2((startIndex+endIndex)/2),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8, 'FontWeight', 'bold', 'Interpreter','latex');
hold off
ylim([depth(startIndex) depth(endIndex)])
set(gca, "YDir","reverse")
set(gca, 'YTickLabel', []); % Remove labels only
set(ax,'TickLabelInterpreter','latex')
grid on, grid minor
box on


nexttile(4+(i-1)*21,[1,1])
y = depthSect;
x = IcSect;
pd = fitdist(x, 'Normal'); % Fit a normal distribution
mu = pd.mu; % Mean
sigma = pd.sigma; % Standard deviation
x_gauss = linspace(min(Ic), max(Ic), 50); % Fine grid for Gaussian curve
f_gauss = (1 / (sqrt(2 * pi) * sigma)) * exp(-((x_gauss - mu).^2) / (2 * sigma^2)); % Gaussian PDF
plot(x_gauss, f_gauss, 'r-', 'LineWidth', chartlinewidth);
title(sprintf('$\\mu$ = %.2f', mu), 'FontSize', 8, 'Interpreter', 'latex');
xlim([2.0 3.0])
set(gca, 'XTickLabel', []); % Remove labels only
set(gca, 'YTickLabel', []); % Remove labels only
grid on, grid minor
box on

ax = nexttile(11+(i-1)*21,[2,1]);
y = depthSect;
x = IcSect;
plot(x, y, 'k-', 'LineWidth', chartlinewidth);
ax.XLim = [2.0 3.0];
hold on
text(Ic((startIndex+endIndex)/2), (ax.YLim(1)+ax.YLim(2))/2, num2str(round(Ic((startIndex+endIndex)/2),2)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 8, 'FontWeight', 'bold','Interpreter','latex');
hold off
ylim([depth(startIndex) depth(endIndex)])
set(gca, "YDir","reverse")
set(gca, 'YTickLabel', []); % Remove labels only
set(ax,'TickLabelInterpreter','latex')
grid on, grid minor
box on


nexttile(5+(i-1)*21,[3,3])
y = u;
x = t;
semilogx(x, y, 'k-', 'LineWidth', chartlinewidth);
ylabel('$u_2$ [kPa]','FontSize',10,'Interpreter','latex');
xlabel('Time [sec]','FontSize',10,'Interpreter','latex');
DepthLeg = vpa(strrep(leg,'d',''))/100;
lege = sprintf('d=%.2f m',DepthLeg);
legend(lege, 'FontSize',8, 'Location','best', 'Interpreter','latex', "color","white");
set(gca, "YDir","normal")
xLabelHandle = gca().XLabel;
xLabelHandle.Units = 'normalized';
xLabelHandle.Position = [0.87, 0.17, 0]; % Center horizontally, slightly below the x-axis
yLabelHandle = gca().YLabel;
yLabelHandle.Units = 'normalized';
yLabelHandle.Position = [0.08, 0.76, 0]; % Center horizontally, slightly below the x-axis
set(gca,'TickLabelInterpreter','latex')
grid on, grid minor
box on
end
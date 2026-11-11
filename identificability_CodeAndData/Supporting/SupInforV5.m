% SupInforV5.m
% Input data: ../Data/ (columns Day, NewCase; Day 0 = 21 January 2020)
% Output files: ../Figures/Sup/FigureS1_R_Degeneracy_Multi.png/.pdf
%               ../Figures/Sup/FigureS2_R0_Degeneracy_Multi.png/.pdf
%               ../Figures/Sup/FigureS2_RealData_SideBySide.png/.pdf
%               ../Figures/Sup/FigureS4_Synthetic_TwoParamFit.png/.pdf
%               ../Figures/Sup/FigureS5_RealData_TwoParamFit.png/.pdf
%               ../Results/betaAB_TwoParamResults.csv, betaAB_table.tex
% Local functions (end of file): SEAIR_model, model_incidence, SSE_beta,
%               SSE_beta2, fit_two_beta, aic_gauss, estimate_growth,
%               poisson_growth_nll, get_window

clear; clc; close all;

set(groot,'defaultTextInterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

projectDir = pwd;

while ~exist(fullfile(projectDir,'Data'),'dir')
    parentDir = fileparts(projectDir);

    if strcmp(parentDir,projectDir)
        error('Project root not found. Could not locate the Data folder.');
    end

    projectDir = parentDir;
end

dataDir = fullfile(projectDir,'Data');
figDir  = fullfile(projectDir,'Figures');
supDir  = fullfile(figDir,'Sup');
resDir  = fullfile(projectDir,'Results');

if ~exist(supDir,'dir')
    mkdir(supDir);
end

if ~exist(resDir,'dir')
    mkdir(resDir);
end

fprintf('\nProject root: %s\n',projectDir);
fprintf('Data directory: %s\n',dataDir);
fprintf('Figure directory: %s\n',figDir);
fprintf('Supporting figures: %s\n',supDir);
fprintf('Results directory: %s\n\n',resDir);




dataFile = fullfile(dataDir, 'real_data_window.csv');


T = readtable(dataFile);

y = double(T.new_cases(:));


if isdatetime(T.date)
    dates = T.date;
else
    dates = datetime(string(T.date), ...
        'InputFormat', 'yyyy-MM-dd');
end


nDays = numel(y);

if nDays ~= 70
    warning('Expected 70 observations, but found %d.', nDays);
end


cumCases = cumsum(y);
threshold = 5;

j = find(cumCases >= threshold, 1, 'first');

if isempty(j)
    error('Cumulative cases do not reach the threshold of five.');
end

t = (1:nDays)';

tFit = t(j:end);
yFit = y(j:end);


nll = @(theta) sum( ...
    exp(theta(1) + theta(2).*tFit) ...
    - yFit .* (theta(1) + theta(2).*tFit) );


theta0 = [log(max(mean(yFit), 1e-3)), 0.1];

opts = optimset( ...
    'Display', 'off', ...
    'TolX', 1e-9, ...
    'TolFun', 1e-12, ...
    'MaxIter', 50000, ...
    'MaxFunEvals', 100000);

thetaHat = fminsearch(nll, theta0, opts);

aHat = thetaHat(1);
rHat = thetaHat(2);



fprintf('\n--- Poisson growth-rate estimation ---\n');

fprintf('Analysis window: %s to %s\n', ...
    datestr(dates(1), 'dd mmm yyyy'), ...
    datestr(dates(end), 'dd mmm yyyy'));

fprintf('Total observations: %d\n', nDays);

fprintf('Cumulative case threshold: %d\n', threshold);

fprintf('Threshold reached on day %d (%s)\n', ...
    j, datestr(dates(j), 'dd mmm yyyy'));

fprintf('Observations used in regression: %d\n', numel(yFit));

fprintf('Regression time indices: %d to %d\n', ...
    tFit(1), tFit(end));

fprintf('Estimated intercept (a): %.8f\n', aHat);

fprintf('Estimated growth rate (r): %.8f per day\n', rHat);

fprintf('Reported growth rate: %.4f per day\n', rHat);






mu    = 1/(80.3*365);
muS   = 0.00159;
p     = 0.956;
eta   = 5.5;
epsi  = 1/2.381;
omega = 1/(eta - 1/epsi);
P = struct('mu',mu,'muS',muS,'p',p,'epsi',epsi,'omega',omega);

S0 = 331002647;
IS0 = 10;  IA0 = IS0;  E0 = 2.5*(IA0+IS0);
Y0 = [S0 E0 IA0 IS0 0 0];

colA = [230 159 0]/255;
colB = [0 114 178]/255;
colA4 = [  0 158 115]/255;
colB4 = [204 121 167]/255;
nu_B = 0.05;
nu_A = 0.30;

nu_grid = linspace(0.03,0.40,500);

starts = [0.2 0.001; 0.2 0.5; 0.2 2; 0.2 4; 0.5 0.001; 0.5 0.5; 0.5 2; 0.5 4; 0.8 0.001; 0.8 0.5; 0.8 2; 0.8 4];

R0_of = @(bA,nu,bS) epsi*(bA*(muS+nu)+bS*(1-p)*omega) / ((mu+epsi)*(mu+omega)*(muS+nu));


nDays = 70;

TblReal = readtable(fullfile(dataDir,'real_data_window.csv'), ...
    'TextType','string');

if isdatetime(TblReal.date)
    realDates = TblReal.date;
else
    realDates = datetime(TblReal.date,'InputFormat','yyyy-MM-dd');
end

[realDates,order] = sort(realDates);
TblReal = TblReal(order,:);
y_real = double(TblReal.new_cases(:));
t_real = (1:nDays)';
r_hat = estimate_growth(y_real,5);



r_values = [0.05, 0.10, 0.15, r_hat, 0.25, 0.30];

r_labels = {
    '$r = 0.05$'
    'r = 0.10'
    'r = 0.15'
    sprintf('$\\hat{r}$ = %.4f, estimated from real data',r_hat)
    'r = 0.25'
    'r = 0.30'
};

viridis_anchors = [
    0.267004 0.004874 0.329415
    0.278826 0.175490 0.483397
    0.229739 0.322361 0.545706
    0.172719 0.448791 0.557885
    0.127568 0.566949 0.550556
    0.157851 0.683765 0.501686
    0.369214 0.788888 0.382914
    0.678489 0.863742 0.189503
    0.993248 0.906157 0.143936
];

viridis_positions = linspace(0,1,size(viridis_anchors,1));
color_positions = linspace(0,1,numel(r_values));
colors = interp1(viridis_positions,viridis_anchors,color_positions);


fig1 = figure('Color','w','Units','inches','Position',[1 1 6.5 4.8]);
hold on

hLines = gobjects(1,numel(r_values));
all_bc = [];

for i = 1:numel(r_values)

    r = r_values(i);

    
    bc = ((mu+epsi+r).*(mu+omega+r).*(muS+nu_grid+r)) ./ ...
         (epsi.*((muS+nu_grid+r)+(1-p).*omega));

    hLines(i) = plot(nu_grid,bc,'-', ...
        'LineWidth',2,'Color',colors(i,:));

    all_bc = [all_bc, bc(:)'];

    
    bA_r = ((mu+epsi+r)*(mu+omega+r)*(muS+nu_A+r)) / ...
           (epsi*((muS+nu_A+r)+(1-p)*omega));

    bB_r = ((mu+epsi+r)*(mu+omega+r)*(muS+nu_B+r)) / ...
           (epsi*((muS+nu_B+r)+(1-p)*omega));

    plot(nu_A,bA_r,'s','MarkerSize',10, ...
        'MarkerFaceColor','k','MarkerEdgeColor','k');

    plot(nu_B,bB_r,'o','MarkerSize',10, ...
        'MarkerFaceColor','k','MarkerEdgeColor','k');

    all_bc = [all_bc, bA_r, bB_r];

end


ymin = min(all_bc);
ymax = max(all_bc);
ylim([ymin, ymax + 0.25*(ymax-ymin)]);


xlabel('\nu (recovery rate)','Interpreter','tex');
ylabel('\beta (transmission rate)','Interpreter','tex');


lgd = legend(hLines,r_labels, ...
    'Location','northwest', ...
    'NumColumns',2, ...
    'Interpreter','latex', ...
    'FontSize',10);

lgd.Box = 'on';
lgd.Color = 'white';


box on
grid on
set(gca,'FontSize',13)

exportgraphics(fig1, ...
    fullfile(figDir,'FigureS1_R_Degeneracy_Multi.png'), ...
    'Resolution',600);


%exportgraphics(fig1,fullfile(figDir,'FigureS1_R_Degeneracy_Multi.pdf'),'ContentType','vector');

R0_values = [1.5, 2.0, 2.5, 3.0, 3.5];

R0_labels = {
    '$\mathcal{R}_0 = 1.5$'
    '$\mathcal{R}_0 = 2.0$'
    '$\mathcal{R}_0 = 2.5$'
    '$\mathcal{R}_0 = 3.0$'
    '$\mathcal{R}_0 = 3.5$'
};


color_positions = linspace(0,1,numel(R0_values));
colors_R0 = interp1(viridis_positions,viridis_anchors,color_positions);


fig2 = figure('Color','w','Units','inches','Position',[1 1 6.5 4.8]);
hold on

hLines2 = gobjects(1,numel(R0_values));
all_bc_R0 = [];

for i = 1:numel(R0_values)

    R0t = R0_values(i);

   
    bc = (R0t.*(mu+epsi).*(mu+omega).*(muS+nu_grid)) ./ ...
         (epsi.*((muS+nu_grid)+(1-p).*omega));

    hLines2(i) = plot(nu_grid,bc,'-', ...
        'LineWidth',2,'Color',colors_R0(i,:));

    all_bc_R0 = [all_bc_R0, bc(:)'];

    
    bA_R = (R0t*(mu+epsi)*(mu+omega)*(muS+nu_A)) / ...
           (epsi*((muS+nu_A)+(1-p)*omega));

    bB_R = (R0t*(mu+epsi)*(mu+omega)*(muS+nu_B)) / ...
           (epsi*((muS+nu_B)+(1-p)*omega));

    plot(nu_A,bA_R,'s','MarkerSize',10, ...
        'MarkerFaceColor','k','MarkerEdgeColor','k');

    plot(nu_B,bB_R,'o','MarkerSize',10, ...
        'MarkerFaceColor','k','MarkerEdgeColor','k');

    all_bc_R0 = [all_bc_R0, bA_R, bB_R];

end


ymin = min(all_bc_R0);
ymax = max(all_bc_R0);
ylim([ymin, ymax + 0.25*(ymax-ymin)]);


xlabel('\nu (recovery rate)','Interpreter','tex');
ylabel('\beta (transmission rate)','Interpreter','tex');


lgd2 = legend(hLines2,R0_labels, ...
    'Location','northwest', ...
    'NumColumns',2, ...
    'Interpreter','latex', ...
    'FontSize',10);

lgd2.Box = 'on';
lgd2.Color = 'white';


box on
grid on
set(gca,'FontSize',13)


exportgraphics(fig2, ...
    fullfile(figDir,'FigureS2_R0_Degeneracy_Multi.png'), ...
    'Resolution',600);

%exportgraphics(fig2,fullfile(figDir,'FigureS2_R0_Degeneracy_Multi.pdf'),'ContentType','vector');

tgrid_r = (0:nDays)';
optsFB = optimset('TolX',1e-6);
[bhA,sseA] = fminbnd(@(b) SSE_beta(b,nu_A,tgrid_r,y_real,Y0,P),0.01,3,optsFB);
[bhB,sseB] = fminbnd(@(b) SSE_beta(b,nu_B,tgrid_r,y_real,Y0,P),0.01,3,optsFB);
AICa = aic_gauss(sseA,nDays,1);
AICb = aic_gauss(sseB,nDays,1);
incA = model_incidence(bhA,bhA,nu_A,tgrid_r,Y0,P);
incB = model_incidence(bhB,bhB,nu_B,tgrid_r,Y0,P);

fig3 = figure('Color','w','Units','inches','Position',[1 1 10 4.2]);
subplot(1,2,1);
plot(t_real,y_real,'x','MarkerSize',8,'LineWidth',1.5,'Color','k');
hold on
plot(t_real,incA,'LineWidth',2,'Color',colA);
xlabel('Days'); ylabel('New reported cases per day');
title('High Transmission (A)','FontWeight','bold');
text(0.05,0.9,sprintf('AIC = %.1f',AICa),'Units','normalized','FontWeight','bold', ...
     'FontSize',12,'BackgroundColor','w','EdgeColor','k');
box on; grid on; set(gca,'FontSize',13)
subplot(1,2,2);
plot(t_real,y_real,'x','MarkerSize',8,'LineWidth',1.5,'Color','k');
hold on
plot(t_real,incB,'LineWidth',2,'Color',colB);
xlabel('Days'); ylabel('New reported cases per day');
title('Low Recovery (B)','FontWeight','bold');
text(0.05,0.9,sprintf('AIC = %.1f',AICb),'Units','normalized','FontWeight','bold', ...
     'FontSize',12,'BackgroundColor','w','EdgeColor','k');
box on; grid on; set(gca,'FontSize',13)

exportgraphics(fig3,fullfile(figDir,'FigureS2_RealData_SideBySide.png'),'Resolution',600);
%exportgraphics(fig3,fullfile(figDir,'FigureS2_RealData_SideBySide.pdf'),'ContentType','vector');


beta_T = 0.4897;
nu_T = 0.12;

syntheticFile = fullfile(dataDir,'Synthetic_Observations.csv');

if ~isfile(syntheticFile)
    error(['Synthetic_Observations.csv not found. ' ...
           'Run RunnerSyntheticV5.m first.']);
end

SynTbl = readtable(syntheticFile);

day_syn = double(SynTbl.day(:));
true_incidence_file = double(SynTbl.true_incidence(:));
y_syn = double(SynTbl.observed(:));

nSyn = numel(y_syn);



tspan = (0:nSyn)';
t_obs_s = day_syn;

inc_true_model = model_incidence(beta_T,beta_T,nu_T,tspan,Y0,P);

inc_true = true_incidence_file;




[thA_s,sseA_s] = fit_two_beta(nu_A,tspan,y_syn,Y0,P,starts);
[thB_s,sseB_s] = fit_two_beta(nu_B,tspan,y_syn,Y0,P,starts);
AICA_s = aic_gauss(sseA_s,70,2);
AICB_s = aic_gauss(sseB_s,70,2);

fprintf('Synthetic, A: beta_A=%.3f, beta_S=%.3f, AIC=%.1f\n',thA_s(1),thA_s(2),AICA_s);
fprintf('Synthetic, B: beta_A=%.3f, beta_S=%.3f, AIC=%.1f\n',thB_s(1),thB_s(2),AICB_s);
fprintf('AIC_A - AIC_B = %.3f\n',AICA_s-AICB_s);

incA_s = model_incidence(thA_s(1),thA_s(2),nu_A,tspan,Y0,P);
incB_s = model_incidence(thB_s(1),thB_s(2),nu_B,tspan,Y0,P);


fig4 = figure('Color','w','Units','inches','Position',[1 1 6.5 4.8]);

scatter(t_obs_s,y_syn,35,'k','filled');
hold on

plot(t_obs_s,inc_true,'k--','LineWidth',2);
plot(t_obs_s,incA_s,'LineWidth',2,'Color',colA4);
plot(t_obs_s,incB_s,'LineWidth',2,'Color',colB4);

xlabel('Days');
ylabel('New symptomatic cases per day');

legend('Synthetic observations','True trajectory', ...
    sprintf('A: \\beta_A=%.3f, \\beta_S=%.3f, AIC=%.1f', ...
        thA_s(1),thA_s(2),AICA_s), ...
    sprintf('B: \\beta_A=%.3f, \\beta_S=%.3f, AIC=%.1f', ...
        thB_s(1),thB_s(2),AICB_s), ...
    'Location','northwest','FontSize',10);

box on
grid on
set(gca,'FontSize',13)

exportgraphics(fig4, ...
    fullfile(supDir,'FigureS4_Synthetic_TwoParamFit.png'), ...
    'Resolution',600);

%exportgraphics(fig4,fullfile(figDir,'FigureS4_Synthetic_TwoParamFit.pdf'),'ContentType','vector');

[thA_r,sseA_r] = fit_two_beta(nu_A,tgrid_r,y_real,Y0,P,starts);
[thB_r,sseB_r] = fit_two_beta(nu_B,tgrid_r,y_real,Y0,P,starts);
AICA_r = aic_gauss(sseA_r,nDays,2);
AICB_r = aic_gauss(sseB_r,nDays,2);

fprintf('Real, A: beta_A=%.3f, beta_S=%.3f, AIC=%.1f\n',thA_r(1),thA_r(2),AICA_r);
fprintf('Real, B: beta_A=%.3f, beta_S=%.3f, AIC=%.1f\n',thB_r(1),thB_r(2),AICB_r);
fprintf('AIC_A - AIC_B = %.3f\n',AICA_r-AICB_r);

incA_r = model_incidence(thA_r(1),thA_r(2),nu_A,tgrid_r,Y0,P);
incB_r = model_incidence(thB_r(1),thB_r(2),nu_B,tgrid_r,Y0,P);

fig5 = figure('Color','w','Units','inches','Position',[1 1 6.5 4.8]);
plot(t_real,y_real,'x','MarkerSize',8,'LineWidth',1.5,'Color','k');
hold on
plot(t_real,incA_r,'LineWidth',2,'Color',colA);
plot(t_real,incB_r,'LineWidth',2,'Color',colB);
xlabel('Days'); ylabel('New reported cases per day');
legend('Observations', ...
    sprintf('A: \\beta_A=%.3f, \\beta_S=%.3f, AIC=%.1f',thA_r(1),thA_r(2),AICA_r), ...
    sprintf('B: \\beta_A=%.3f, \\beta_S=%.3f, AIC=%.1f',thB_r(1),thB_r(2),AICB_r), ...
    'Location','northwest','FontSize',10);
box on; grid on
set(gca,'FontSize',13)
exportgraphics(fig5,fullfile(figDir,'FigureS5_RealData_TwoParamFit.png'),'Resolution',600);
%exportgraphics(fig5,fullfile(figDir,'FigureS5_RealData_TwoParamFit.pdf'),'ContentType','vector');

Dataset    = {'Synthetic';'Synthetic';'COVID-19';'COVID-19'};
Hypothesis = {'High Transmission';'Low Recovery';'High Transmission';'Low Recovery'};
betaA_hat  = [thA_s(1); thB_s(1); thA_r(1); thB_r(1)];
betaS_hat  = [thA_s(2); thB_s(2); thA_r(2); thB_r(2)];
AIC_val    = [AICA_s; AICB_s; AICA_r; AICB_r];
nus        = [nu_A; nu_B; nu_A; nu_B];

fid = fopen(fullfile(resDir,'betaAB_TwoParamResults.csv'),'w');
fprintf(fid,'dataset,hypothesis,beta_A_hat,beta_S_hat,R0,AIC\n');
for i = 1:4
    fprintf(fid,'%s,%s,%.4f,%.4f,%.4f,%.2f\n',Dataset{i},Hypothesis{i}, ...
            betaA_hat(i),betaS_hat(i),R0_of(betaA_hat(i),nus(i),betaS_hat(i)),AIC_val(i));
end
fclose(fid);

fid = fopen(fullfile(resDir,'betaAB_table.tex'),'w');
fprintf(fid,'\\begin{table}[!ht]\n\\centering\n');
fprintf(fid,['\\caption{Estimated transmission parameters when $\\beta_A$ and $\\beta_S$ are fit independently ' ...
             '($\\beta_A\\ne\\beta_S$) to daily symptomatic incidence, with $\\nu$ fixed at its characteristic value ' ...
             'for each hypothesis. Estimates are constrained to $0 \\le \\beta \\le 5$.}\n']);
fprintf(fid,'\\label{tab:betaAB}\n');
fprintf(fid,'\\begin{tabular}{lcccc}\n\\toprule\n');
fprintf(fid,'Dataset & Hypothesis & $\\hat{\\beta}_A$ & $\\hat{\\beta}_S$ & AIC \\\\\n\\midrule\n');
for i = 1:4
    fprintf(fid,'%s & %s & %.3f & %.3f & %.1f \\\\\n',Dataset{i},Hypothesis{i},betaA_hat(i),betaS_hat(i),AIC_val(i));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
fclose(fid);
%%
function dydt = SEAIR_model(~,y,q)
S  = y(1); E = y(2); IA = y(3); IS = y(4); R = y(5);
N = S+E+IA+IS+R;
lambda = (q.betaA*IA + q.betaS*IS)/N;
onset  = (1-q.p)*q.omega*IA;
dS  = q.mu*(N-S) - lambda*S;
dE  = lambda*S - (q.epsi+q.mu)*E;
dIA = q.epsi*E - (q.omega+q.mu)*IA;
dIS = onset - (q.nu+q.muS)*IS;
dR  = q.p*q.omega*IA + q.nu*IS - q.mu*R;
dC  = onset;
dydt = [dS;dE;dIA;dIS;dR;dC];
end

function inc = model_incidence(betaA,betaS,nu,tgrid,Y0,P)
q = P; q.betaA = betaA; q.betaS = betaS; q.nu = nu;
opts = odeset('RelTol',1e-8,'AbsTol',1e-6);
[~,Y] = ode45(@(t,y) SEAIR_model(t,y,q), tgrid, Y0, opts);
inc = diff(Y(:,6));
end

function SSE = SSE_beta(beta,nu,tgrid,yobs,Y0,P)
try
    inc = model_incidence(beta,beta,nu,tgrid,Y0,P);
    SSE = sum((yobs(:)-inc(:)).^2);
    if ~isfinite(SSE), SSE = 1e20; end
catch
    SSE = 1e20;
end
end

function SSE = SSE_beta2(theta,nu,tgrid,yobs,Y0,P)
try
    inc = model_incidence(theta(1),theta(2),nu,tgrid,Y0,P);
    SSE = sum((yobs(:)-inc(:)).^2);
    if ~isfinite(SSE), SSE = 1e20; end
catch
    SSE = 1e20;
end
end

function [theta,SSEmin] = fit_two_beta(nu,tgrid,yobs,Y0,P,starts)
lo = 0;  hi = 5;
toBeta = @(u) lo + (hi-lo)*(sin(u)+1)/2;
toU    = @(b) asin(max(min(2*(b-lo)/(hi-lo)-1,1),-1));
opts = optimset('Display','off','TolX',1e-6,'TolFun',1e-9,'MaxFunEvals',600,'MaxIter',600);
SSEmin = Inf;  theta = [NaN NaN];
for s = 1:size(starts,1)
    obj = @(u) SSE_beta2(toBeta(u),nu,tgrid,yobs,Y0,P);
    [u,fval] = fminsearch(obj,toU(starts(s,:)),opts);
    if fval < SSEmin
        SSEmin = fval;
        theta  = toBeta(u);
    end
end
end

function AIC = aic_gauss(SSE,n,k)
logLik = -0.5*n*(log(2*pi*SSE/n)+1);
AIC = 2*k - 2*logLik;
end

function [r,j] = estimate_growth(y,thr)
y = y(:);  n = numel(y);
j = find(cumsum(y) >= thr,1,'first');
if isempty(j), j = 1; end
t = (j:n)';  yy = y(j:n);
opts = optimset('Display','off','TolX',1e-9,'TolFun',1e-10,'MaxFunEvals',8000,'MaxIter',8000);
th = fminsearch(@(q) poisson_growth_nll(q,t,yy),[0 0.05],opts);
r = th(2);
end

function nll = poisson_growth_nll(theta,t,y)
a = theta(1); r = theta(2);
lam = max(exp(a + r*t),1e-10);
nll = -sum(y.*log(lam) - lam - gammaln(y+1));
end

function y = get_window(counts,startDay,nDays)
days = (startDay:startDay+nDays-1)';
if days(end) > numel(counts)-1
    error('Window ends after the last day in Data120.csv.');
end
y = zeros(nDays,1);
ok = days >= 0;
y(ok) = counts(days(ok)+1);
end

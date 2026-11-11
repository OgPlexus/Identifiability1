% RunnerRealDataV5.m
% Input data: ../Data/Data120.csv (columns Day, NewCase; Day 0 = 21 January 2020)
% Output files: ../Figures/Figure3A_RealData_Degeneracy.png/.pdf
%               ../Figures/Figure3B_RealData_Fits.png/.pdf
%               ../Figures/Sup/RealDataOnly.png
%               ../Results/RealDataSummaryTable.csv/.tex
%               ../Results/RealDataWindowSensitivity.csv
% Local functions (end of file): SEAIR_model, model_incidence, SSE_beta,
%               fit_real, aic_gauss, estimate_growth, poisson_growth_nll,
%               get_window, write_table
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

%

nDays = 70;
IS0   = 10;
burnInThreshold = 5;

mu    = 1/(80.3*365);
muS   = 0.00159;
p     = 0.956;
eta   = 5.5;
epsi  = 1/2.381;
omega = 1/(eta - 1/epsi);
P = struct('mu',mu,'muS',muS,'p',p,'epsi',epsi,'omega',omega);

S0 = 331002647;
colA = [230 159 0]/255;
colB = [0 114 178]/255;
nu_B = 0.05;
nu_A = 0.30;

beta_closed = @(nu,r) ((mu+epsi+r)*(mu+omega+r)*(muS+nu+r)) / (epsi*((muS+nu+r)+(1-p)*omega));
R0_of = @(b,nu) epsi*(b*(muS+nu)+b*(1-p)*omega) / ((mu+epsi)*(mu+omega)*(muS+nu));


TblReal = readtable(fullfile(dataDir,'real_data_window.csv'), ...
    'TextType','string');


if isdatetime(TblReal.date)
    realDates = TblReal.date;
else
    realDates = datetime(TblReal.date,'InputFormat','yyyy-MM-dd');
end

[realDates,order] = sort(realDates);
TblReal = TblReal(order,:);
y_obs = double(TblReal.new_cases(:));
t_obs = (1:nDays)';
[r_hat,j5] = estimate_growth(y_obs,burnInThreshold);
fprintf('r_hat = %.4f per day (70-day real-data window)\n',r_hat);


rstar = r_hat;
beta_A_closed = beta_closed(nu_A,rstar);
beta_B_closed = beta_closed(nu_B,rstar);

nu_grid = linspace(0.03,0.40,500);
beta_curve = ((mu+epsi+rstar).*(mu+omega+rstar).*(muS+nu_grid+rstar)) ./ ...
             (epsi.*((muS+nu_grid+rstar)+(1-p).*omega));

figure('Color','w','Units','inches','Position',[1 1 6 4.5]);
plot(nu_grid,beta_curve,'k','LineWidth',2);
hold on
plot(nu_A,beta_A_closed,'s','MarkerSize',12,'MarkerFaceColor',colA,'MarkerEdgeColor','k');
plot(nu_B,beta_B_closed,'o','MarkerSize',12,'MarkerFaceColor',colB,'MarkerEdgeColor','k');
text(nu_A+0.01,beta_A_closed,'(A)','FontSize',12,'FontWeight','bold');
text(nu_B+0.01,beta_B_closed,'(B)','FontSize',12,'FontWeight','bold');
xlabel('\nu (recovery rate)');
ylabel('\beta (transmission rate)');
text(0.215,0.763,sprintf('$\\mathbf{\\hat{r} = %.4f }$ \\textbf{per day,}\n \\textbf{estimated from data}',rstar), ...
     'FontSize',12,'Interpreter','latex');
legend({'Degeneracy curve','High transmission (A)','Low recovery (B)'},'Location','southeast');
box on; grid on
set(gca,'FontSize',14)

exportgraphics(gcf,fullfile(figDir,'Figure3A_RealData_Degeneracy.png'),'Resolution',600);
%exportgraphics(gcf,fullfile(figDir,'Figure3A_RealData_Degeneracy.pdf'),'ContentType','vector');

f = fit_real(y_obs,nu_A,nu_B,S0,IS0,P);
R0_A = R0_of(f.betaA,nu_A);
R0_B = R0_of(f.betaB,nu_B);

fprintf('A: beta_hat = %.4f, R0 = %.4f, AIC = %.1f\n',f.betaA,R0_A,f.AICA);
fprintf('B: beta_hat = %.4f, R0 = %.4f, AIC = %.1f\n',f.betaB,R0_B,f.AICB);
fprintf('AIC_A - AIC_B = %.2f (n = %d)\n',f.AICA-f.AICB,f.n);

fig = figure('Color','w','Units','inches','Position',[1 1 6 4.5]);
hold on
scatter(t_obs,y_obs,32,'o','MarkerFaceColor','k','MarkerEdgeColor','w','LineWidth',0.5);
xlabel('Days','FontSize',13)
ylabel('New reported cases per day','FontSize',13)
xlim([0 nDays+1])
ylim([0 1.05*max(y_obs)])
box on; grid on
set(gca,'FontSize',14,'LineWidth',1.2,'TickDir','out','Layer','top')
legend({'Observed data'},'Location','northwest','Box','on','FontSize',11)

exportgraphics(fig,fullfile(figDir,'Sup','RealDataOnly.png'),'Resolution',600);

figure('Color','w','Units','inches','Position',[1 1 6 4.5]);
plot(t_obs,y_obs,'x','MarkerSize',8,'LineWidth',1.5,'Color','k');
hold on
plot(t_obs,f.incA,'LineWidth',2,'Color',colA);
plot(t_obs,f.incB,'LineWidth',2,'Color',colB);
xlabel('Days');
ylabel('New reported cases per day');
legend('Observed data', ...
       sprintf('High transmission (A), AIC=%.1f',f.AICA), ...
       sprintf('Low recovery (B), AIC=%.1f',f.AICB),'Location','northwest');
box on; grid on
set(gca,'FontSize',14)

exportgraphics(gcf,fullfile(figDir,'Figure3B_RealData_Fits.png'),'Resolution',600);
%exportgraphics(gcf,fullfile(figDir,'Figure3B_RealData_Fits.pdf'),'ContentType','vector');

header = {'Case','r_hat','nu','beta_closed_form','beta_fit','R0','AIC'};
rows = {'High Transmission (A)',r_hat,nu_A,beta_A_closed,f.betaA,R0_A,f.AICA;
        'Low Recovery (B)',     r_hat,nu_B,beta_B_closed,f.betaB,R0_B,f.AICB};
write_table(fullfile(resDir,'RealDataSummaryTable'),header,rows, ...
    ['Real data summary. Closed form beta is the value predicted by the degeneracy curve at the stated nu and the estimated growth rate. ' ...
     'Fitted beta is recovered by fitting daily symptomatic incidence to the real data (n = 70). ' ...
     'AIC values are comparable only within this dataset.'],'tab:realdata_summary');
%%
TblOriginal = readtable( ...
    fullfile(dataDir,'owid_archived_ecdc_us_daily_cases.csv'), ...
    'TextType','string');

originalDates = datetime(TblOriginal.date,'InputFormat','yyyy-MM-dd');
originalCases = double(TblOriginal.new_cases(:));


[originalDates,order] = sort(originalDates);
originalCases = originalCases(order);


referenceDate = datetime(2020,1,21);
counts = originalCases(originalDates >= referenceDate);


startDays = [-2 -1 0 1 2 4];
sweep = zeros(numel(startDays),5);

for i = 1:numel(startDays)

    yi = get_window(counts,startDays(i),nDays);

    ri = estimate_growth(yi,burnInThreshold);
    fi = fit_real(yi,nu_A,nu_B,S0,IS0,P);

    sweep(i,:) = [startDays(i),ri,fi.betaA,fi.betaB, ...
                  fi.AICA-fi.AICB];
end


fprintf('\nReal-data window sensitivity analysis\n');
fprintf('Start day, r_hat, beta_A, beta_B, AIC_A - AIC_B\n');
fprintf('%4d  %.4f  %.4f  %.4f  %.2f\n',sweep');


fid = fopen(fullfile(resDir,'RealDataWindowSensitivity.csv'),'w');

fprintf(fid,'start_day,r_hat,beta_A,beta_B,AIC_A_minus_AIC_B\n');
fprintf(fid,'%d,%.5f,%.5f,%.5f,%.4f\n',sweep');

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

function f = fit_real(y,nu_A,nu_B,S0,IS0,P)
y = y(:);  n = numel(y);
tgrid = (0:n)';
IA0 = IS0;  E0 = 2.5*(IA0+IS0);
Y0 = [S0 E0 IA0 IS0 0 0];
optsFB = optimset('TolX',1e-6);
[bA,sA] = fminbnd(@(b) SSE_beta(b,nu_A,tgrid,y,Y0,P),0.01,3,optsFB);
[bB,sB] = fminbnd(@(b) SSE_beta(b,nu_B,tgrid,y,Y0,P),0.01,3,optsFB);
f.betaA = bA;  f.betaB = bB;  f.n = n;
f.AICA = aic_gauss(sA,n,1);
f.AICB = aic_gauss(sB,n,1);
f.incA = model_incidence(bA,bA,nu_A,tgrid,Y0,P);
f.incB = model_incidence(bB,bB,nu_B,tgrid,Y0,P);
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

function write_table(base,header,rows,captionStr,labelStr)
fid = fopen([base '.csv'],'w');
fprintf(fid,'%s\n',strjoin(header,','));
for i = 1:size(rows,1)
    fprintf(fid,'%s\n',strjoin(cellfun(@fmt_cell,rows(i,:),'UniformOutput',false),','));
end
fclose(fid);
fid = fopen([base '.tex'],'w');
fprintf(fid,'\\begin{table}[!ht]\n\\centering\n\\begin{tabular}{%s}\n\\toprule\n',repmat('c',1,numel(header)));
fprintf(fid,'%s \\\\\n\\midrule\n',strjoin(strrep(header,'_',' '),' & '));
for i = 1:size(rows,1)
    fprintf(fid,'%s \\\\\n',strjoin(cellfun(@fmt_cell,rows(i,:),'UniformOutput',false),' & '));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\caption{%s}\n\\label{%s}\n\\end{table}\n',captionStr,labelStr);
fclose(fid);
end

function s = fmt_cell(v)
if ischar(v)
    s = v;
elseif isnan(v)
    s = '--';
else
    s = sprintf('%.4f',v);
end
end

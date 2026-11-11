% RunnerSyntheticV5.m
% Input data: none (synthetic data are generated in this script)
% Output files: ../Figures/Figure2A_Synthetic_Degeneracy.png/.pdf
%               ../Figures/Figure2B_Synthetic_Fits.png/.pdf
%               ../Figures/Sup/SyntheticData.png
%               ../Results/SyntheticSummaryTable.csv/.tex
%               ../Results/SyntheticReplicates.csv
% Local functions (end of file): SEAIR_model, model_incidence, SSE_beta,
%               aic_gauss, write_table

clear; clc; close all;

set(groot,'defaultTextInterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

% Find project root
projectDir = pwd;

while ~exist(fullfile(projectDir,'Data'),'dir')
    parentDir = fileparts(projectDir);

    if strcmp(parentDir,projectDir)
        error('Could not locate project root.');
    end

    projectDir = parentDir;
end

% Define output directories
dataDir = fullfile(projectDir,'Data');
figDir = fullfile(projectDir,'Figures');
supDir = fullfile(figDir,'Sup');
resDir = fullfile(projectDir,'Results');

% Create directories
if ~exist(supDir,'dir'), mkdir(supDir); end
if ~exist(resDir,'dir'), mkdir(resDir); end

fprintf('Project root: %s\n',projectDir);
fprintf('Figures: %s\n',figDir);
fprintf('Results: %s\n',resDir);
%

nRep = 200;

mu    = 1/(80.3*365);
muS   = 0.00159;
p     = 0.956;
eta   = 5.5;
epsi  = 1/2.381;
omega = 1/(eta - 1/epsi);
P = struct('mu',mu,'muS',muS,'p',p,'epsi',epsi,'omega',omega);

colA = [  0 158 115]/255;
colB = [204 121 167]/255;

rstar = 0.10;
nu_B = 0.05;   nu_T = 0.12;   nu_A = 0.30;
beta_T = 0.4897;

beta_closed = @(nu) ((mu+epsi+rstar)*(mu+omega+rstar)*(muS+nu+rstar)) / ...
                    (epsi*((muS+nu+rstar)+(1-p)*omega));
R0_of = @(beta,nu) epsi*(beta*(muS+nu)+beta*(1-p)*omega) / ((mu+epsi)*(mu+omega)*(muS+nu));

beta_B = beta_closed(nu_B);
beta_A = beta_closed(nu_A);

nu_grid    = linspace(0.03,0.40,500);
beta_curve = ((mu+epsi+rstar).*(mu+omega+rstar).*(muS+nu_grid+rstar)) ./ ...
             (epsi.*((muS+nu_grid+rstar)+(1-p).*omega));
%%
figure('Color','w','Units','inches','Position',[1 1 6 4.5]);
plot(nu_grid,beta_curve,'k','LineWidth',2);
hold on
plot(nu_A,beta_A,'s','MarkerSize',12,'MarkerFaceColor',colA,'MarkerEdgeColor','k');
plot(nu_B,beta_B,'o','MarkerSize',12,'MarkerFaceColor',colB,'MarkerEdgeColor','k');
text(nu_A+0.01,beta_A-0.001,'(A)','FontSize',12,'FontWeight','bold');
text(nu_B+0.01,beta_B,'(B)','FontSize',12,'FontWeight','bold');
text(0.2,0.495,sprintf('Every point on this curve gives \n r = %.2f per day',rstar), ...
     'FontSize',11,'FontWeight','bold');
xlabel('\nu (recovery rate)','Interpreter','tex');
ylabel('\beta (transmission rate)','Interpreter','tex');
legend({'Degeneracy curve','High transmission (A)','Low recovery (B)'},'Location','southeast');
box on; grid on
set(gca,'FontSize',14)
%
exportgraphics(gcf,fullfile(figDir,'Figure2A_Synthetic_Degeneracy.png'),'Resolution',600);
%exportgraphics(gcf,fullfile(figDir,'Figure2A_Synthetic_Degeneracy.pdf'),'ContentType','vector');
%%
IS0 = 10;  IA0 = IS0;  E0 = 2.5*(IA0+IS0);  S0 = 331002647;
Y0 = [S0 E0 IA0 IS0 0 0];

Tmax  = 70;


syntheticFile = fullfile(dataDir,'Synthetic_Observations.csv');

tspan = (0:Tmax)';
t_obs = (1:Tmax)';

if isfile(syntheticFile)

    SynTbl = readtable(syntheticFile);
    inc_true = SynTbl.true_incidence(:);
    y_obs = SynTbl.observed(:);

    fprintf('Existing synthetic dataset loaded: %s\n',syntheticFile);

else

    inc_true = model_incidence(beta_T,beta_T,nu_T,tspan,Y0,P);

    rng(2026,'twister');
    y_obs = poissrnd(max(inc_true,1e-6));

    SyntheticData = table( ...
        (1:Tmax)', inc_true(:), double(y_obs(:)), ...
        'VariableNames',{'day','true_incidence','observed'});

    writetable(SyntheticData,syntheticFile);

    fprintf('New synthetic dataset generated and saved: %s\n',syntheticFile);

end

bmin = 0.01;  bmax = 3;
optsFB = optimset('TolX',1e-6);
[beta_hat_A,SSE_A] = fminbnd(@(b) SSE_beta(b,nu_A,tspan,y_obs,Y0,P),bmin,bmax,optsFB);
[beta_hat_B,SSE_B] = fminbnd(@(b) SSE_beta(b,nu_B,tspan,y_obs,Y0,P),bmin,bmax,optsFB);

inc_fit_A = model_incidence(beta_hat_A,beta_hat_A,nu_A,tspan,Y0,P);
inc_fit_B = model_incidence(beta_hat_B,beta_hat_B,nu_B,tspan,Y0,P);

n = numel(y_obs);
AIC_A = aic_gauss(SSE_A,n,1);
AIC_B = aic_gauss(SSE_B,n,1);

R0_T = R0_of(beta_T,nu_T);
R0_A = R0_of(beta_hat_A,nu_A);
R0_B = R0_of(beta_hat_B,nu_B);

fprintf('A: beta_hat = %.4f, R0 = %.4f, AIC = %.1f\n',beta_hat_A,R0_A,AIC_A);
fprintf('B: beta_hat = %.4f, R0 = %.4f, AIC = %.1f\n',beta_hat_B,R0_B,AIC_B);
fprintf('AIC_A - AIC_B = %.2f (n = %d)\n',AIC_A-AIC_B,n);
%%


rng(12345,'twister');

optsRep = optimset('Display','off','TolX',1e-7);

dAIC_rep = zeros(nRep,1);
bA_rep = zeros(nRep,1);
bB_rep = zeros(nRep,1);

for i = 1:nRep
    yr = poissrnd(max(inc_true(:),1e-6));
    [bA,sA] = fminbnd( ...
        @(b) SSE_beta(b,nu_A,tspan,yr,Y0,P), ...
        0.01,3,optsRep);
    [bB,sB] = fminbnd( ...
        @(b) SSE_beta(b,nu_B,tspan,yr,Y0,P), ...
        0.01,3,optsRep);
    dAIC_rep(i) = aic_gauss(sA,n,1) - aic_gauss(sB,n,1);

    bA_rep(i) = bA;
    bB_rep(i) = bB;

end

mean_bA = mean(bA_rep);
sd_bA = std(bA_rep,1);

mean_bB = mean(bB_rep);
sd_bB = std(bB_rep,1);

mean_dAIC = mean(dAIC_rep);
sd_dAIC = std(dAIC_rep,1);

count_dAIC_lt2 = sum(abs(dAIC_rep)<2);
count_dAIC_lt4 = sum(abs(dAIC_rep)<4);
count_dAIC_gt10 = sum(abs(dAIC_rep)>10);

count_A_preferred = sum(dAIC_rep<0);


fprintf('\nReplicate analysis (%d realizations)\n',nRep);

fprintf('High Transmission: beta_hat = %.4f +/- %.4f\n', ...
    mean_bA,sd_bA);

fprintf('Low Recovery: beta_hat = %.4f +/- %.4f\n', ...
    mean_bB,sd_bB);

fprintf('AIC_A - AIC_B: mean = %.2f, SD = %.2f\n', ...
    mean_dAIC,sd_dAIC);

fprintf('|dAIC| < 2: %d of %d\n',count_dAIC_lt2,nRep);
fprintf('|dAIC| < 4: %d of %d\n',count_dAIC_lt4,nRep);
fprintf('|dAIC| > 10: %d of %d\n',count_dAIC_gt10,nRep);

fprintf('Lower AIC under High Transmission: %d of %d\n', ...
    count_A_preferred,nRep);


fid = fopen(fullfile(resDir,'SyntheticReplicates.csv'),'w');

fprintf(fid,'replicate,beta_hat_A,beta_hat_B,AIC_A_minus_AIC_B\n');

for i = 1:nRep

    fprintf(fid,'%d,%.6f,%.6f,%.6f\n', ...
        i,bA_rep(i),bB_rep(i),dAIC_rep(i));

end

fclose(fid);


fid = fopen(fullfile(resDir,'SyntheticReplicateSummary.csv'),'w');

fprintf(fid,'statistic,value\n');
fprintf(fid,'nReplicates,%d\n',nRep);
fprintf(fid,'beta_A_mean,%.8f\n',mean_bA);
fprintf(fid,'beta_A_sd,%.8f\n',sd_bA);
fprintf(fid,'beta_B_mean,%.8f\n',mean_bB);
fprintf(fid,'beta_B_sd,%.8f\n',sd_bB);
fprintf(fid,'dAIC_mean,%.8f\n',mean_dAIC);
fprintf(fid,'dAIC_sd,%.8f\n',sd_dAIC);
fprintf(fid,'abs_dAIC_below_2,%d\n',count_dAIC_lt2);
fprintf(fid,'abs_dAIC_below_4,%d\n',count_dAIC_lt4);
fprintf(fid,'abs_dAIC_above_10,%d\n',count_dAIC_gt10);
fprintf(fid,'A_preferred,%d\n',count_A_preferred);

fclose(fid);

fprintf('\nReplicate results saved in:\n');
fprintf('%s\n',resDir);

%%
fig = figure('Color','w','Units','inches','Position',[1 1 6 4.5]);
hold on
plot(t_obs,inc_true,'-','Color',[0.25 0.25 0.25],'LineWidth',2.5);
scatter(t_obs,y_obs,32,'o','MarkerFaceColor','k','MarkerEdgeColor','w','LineWidth',0.5);
xlabel('Days','FontSize',13)
ylabel('New symptomatic cases per day','FontSize',13)
xlim([0 Tmax])
ylim([0 1.05*max([inc_true(:);y_obs(:)])])
box on; grid on
set(gca,'FontSize',14,'LineWidth',1.2,'TickDir','out','Layer','top')
legend({'True trajectory','Synthetic observations'},'Location','northwest','Box','off','FontSize',11)
exportgraphics(fig,fullfile(figDir,'Sup','SyntheticData.png'),'Resolution',600);

figure('Color','w','Units','inches','Position',[1 1 6 4.5]);
scatter(t_obs,y_obs,35,'k','filled');
hold on
plot(t_obs,inc_true,'k--','LineWidth',2);
plot(t_obs,inc_fit_A,'LineWidth',2,'Color',colA);
plot(t_obs,inc_fit_B,'LineWidth',2,'Color',colB);
xlabel('Days');
ylabel('New symptomatic cases per day');
legend('Synthetic observations','True trajectory', ...
       sprintf('High transmission (A), AIC=%.1f',AIC_A), ...
       sprintf('Low recovery (B), AIC=%.1f',AIC_B),'Location','northwest');
box on; grid on
set(gca,'FontSize',14)
%
exportgraphics(gcf,fullfile(figDir,'Figure2B_Synthetic_Fits.png'),'Resolution',600);
%exportgraphics(gcf,fullfile(figDir,'Figure2B_Synthetic_Fits.pdf'),'ContentType','vector');
%
header = {'Case','r_target','nu','beta_closed_form','beta_fit','R0','AIC'};
rows = {'Truth',                rstar,nu_T,beta_T,NaN,       R0_T,NaN;
        'High Transmission (A)',rstar,nu_A,beta_A,beta_hat_A,R0_A,AIC_A;
        'Low Recovery (B)',     rstar,nu_B,beta_B,beta_hat_B,R0_B,AIC_B};
write_table(fullfile(resDir,'SyntheticSummaryTable'),header,rows, ...
    ['Synthetic data summary. Closed form beta is the value predicted by the degeneracy curve at the stated nu and r target. ' ...
     'Fitted beta is recovered by fitting daily symptomatic incidence to the noisy synthetic data (n = 70). ' ...
     'AIC values are comparable only within this dataset.'],'tab:synthetic_summary');
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

function AIC = aic_gauss(SSE,n,k)
logLik = -0.5*n*(log(2*pi*SSE/n)+1);
AIC = 2*k - 2*logLik;
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

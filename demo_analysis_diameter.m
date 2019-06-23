% Random Walker segmentation (RaW)
% Quantification of axonal diamter
% author: Hong-Hsi Lee, 2018

clear

% Setup your directory, remember to change it if necessary
root = '.';

% Setup other directories
rootdata = fullfile(root,'data');
addpath(genpath(fullfile(root,'tools')));
addpath(genpath(fullfile(root,'lib')));
target = fullfile(root,'result'); mkdir(target);

%% Load segmentation
load(fullfile(target,'fiberfill.mat'))
load(fullfile(target,'proofread2_fiberlabel.mat'))
load(fullfile(target,'myelinsheath.mat'))

%% Diameter variation along each fiber
clear fb
as = analyzeseg();
vox = [100 100 100]*1e-3;  % voxel size, µm
sigma = 1;  % smoothing kernel width to smooth the skeleton, µm
diametermax = 10;  % The possible maximal diameter, µm
tic;
parfor i = 1:numel(I)
    fiber = fiberfill == I(i);
    myelin = myelins == I(i);
    fb(i) = as.fiberdiameter(fiber,myelin,sigma,diametermax,vox);
    i
end
toc;
save(fullfile(target,'diameter.mat'),'fb');

%% Simulate the time-dependence of diameter
load(fullfile(target,'diameter.mat'))
vox = [100 100 100]*1e-3;  % voxel size, µm
t = 1:100;              % diffusion time
D = 2;                  % intrinsic diffusivity
sigma = sqrt(D*t/2);    % smoothing kernel width
cutlength = 1;          % truncate 1 µm at both ends

as = analyzeseg();
dhist = zeros(numel(t),1);
dmri = zeros(numel(t),1);
parfor i = 1:numel(t)
    diameters = [];
    for j = 1:numel(fb)
        [diameter,zaxis] = as.smoothdiameter(fb(j).diameter,fb(j).zaxis,sigma(i),cutlength,vox);
        diameters = cat(1,diameters,diameter(:));
    end
    dhist(i) = mean(diameters);
    dmri(i) = as.diameternuman(diameters);
end
save(fullfile(target,'diameter_time_dependence.mat'),'t','dhist','dmri');

%% Plot diameter variation along each fiber
load(fullfile(target,'diameter.mat'));

vox = [100 100 100]*1e-3;  % voxel size, µm
t = [1 10 100];         % diffusion time
D = 2;                  % intrinsic diffusivity
sigma = sqrt(D*t/2);    % smoothing kernel width
cutlength = 1;          % truncate 1 µm at both ends

as = analyzeseg();
figure('unit','inch','position',[0 0 16 10])
hh = gobjects(3,1); lgdtext = strings(3,1);
for i = 1:numel(t)
    subplot(2,3,i); hold on
    diameters = [];
    for j = 1:numel(fb)
        [diameter,zaxis] = as.smoothdiameter(fb(j).diameter,fb(j).zaxis,sigma(i),cutlength,vox);
        h = plot(zaxis,diameter);
        diameters = cat(1,diameters,diameter(:));
    end
    set(gca,'xtick',0:10:30,'ytick',0:4)
    ax = ancestor(h(1),'axes');
    xrule = ax.XAxis; xrule.FontSize = 12;
    yrule = ax.YAxis; yrule.FontSize = 12;
    box on; pbaspect([1 1 1])
    xlim([0 30]); ylim([0 4])
    xlabel('$z_{\rm axon}$ ($\mu$m)','interpreter','latex','fontsize',20)
    ylabel('$2r$ ($\mu$m)','interpreter','latex','fontsize',20)
    title(sprintf('$t$ = %u ms',t(i)),'interpreter','latex','fontsize',20)
    
    subplot(2,3,4); hold on
    edges = linspace(0,4,50);
    centers = edges(1:end-1)/2 + edges(2:end)/2;
    N = histcounts(diameters,edges);
    N = N/sum(N)/mean(diff(edges));
    hh(i) = plot(centers,N);
    set(hh(i),'linewidth',2);
    lgdtext{i} = sprintf('$t$ = %u ms',t(i));
    set(gca,'xtick',0:4,'ytick',0:0.5:2)
    ax = ancestor(hh(1),'axes');
    xrule = ax.XAxis; xrule.FontSize = 12;
    yrule = ax.YAxis; yrule.FontSize = 12;
    box on; pbaspect([1 1 1]); grid on
    xlim([0 4]); ylim([0 2])
    xlabel('$2r$ ($\mu$m)','interpreter','latex','fontsize',20)
    ylabel('PDF ($\mu$m$^{-1}$)','interpreter','latex','fontsize',20)
end
legend(hh,lgdtext,'interpreter','latex','fontsize',20)

load(fullfile(target,'diameter_time_dependence.mat'));
subplot(2,3,5)
hold on
h1 = plot(t,dhist,'--k'); set(h1,'linewidth',2);
h2 = plot(t,dmri,'-k'); set(h2,'linewidth',2);
set(gca,'xtick',0:20:100,'ytick',0.9:0.1:1.6)
ax = ancestor(h1,'axes');
xrule = ax.XAxis; xrule.FontSize = 12;
yrule = ax.YAxis; yrule.FontSize = 12;
box on; pbaspect([1 1 1]); grid on
xlim([0 100]); ylim([0.9 1.6])
xlabel('$t$ (ms)','interpreter','latex','fontsize',20)
ylabel('inner diameter ($\mu$m)','interpreter','latex','fontsize',20)
legend([h1 h2],{'$2\langle r \rangle$','$2r_{\rm eff}$'},'interpreter','latex','fontsize',20)

mkdir(fullfile(target,'diameter'))
savefig(fullfile(target,'diameter','diameter_time_dependence.fig'))

%% Report size-related metrics in different definitions
load(fullfile(target,'diameter.mat'))
vox = [100 100 100]*1e-3;  % voxel size, µm
cutlength = 1;             % truncate 1 µm at both ends
cutnum = ceil(cutlength/vox(3));
diameter = [];
shortaxis = [];
longaxis = [];
incircle = [];
gratio = [];
outerdiameter = [];
for i = 1:numel(fb)
    diameter = cat(1,diameter,fb(i).diameter(cutnum+1:end-cutnum));
    shortaxis = cat(1,shortaxis,fb(i).shortaxis(cutnum+1:end-cutnum));
    longaxis = cat(1,longaxis,fb(i).longaxis(cutnum+1:end-cutnum));
    incircle = cat(1,incircle,fb(i).incircle(cutnum+1:end-cutnum));
    gratio = cat(1,gratio,fb(i).gratio(cutnum+1:end-cutnum));
    outerdiameter = cat(1,outerdiameter,fb(i).outerdiameter(cutnum+1:end-cutnum));
end
Iexclude = logical(isnan(gratio) + isinf(gratio) + (outerdiameter==0) + (diameter == outerdiameter) );
diameter(Iexclude) = [];
outerdiameter(Iexclude) = [];
gratio(Iexclude) = [];
eccentricity = sqrt(1-shortaxis(longaxis>0).^2./longaxis(longaxis>0).^2);

fprintf('Equivalent circle diameter = %.2f ± %.2f (mean ± SD).\n',mean(diameter),std(diameter));
fprintf('Equivalent circle diameter = %.2f ± %.2f (median ± IQR).\n',median(diameter),iqr(diameter));

fprintf('Short axis = %.2f ± %.2f (mean ± SD).\n',mean(shortaxis),std(shortaxis));
fprintf('Short axis = %.2f ± %.2f (median ± IQR).\n',median(shortaxis),iqr(shortaxis));

fprintf('Long axis = %.2f ± %.2f (mean ± SD).\n',mean(longaxis),std(longaxis));
fprintf('Long axis = %.2f ± %.2f (median ± IQR).\n',median(longaxis),iqr(longaxis));

fprintf('Inscribed circle diameter = %.2f ± %.2f (mean ± SD).\n',mean(incircle),std(incircle));
fprintf('Inscribed circle diameter = %.2f ± %.2f (median ± IQR).\n',median(incircle),iqr(incircle));

fprintf('Outer diameter = %.2f ± %.2f (mean ± SD).\n',mean(outerdiameter),std(outerdiameter));
fprintf('Outer diameter = %.2f ± %.2f (median ± IQR).\n',median(outerdiameter),iqr(outerdiameter));

fprintf('g-ratio = %.2f ± %.2f (mean ± SD).\n',mean(gratio),std(gratio));
fprintf('g-ratio = %.2f ± %.2f (median ± IQR).\n',median(gratio),iqr(gratio));

fprintf('Eccentricity = %.2f ± %.2f (mean ± SD).\n',mean(eccentricity),std(eccentricity));
fprintf('Eccentricity = %.2f ± %.2f (median ± IQR).\n',median(eccentricity),iqr(eccentricity));

% Fit g-ratio to linear log model
as = analyzeseg();
C = as.loglinearfit(diameter,outerdiameter);
di = 0.1:0.1:4;
gratio_theory = as.loglinearmodel(di,C);

%% Plot inner diameter, outer diameter, and g-ratio histogram
close all
figure('unit','inch','position',[0 0 10 10])
subplot(221);
edges = linspace(0,5,50);
centers = edges(1:end-1)/2 + edges(2:end)/2;
N = histcounts(outerdiameter,edges);
N = N/sum(N)/mean(diff(edges));
h = bar(centers,N,1);
xlim([0 4]); ylim([0 1.2]); pbaspect([1 1 1]);
set(gca,'xtick',0:5,'ytick',0:0.2:1.2)
ax = ancestor(h,'axes');
xrule = ax.XAxis; xrule.FontSize = 16;
yrule = ax.YAxis; yrule.FontSize = 16;
box on; grid on;
xlabel('outer diameter ($\mu$m)','interpreter','latex','fontsize',20);
ylabel('PDF ($\mu$m$^{-1}$)','interpreter','latex','fontsize',20)

subplot(222);
N = histcounts(diameter,edges);
N = N/sum(N)/mean(diff(edges));
h = bar(centers,N,1);
xlim([0 4]); ylim([0 1.2]); pbaspect([1 1 1]);
set(gca,'xtick',0:5,'ytick',0:0.2:1.2)
ax = ancestor(h,'axes');
xrule = ax.XAxis; xrule.FontSize = 16;
yrule = ax.YAxis; yrule.FontSize = 16;
box on; grid on;
xlabel('inner diameter ($\mu$m)','interpreter','latex','fontsize',20);
ylabel('PDF ($\mu$m$^{-1}$)','interpreter','latex','fontsize',20)

subplot(223);
edges = linspace(0,1,50);
centers = edges(1:end-1)/2 + edges(2:end)/2;
N = histcounts(gratio,edges);
N = N/sum(N)/mean(diff(edges));
h = bar(centers,N,1);
xlim([0 1]); ylim([0 5]); pbaspect([1 1 1]);
set(gca,'xtick',0:0.2:1,'ytick',0:6)
ax = ancestor(h,'axes');
xrule = ax.XAxis; xrule.FontSize = 16;
yrule = ax.YAxis; yrule.FontSize = 16;
box on; grid on;
xlabel('g-ratio','interpreter','latex','fontsize',20);
ylabel('PDF','interpreter','latex','fontsize',20)

subplot(224);
hist2(diameter,gratio,100,[0 5],[0 1],'k');
box on; grid on; xlim([0 4]); ylim([0 1])
pbaspect([1 1 1])
ax = ancestor(gca,'axes');
xrule = ax.XAxis; xrule.FontSize = 16;
yrule = ax.YAxis; yrule.FontSize = 16;
set(gca,'xtick',0:5,'ytick',0:0.2:1)
xlabel('inner diameter ($\mu$m)','interpreter','latex','fontsize',20);
ylabel('g-ratio','interpreter','latex','fontsize',20);

hold on
h = plot(di,gratio_theory,'-r'); set(h,'linewidth',1.5)
legend(h,{'log-linear'},'location','southeast','fontsize',22,'interpreter','latex')
savefig(fullfile(target,'diameter','diameter_gratio_histogram.fig'))

as = analyzeseg();
[gmri,sigmagmri] = as.gratiomri(gratio,outerdiameter);
fprintf('g_MRI = %.2f ± %.2f\n',gmri,sigmagmri);

%% Plot histogram of diameters in different definitions
lgdtext = {{'Equivalent circle';'diameter ($\mu$m)'},{'Short axis';'length ($\mu$m)'},...
    {'Long axis';'length ($\mu$m)'},{'Inscribed circle';'diameter ($\mu$m)'}};
Nbin = 50;
edges = linspace(0,5,Nbin);
centers = edges(1:end-1)/2 + edges(2:end)/2;
pd = struct([]);
pd(1).diameter = diameter;
pd(2).diameter = shortaxis;
pd(3).diameter = longaxis;
pd(4).diameter = incircle;
for i = 1:4
    di = pd(i).diameter;
    N = histcounts(di,edges);
    pd(i).N = N/sum(N)/mean(diff(edges));
    pd(i).gamma = fitdist(di,'Gamma');
    pd(i).gev = fitdist(di,'GeneralizedExtremeValue');
    pd(i).mean = mean(di);
    pd(i).median = median(di);
end

cmap = colormap('lines'); close gcf
edges2 = linspace(0,5,100);
figure('unit','inch','position',[0 0 20 10])
for i = 1:4
    subplot(2,4,i); hold on
    hbar = bar(centers,pd(i).N,1); set(hbar,'linewidth',1.5,'facecolor','k','facealpha',0.2);
    
    ygev = pdf(pd(i).gev,edges2);
    ygamma = pdf(pd(i).gamma,edges2);
    hgev = plot(edges2,ygev,'-r'); set(hgev,'linewidth',2.5,'color',cmap(1,:))
    hgamma = plot(edges2,ygamma,'-b'); set(hgamma,'linewidth',2.5,'color',cmap(7,:))
    
    legend([hgev,hgamma], {'GEV','Gamma'},'interpreter','latex','fontsize',20,'linewidth',1.5)
    
    box on; grid on; xlim([0 5]); ylim([0 1.6]); pbaspect([1 1 1])
    ax = ancestor(gca,'axes');
    xrule = ax.XAxis; xrule.FontSize = 16;
    yrule = ax.YAxis; yrule.FontSize = 16;
    set(gca,'xtick',0:5,'ytick',0:0.5:1.5,'linewidth',1.5)
    xlabel(lgdtext{i},'interpreter','latex','fontsize',20)
    ylabel('PDF ($\mu$m$^{-1}$)','interpreter','latex','fontsize',20)
end

for i = 1:4
    subplot(2,4,i+4); hold on
    hbar = plot(centers,pd(i).N,'.k'); set(hbar,'markersize',20);
    
    ygev = pdf(pd(i).gev,edges2);
    ygamma = pdf(pd(i).gamma,edges2);
    hgev = plot(edges2,ygev,'-r'); set(hgev,'linewidth',2.5,'color',cmap(1,:))
    hgamma = plot(edges2,ygamma,'-b'); set(hgamma,'linewidth',2.5,'color',cmap(7,:))
        
    legend([hgev,hgamma], {'GEV','Gamma'},'interpreter','latex','fontsize',20,'linewidth',1.5)
    box on; grid on; xlim([0 5]); ylim([1e-4 1e1]); pbaspect([1 1 1])
    ax = ancestor(gca,'axes');
    xrule = ax.XAxis; xrule.FontSize = 16;
    yrule = ax.YAxis; yrule.FontSize = 16;
    set(gca,'yscale','log','xtick',0:5,'ytick',10.^(-4:2:2),'linewidth',1.5)
    xlabel(lgdtext{i},'interpreter','latex','fontsize',20)
    ylabel('PDF ($\mu$m$^{-1}$)','interpreter','latex','fontsize',20)
end
savefig(fullfile(target,'diameter','diameter_definitions.fig'))


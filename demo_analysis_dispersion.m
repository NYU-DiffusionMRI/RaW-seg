% Random Walker segmentation (RaW)
% Quantification of orientation dispersion
% author: Hong-Hsi Lee, 2018

clear

% Setup your directory, remember to change it if necessary
root = '.';

% Setup other directories
rootdata = fullfile(root,'data');
addpath(genpath(fullfile(root,'tools')));
addpath(genpath(fullfile(root,'lib')));
target = fullfile(root,'result'); mkdir(target);

%% Fiber orientation distribution (FOD)
mkdir(fullfile(target,'FOD'))
load(fullfile(target,'fibers.mat'))
load(fullfile(target,'proofread2_fiberlabel.mat'))

% Calculate fiber skeleton
vox = [24,24,100]*1e-3;  % voxel size of the EM image, µm
[nx,ny,~] = size(fibers);
skeleton = struct([]);
as = analyzeseg();
parfor i = 1:numel(I)
    fiberi = fibers == I(i);
    cm = as.centermass(fiberi,vox);
    skeleton(i).cm = cm;
end
save(fullfile(target,'skeleton.mat'),'skeleton')

%% Plot time-dependence of the fiber skeleton
load(fullfile(target,'skeleton.mat'))
vox = [24,24,100]*1e-3;  % voxel size of the EM image, µm
t = [1,10,100];          % diffusion time
D = 2;                   % intrinsic diffusivity
sigma = sqrt(D*t/2);     % smoothing kernel width
as = analyzeseg();
fod = struct([]);
for i = 1:numel(t)
    % Plot 3D skeletons
    figure('unit','inch','position',[0 0 10 10])
    hold on; grid on; axis equal off; view(3);
    tangent = [];
    for j = 1:numel(skeleton)
       [tg,c] = as.smoothtangent(skeleton(j).cm,sigma(i),vox(3));
       tangent = cat(1,tangent,tg);
       h = plot3(c(:,1),c(:,2),c(:,3));
       set(h,'linewidth',1);
    end
    pause(1)
    savefig(fullfile(target,'FOD',sprintf('3dline_view3_t=%ums.fig',t(i))))
    
    pause(1)
    % View 3D skeletons from the top
    view([90 -90])
    savefig(fullfile(target,'FOD',sprintf('3dline_view2_t=%ums.fig',t(i))))
    
    % Plot FOD projected on a triangulated sphere
    figure('unit','inch','position',[0 0 10 10]); hold on
    fod(i).diffusiontime = t(i);
    fod(i).tangent = tangent;
    fod(i).fod = as.fodsphere([tangent;-tangent],'density','low','range',[0 10]);
    [U,~,~] = svd(tangent'*tangent/size(tangent,1));
    v = U(:,3); v = sign(v(1))*normc([v(1) v(2) 0].');
    fod(i).viewangle = v;
    h = quiver3(v(1)*1.1,v(2)*1.1,0,-v(1)*0.1,-v(2)*0.1,0,'k');
    set(h,'linewidth',3,'maxheadsize',10)
    savefig(fullfile(target,'FOD',sprintf('odf_hist_t=%ums.fig',t(i))))
end

% Decompose FOD by using spherical harmonics and plot the 3D glyph
nord = 10;  % order of spherical harmonics
as = analyzeseg();
for i = 1:numel(t)
    figure('unit','inch','position',[0 0 10 10])
    as.fodsht(fod(i).fod,nord,'smooth','false','density','low','range',[0 10],'colorbar','off','glyph','on');
    xlim([-1 1]*1e-3); ylim([-1 1]*1e-3); zlim([-2 2]*1e-3)
    view(fod(i).viewangle)
    pause(1)
    savefig(fullfile(target,'FOD',sprintf('glyph_sht_t=%ums.fig',t(i))))
end
%%
mkdir(fullfile(target,'dispersionangle'))
% Dispersion angle calculated by projecting fiber segments to 2d planes
figure('unit','inch','position',[0 0 4.5 4.5]); hold on
nphi = 30;  % # azimuthal angles
t = [1,10,100];  % diffusion time
as = analyzeseg();
h = gobjects(3,1); lgdtext = strings(3,1);
for i = 1:numel(t)
    [dispang,phi] = as.dispersion2d(fod(i).tangent,nphi);
    fod(i).dispang2d = dispang;
    fod(i).phi2d = phi;
    h(i) = plot(phi,dispang);
    lgdtext{i} = sprintf('$t$ = %u ms',t(i));
    set(gca,'xtick',-180:60:180)
end
ax = ancestor(gca,'axes');
xrule = ax.XAxis; xrule.FontSize = 12;
yrule = ax.YAxis; yrule.FontSize = 12;
set(h,'linewidth',2)
xlim([-180 180]); ylim([0 35]);
xlabel('$\phi(^\circ)$','interpreter','latex','fontsize',16)
ylabel('$\theta_{2d}(\phi)(^\circ)$','interpreter','latex','fontsize',16)
pbaspect([1 1 1]); box on; grid on
legend(h,lgdtext,'interpreter','latex','fontsize',12,'location','southwest')
savefig(fullfile(target,'dispersionangle','theta2D_vs_phi_3times.fig'))

% Dispersion angle of 3d fiber segments
figure('unit','inch','position',[0 0 4.5 4.5]); hold on
nphi = 30;  % # azimuthal angles
t = [1,10,100];  % diffusion time
as = analyzeseg();
h = gobjects(3,1); lgdtext = strings(3,1);
for i = 1:numel(t)
    [dispang,phi] = as.dispersion3d(fod(i).tangent,nphi);
    fod(i).dispang3d = dispang;
    fod(i).phi3d = phi;
    h(i) = plot(phi,dispang);
    lgdtext{i} = sprintf('$t$ = %u ms',t(i));
    set(gca,'xtick',-180:60:180)
end
ax = ancestor(gca,'axes');
xrule = ax.XAxis; xrule.FontSize = 12;
yrule = ax.YAxis; yrule.FontSize = 12;
set(h,'linewidth',2)
xlim([-180 180]); ylim([0 35]);
xlabel('$\phi(^\circ)$','interpreter','latex','fontsize',16)
ylabel('$\theta_{\rm eff}(\phi)(^\circ)$','interpreter','latex','fontsize',16)
pbaspect([1 1 1]); box on; grid on
legend(h,lgdtext,'interpreter','latex','fontsize',12,'location','southwest')
savefig(fullfile(target,'dispersionangle','thetaMR_vs_phi_3times.fig'))

%%
% Calculate time-dependent dispersion angles and rotational invariants
vox = [24,24,100]*1e-3; % voxel size of the EM image, µm
t = 1:100;              % diffusion time
D = 2;                  % intrinsic diffusivity
sigma = sqrt(D*t/2);    % smoothing kernel width
nord = 10;
dispang2d = zeros(numel(t),1);
dispang3d = zeros(numel(t),1);
pl = zeros(floor(nord/2)+1,numel(t));
as = analyzeseg();
tic;
parfor i = 1:numel(t)
    close all
    tangent = [];
    for j = 1:numel(skeleton)
       [tg,~] = as.smoothtangent(skeleton(j).cm,sigma(i),vox(3));
       tangent = cat(1,tangent,tg);
    end
    dispang2d(i) = as.dispersionangle2d(tangent);
    dispang3d(i) = as.dispersionangle3d(tangent);
    
    fodi = as.fodsphere([tangent;-tangent],'density','low','range',[0 10]);
    pl(:,i) = as.rotinv(fodi,nord);
    i
end
dispangp2 = as.dispersionanglep2(pl(2,:));
toc;
save(fullfile(target,'FOD','dispang_rotinv_time_dependence.mat'),...
    'dispang2d','dispang3d','dispangp2','pl','t');

%%
load(fullfile(target,'FOD','dispang_rotinv_time_dependence.mat'))
% Plot time-dependent rotational invariants
figure('unit','inch','position',[0 0 4.5 4.5]); hold on
as = analyzeseg();
h = gobjects(3,1); lgdtext = strings(3,1);
shp = {'-k','--k',':k'};
for i = 1:3
    h(i) = plot(t,pl(i+1,:),shp{i});
    lgdtext{i} = sprintf('$p_%u$',2*i);
    set(gca,'xtick',0:20:100,'ytick',0:0.2:1)
end
ax = ancestor(gca,'axes');
xrule = ax.XAxis; xrule.FontSize = 12;
yrule = ax.YAxis; yrule.FontSize = 12;
set(h,'linewidth',2)
xlim([0 100]); ylim([0 1]);
xlabel('$t$ (ms)','interpreter','latex','fontsize',16)
ylabel('$p_l$','interpreter','latex','fontsize',16)
pbaspect([1 1 1]); box on; grid on
legend(h,lgdtext,'interpreter','latex','fontsize',16,'location','southeast')
savefig(fullfile(target,'dispersionangle','rotational_invariant_vs_diffusiontime.fig'))

% Plot time-dependent dispersion angles
figure('unit','inch','position',[0 0 4.5 4.5]); hold on
h = gobjects(3,1);
shp = {'-k',':k','--k'};
as = analyzeseg();
dispangs = [dispang3d(:), dispangp2(:), dispang2d(:)];
lgdtext = {'$\theta_{\rm eff}$','$\theta_{p_2}$','$\theta_{2d}$'};
for i = 1:3
    h(i) = plot(t(:),dispangs(:,i),shp{i});
    set(gca,'xtick',0:20:100,'ytick',10:30)
end
ax = ancestor(gca,'axes');
xrule = ax.XAxis; xrule.FontSize = 12;
yrule = ax.YAxis; yrule.FontSize = 12;
set(h,'linewidth',2)
xlim([0 100]); ylim([16 24]);
xlabel('$t$ (ms)','interpreter','latex','fontsize',16)
ylabel('angles ($^\circ$)','interpreter','latex','fontsize',16)
pbaspect([1 1 1]); box on; grid on
legend(h,lgdtext,'interpreter','latex','fontsize',16,'location','east')
savefig(fullfile(target,'dispersionangle','theta_eff_p2_2d_vs_diffusiontime.fig'))

%%
% Plot time-dependent rotational invariants wrt order l
figure('unit','inch','position',[0 0 4.5 4.5]); hold on
as = analyzeseg();
h = gobjects(3,1); lgdtext = strings(3,1);
shp = {'-k','--k',':k'};
j = 0; l = (0:2:10);
for i = [1,10,100]
    j = j+1;
    h(j) = plot(l(:),pl(:,i),'.-');
    lgdtext{j} = sprintf('$t$ = %u ms',t(i));
    set(gca,'xtick',0:2:10,'ytick',[0.1 0.2:0.2:1],'yscale','log')
end
ax = ancestor(gca,'axes');
xrule = ax.XAxis; xrule.FontSize = 12;
yrule = ax.YAxis; yrule.FontSize = 12;
set(h,'linewidth',2,'markersize',20)
xlim([0 10]); ylim([0.1 1]);
xlabel('$l$','interpreter','latex','fontsize',16)
ylabel('$p_l$','interpreter','latex','fontsize',16)
pbaspect([1 1 1]); box on; grid on
legend(h,lgdtext,'interpreter','latex','fontsize',16,'location','southwest')
savefig(fullfile(target,'dispersionangle','rotational_invariant_3times.fig'))

% Fit rotational invariants to the Poisson model
[lambda,C] = as.rotinvpoisson(pl(2:end,:),l(2:end));
figure('unit','inch','position',[0 0 4.5 4.5]); hold on
h = gobjects(2,1);
h(1) = plot(t,C,'--k');
h(2) = plot(t,lambda,'-k');
lgdtext = {'$C$','$\lambda$'};
ax = ancestor(gca,'axes');
xrule = ax.XAxis; xrule.FontSize = 12;
yrule = ax.YAxis; yrule.FontSize = 12;
set(h,'linewidth',2)
xlim([0 100]); ylim([0.7 1.3]);
xlabel('$t$ (ms)','interpreter','latex','fontsize',16)
pbaspect([1 1 1]); box on; grid on
legend(h,lgdtext,'interpreter','latex','fontsize',16,'location','east')
savefig(fullfile(target,'dispersionangle','power-law-dispersion-parameters.fig'))

%%
% Calculate orientation dispersion index, defined by Bingham distribution
load(fullfile(target,'skeleton.mat'))
vox = [24,24,100]*1e-3; % voxel size of the EM image, µm
t = 1;                  % diffusion time
D = 2;                  % intrinsic diffusivity
sigma = sqrt(D*t/2);    % smoothing kernel width
tangent = [];
for i = 1:numel(skeleton)
   [tg,~] = as.smoothtangent(skeleton(i).cm,sigma,vox(3));
   tangent = cat(1,tangent,tg);
end
as = analyzeseg();
[kappa,odi] = as.bingham(tangent);
fprintf('By fitting FOD to Bingham distribution, ODI_1 = %.3f, ODI_2 = %.3f.\n',odi(1),odi(2))


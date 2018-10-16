classdef analyzeseg < handle
    properties (Constant = true, Access = protected)
        
    end
    
    properties (GetAccess = public, SetAccess = protected)
    end
    
    properties (GetAccess = private, SetAccess = protected)
        
    end
    
    methods (Access = public)
        function this = analyzeseg(varargin)
            
        end
        
        function [tangent,cms,zaxis] = smoothtangent(this,cm,sigma,dz)
            u = normr(cm(end,:)-cm(1,:));
            cmr = this.rotateaxis(cm,u,mean(cm,1),'forward');
            nzq = round(norm(cmr(end,:)-cmr(1,:))/dz);
            zq = linspace(cmr(1,3),cmr(end,3),nzq);
            zq = zq(:);
            xq = interp1(cmr(:,3),cmr(:,1),zq,'spline');
            yq = interp1(cmr(:,3),cmr(:,2),zq,'spline');
            xq = imgaussfilt(xq,sigma/mean(diff(zq)));
            yq = imgaussfilt(yq,sigma/mean(diff(zq)));
            cms = [xq,yq,zq];
            cms = this.rotateaxis(cms,u,mean(cm,1),'backward');
            tangent = normr(diff(cms,1,1));
            zaxis = zq-zq(1);
        end
        
        function dispang = dispersionangle2d(this,tangent)
%DISPERSIONANGLE2D    Dispersion angle of projected fiber segments.
%   DISPERSIONANGLE2D(tangent) returns dispersion angle of projected fiber
%   segments on the 2D planes parallel to the main direction of unit vector
%   tangent (N x 3).
%
% author: Hong-Hsi Lee, 2018
            dispang = this.dispersion2d(tangent,360);
            dispang = rms(dispang);
        end
        
        function BWalign = bwalign(this,BW,center,u,R)
%BWALIGN    Align binary map along a given axis.
%   BWALIGN(BW,cm,u,R) returns the aligned binary map along the direction u,
%   with a boundary box in the width of R. The input BW is a 3D binary
%   map, and u and R are 3-element vectors.

            list = find(BW);
            [I,J,K] = ind2sub(size(BW),list);
            if isempty(center)
                center = [mean(I) mean(J) mean(K)];
            end
            rx = floor(R(1)/2); ry = floor(R(2)/2); rz = floor(R(3)/2);
            [~,M] = this.rotateaxis([],u,center,'backward');

            rx = (center(1)-rx):(center(1)-rx+R(1)-1+1e-10);
            ry = (center(2)-ry):(center(2)-ry+R(2)-1+1e-10);
            rz = (center(3)-rz):(center(3)-rz+R(3)-1+1e-10);
            BWalign = affineTrans(single(BW), M, rx, ry, rz,'lanzcos3') > 0.5;
        end
        
        function fb = fiberdiameter(this,fiber,myelin,sigma,diametermax,vox)
            fov = ceil(diametermax/vox(1));
            
            cm = this.centermass(fiber,vox);
            [tangent,cm,zaxis] = this.smoothtangent(cm,sigma,vox(3));
            tangent = [tangent;tangent(end,:)];
            fb.centermass = cm;
            fb.tangent = tangent;
            fb.zaxis = zaxis;

            nslice = size(cm,1);
            diameter = zeros(nslice,1);
            shortaxis = zeros(nslice,1);
            longaxis = zeros(nslice,1);
            incircle = zeros(nslice,1);
            outerdiameter = zeros(nslice,1);
            for i = 1:nslice
                fiberi = this.bwalign(fiber,cm(i,:)./vox,tangent(i,:),[fov fov 1]);
                if nnz(fiberi)==0, continue; end
                diameter(i) = sqrt(sum(fiberi(:))/pi)*2;
                stats = regionprops(fiberi,'minoraxislength','majoraxislength','area');
                [~,I] = max([stats.Area]);
                shortaxis(i) = stats(I).MinorAxisLength;
                longaxis(i) = stats(I).MajorAxisLength;
                fiberi = imresize(fiberi,10,'nearest');
                D = bwdist(~fiberi);
                incircle(i) = max(D(:))/10*2;

                myelini = this.bwalign(logical(myelin+fiber),cm(i,:)./vox,tangent(i,:),[fov fov 1]);
                outerdiameter(i) = sqrt(sum(myelini(:))/pi)*2;
            end
            fb.diameter = diameter*vox(1);
            fb.shortaxis = shortaxis*vox(1);
            fb.longaxis = longaxis*vox(1);
            fb.incircle = incircle*vox(1);
            fb.outerdiameter = outerdiameter*vox(1);
            fb.gratio = diameter./outerdiameter;
        end
    end
    
    methods (Static)
        function cm = centermass(fiber,vox)
            [nx,ny,~] = size(fiber);
            zlist = find(sum(sum(fiber,1),2));
            cm = zeros(numel(zlist),3);
            for i = 1:numel(zlist)
                fiberi = fiber(:,:,zlist(i));
                [I,J] = ind2sub([nx,ny],find(fiberi));
                cm(i,:) = [mean(I)*vox(1),mean(J)*vox(2),zlist(i)*vox(3)];
            end
        end
        
        function [xr,M] = rotateaxis(x,u,cm,direction)
            if nargin < 4, direction = 'forward'; end
            v = -[-u(2), u(1), 0];
            s = sqrt(sum(v.^2));
            c = u(3);
            V = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
            
            R = eye(4);
            R(1:3, 1:3) = eye(3) + V + V*V * (1-c)/s^2;
            if strcmpi(direction,'backward'), R = R'; end
            
            T1 = eye(4); T1(4, 1:3) = cm - [1 1 1]; T1 = T1';
            T2 = eye(4); T2(4, 1:3) = -cm + [1 1 1]; T2 = T2';
            
            M = T1*R*T2;
            if ~isempty(x)
                x = [x.'; ones(1,size(x,1))];
                xr = M*x;
                xr = xr(1:3,:).';
            else
                xr = [];
            end
        end
        
        function fod = fodsphere(tangent,varargin)
%FODSPHERE    Plot fiber orientation distribution on a triangulated sphere.
%   FODSPHERE(tangent,options) plots fiber orientation distribution based 
%   on unit vectors of fiber tangent (N x 3) over a triangulated sphere.
%   The unit of the histogram is /Steradian.
%   options can be:
%   'Smooth'     -smooth the face color with linear gradient ('on', 'off'),
%                 the default is 'off'
% 
%   'Density'    -the density of triangles over the sphere ('low', 'medium'
%                 , 'high'), the default is 'medium'
% 
%   'Range'      -the numeric range of histogram, a two-element array, the
%                 default is []
% 
%   'Colorbar'   -colorbar of histogram ('on', 'off'), the default is 'off'
% 
% author: Hong-Hsi Lee, 2017
            smoothflag = 'off';
            densityflag = 'medium';
            rangeflag = [];
            colorbarflag = 'off';
            if nargin > 1
                for i = 1:nargin-1
                    if strcmpi('Smooth',varargin{i})
                        smoothflag = varargin{i+1};
                    end
                    if strcmpi('Density',varargin{i})
                        densityflag = varargin{i+1};
                    end
                    if strcmpi('Range',varargin{i})
                        rangeflag = varargin{i+1};
                    end
                    if strcmpi('Colorbar',varargin{i})
                        colorbarflag = varargin{i+1};
                    end
                end
            end
            
            str = which('analyzeseg.m');
            [root,~,~] = fileparts(str);
            if strcmpi('low',densityflag)
                load(fullfile(root,'TR_low.mat'))
            elseif strcmpi('medium',densityflag)
                load(fullfile(root,'TR_medium.mat'))
            elseif strcmpi('high',densityflag)
                load(fullfile(root,'TR_high.mat'))
            end
            pts = TR.Points;
            K = TR.ConnectivityList;

            vi = nearestNeighbor(TR,tangent);
            fod = zeros(size(K,1),1);
            for i = 1:size(tangent,1)
                ti = vertexAttachments(TR,vi(i));
                N = neighbors(TR,ti{:}.');
                N = unique([N(:); ti{:}.']);
                B = cartesianToBarycentric(TR,N,repmat(tangent(i,:),[numel(N),1]));
                I = find( ( B(:,1)>=0 ) .* ( B(:,2)>=0 ) .* ( B(:,3)>=0 ) );
                if isempty(I)
                    continue
                end
                I = N(I(1));
                fod(I) = fod(I) +1;
            end
            solang = triarea/sum(triarea)*4*pi;
            normfod = sum(fod.*solang);
            fod = fod ./ normfod;
            
            if isempty(rangeflag)
                if strcmpi('on',smoothflag)
                    fodv = zeros(size(pts,1),1);
                    for i = 1:size(pts,1)
                        I = find(sum(K==i,2));
                        fodv(i) = sum(fod(I).*triarea(I))/sum(triarea(I));
                    end
                    h = patch('faces',K,'vertices',pts,'facevertexcdata',fodv,'facecolor','interp');
                    fod = fodv;
                elseif strcmpi('off',smoothflag)
                    h = patch('faces',K,'vertices',pts,'facevertexcdata',fod,'facecolor','flat');
                end
                set(h,'edgealpha',0); axis equal off
                colormap hot
                if strcmpi('on',colorbarflag)
                    colorbar
                end
            else
                cmap = colormap('hot');
                if strcmpi('on',smoothflag)
                    fodv = zeros(size(pts,1),1);
                    for i = 1:size(pts,1)
                        I = find(sum(K==i,2));
                        fodv(i) = sum(fod(I).*triarea(I))/sum(triarea(I));
                    end
                    fodv = ceil( (fodv-min(rangeflag))/range(rangeflag) * 64);
                    I = find(fodv<=0);
                    fodv(fodv<=0) = 64; fodv(fodv>64) = 64;
                    col = cmap(round(fodv),:);
                    col(I,:) = repmat([0.4 0.4 0.4],[numel(I),1]);
                    h = patch('faces',K,'vertices',pts,'facevertexcdata',col,'facecolor','interp');
                    fod = fodv;
                elseif strcmpi('off',smoothflag)
                    ffod = ceil( (fod-min(rangeflag))/range(rangeflag) * 64);
                    I = find(ffod<=0); 
                    ffod(ffod<=0) = 64; ffod(ffod>64) = 64;
                    col = cmap(round(ffod),:);
                    col(I,:) = repmat([0.4 0.4 0.4],[numel(I),1]);
                    h = patch('faces',K,'vertices',pts,'facevertexcdata',col,'facecolor','flat');
                end
                set(h,'edgealpha',0); axis equal off
                if strcmpi('on',colorbarflag)
                    colorbar; caxis(rangeflag)
                end
            end
        end
        
        function F = fodsht(fod,nord,varargin)
%FODSHT    Plot FOD by using spherical harmonics.
%   FODSHT(fod,nord,options) plots fiber orientation distribution fod
%   with spherical harmonics up to the order of nord. The unit of the
%   histogram is /Steradian.
%   options can be:
%   'Smooth'     -whether the FOD has been smoothed ('true', 'false'), the
%                 default is 'false'
%
%   'Density'    -the density of triangles over the sphere ('low', 'medium'
%                 , 'high'), the default is 'medium'
%
%   'Range'      -the numeric range of histogram, a two-element array, the
%                 default is []
%
%   'Colorbar'   -colorbar of histogram ('on', 'off'), the default is 'off'
%
%   'Glyph'      -plot 3D glyph ('on', 'off'). 'Range' and 'Colorbar' are
%                 not applicable when 'Glyph' is 'on'. The default is 'off'
%
% author: Hong-Hsi Lee, 2017
            smoothflag = 'false';
            densityflag = 'medium';
            rangeflag = [];
            colorbarflag = 'off';
            glyphflag = 'off';
            if nargin > 2
                for i = 1:nargin-2
                    if strcmpi('Smooth', varargin{i})
                        smoothflag = varargin{i+1};
                    end
                    if strcmpi('Density', varargin{i})
                        densityflag = varargin{i+1};
                    end
                    if strcmpi('Range', varargin{i})
                        rangeflag = varargin{i+1};
                    end
                    if strcmpi('Colorbar', varargin{i})
                        colorbarflag = varargin{i+1};
                    end
                    if strcmpi('Glyph', varargin{i})
                        glyphflag = varargin{i+1};
                    end
                end
            end

            str = which('analyzeseg.m');
            [root,~,~] = fileparts(str);
            if strcmpi('low',densityflag)
                load(fullfile(root,'TR_low.mat'))
            elseif strcmpi('medium',densityflag)
                load(fullfile(root,'TR_medium.mat'))
            elseif strcmpi('high',densityflag)
                load(fullfile(root,'TR_high.mat'))
            end

            if strcmpi('true',smoothflag)
                pts = TR.Points;
                K = TR.ConnectivityList;
                wts = ones(size(pts,1),1) *4*pi/size(pts,1);
            elseif strcmpi('false',smoothflag)
                pts = TRCC.Points;
                K = TRCC.ConnectivityList;
                wts = triarea/sum(triarea)* 4*pi;
            end

            [phi,~] = cart2pol(pts(:,1),pts(:,2));
            theta = acos(pts(:,3));
            dirs = [phi theta];
            F_N = directSHT(nord, fod, dirs, 'complex', wts);

            if strcmpi('on',glyphflag)
                    [X,Y,Z] = sphere(50);
                    [phi,~] = cart2pol(X(:),Y(:));
                    theta = acos(Z(:));
                    dirs2 = [phi theta];
                    F = inverseSHT(F_N, dirs2, 'complex');
                    F = abs(F); F = F/sum(F);

                    rho = reshape(F,size(X));
                    C = abs(cat(3,X,Y,Z));
                    X = rho.*X; Y = rho.*Y; Z = rho.*Z;
                    h = surf(X,Y,Z,C);
                    set(h,'edgealpha',1); axis equal off
            else
                F = inverseSHT(F_N, dirs, 'real');
                if isempty(rangeflag)
                    K = TRCC.ConnectivityList;
                    h = patch('faces',K,'vertices',pts,'facevertexcdata',F,'facecolor','interp');
                    set(h,'edgealpha',0); axis equal off; colormap hot
                    if strcmpi('on',colorbarflag), colorbar; end
                else
                    cmap = colormap('hot');
                    F = ceil( (F-min(rangeflag))/range(rangeflag) * 64);
                    I = F<=0;
                    F(F<=0) = 1; F(F>64) = 64;
                    col = cmap(round(F),:);
                    col(I,:) = zeros(sum(I),3);
                    h = patch('faces',K,'vertices',pts,'facevertexcdata',col,'facecolor','interp');
                    set(h,'edgealpha',0); axis equal off
                    if strcmpi('on',colorbarflag), colorbar; caxis(rangeflag); end
                end
            end
        end
        
        function [dispang,phi,U] = dispersion2d(tangent,nphi)
%DISPERSION2D    Calculate dispersion angles by projecting FOD to 2D plane. 
%   DISPERSION2D(tangent,nphi) project fiber orientation distribution based
%   on unit vector tangent (N x 3) to azimuthal angle phi in nphi
%   directions.
%
%   [dispang,phi,U] = DISPERSION2D(tangent,nphi) returns dispersion angle's
%   standard deviation dispang with respect to phi, and the PCA axis U.
%
% author: Hong-Hsi Lee, 2017
        
            % Flip tangent to upper-half sphere
            tangent = tangent.*repmat(sign(tangent(:,3)),[1,3]);
            % Define mean direction
            [U,~,~] = svd(tangent'*tangent/size(tangent,1));
            % Align tangent to the main direction
            Lx = [0 0 0; 0 0 -1; 0 1 0]; Ly = [0 0 1; 0 0 0; -1 0 0]; Lz = [0 -1 0; 1 0 0; 0 0 0];
            uL = U(2,1)*Lx - U(1,1)*Ly; uL = uL/sqrt(U(1,1)^2+U(2,1)^2);
            thetar = acos(U(3,1));
            tangentr = expm(thetar*uL)*tangent.'; tangentr = tangentr.';
            phi = linspace(-180,180,nphi+1);
            % Calculate projecting angle
            dispang = zeros(nphi,1);
            for i = 1:nphi
                phii = phi(i);
                Rz = [cosd(phii) -sind(phii) 0; sind(phii) cosd(phii) 0; 0 0 1];
                tangentri = tangentr*Rz;
                thetap = acosd(tangentri(:,3)./sqrt(1-tangentri(:,2).^2)).*sign(tangentri(:,1));
                dispang(i) = std(thetap);
            end
            phi = phi(1:end-1);
        end
        
        function [dispang,phi,U] = dispersion3d(tangent,nphi)
%DISPERSION3D    Calculate dispersion angles in 3D by using rms of cosines.
%   DISPERSION3D(tangent,nphi) returns the histogram of unit vector tangent
%   (N x 3) wrt the azimuthal angle into nphi bins.
%
%   [dispang,phi,U] = DISPERSION3D(tangent,nbin) returns the
%   dispersion angle dispang with respect to the azimuthal angle phi.
%   Tangents are aligned to the principle axis U(:,1). The dispersion
%   angle is defined by acosd(rms(cos)) within each bin.
%
% author: Hong-Hsi Lee, 2017

            % Flip tangent to upper-half sphere
            tangent = tangent.*repmat(sign(tangent(:,3)),[1,3]);

            % Mean direction
            [U,~,~] = svd(tangent.'*tangent/size(tangent,1));

            % Align tangent to the main direction
            Lx = [0 0 0; 0 0 -1; 0 1 0]; Ly = [0 0 1; 0 0 0; -1 0 0]; Lz = [0 -1 0; 1 0 0; 0 0 0];
            uL = U(2,1)*Lx - U(1,1)*Ly; uL = uL/sqrt(U(1,1)^2+U(2,1)^2);
            thetar = acos(U(3,1));
            tangentr = expm(thetar*uL)*tangent.'; tangentr = tangentr.';

            % Bin tangents into different azimuthal angle phi
            u = tangentr(:,1); v = tangentr(:,2); w = tangentr(:,3);
            [phii,~] = cart2pol(u,v); phii = phii*180/pi;
            edges = linspace(-180,180,nphi+1).';
            I = zeros(size(tangent,1),1);
            for i = 1:nphi
                I = I + double(phii>=edges(i));
            end

            dispang = zeros(nphi,1);
            for i = 1:nphi
                dispang(i) = acosd(rms(w(I==i)));
            end
            phi = edges(1:end-1);
        end
        
        function dispang = dispersionangle3d(tangent)
%DISPERSIONANGLE3D    Dispersion angle based on rms of cosine.
%   DISPERSIONANGLE3D(tangent) returns dispersion angle calculated by the
%   root-mean-square of the cosine of each unit vector tangent wrt the main
%   direction.
%
% author: Hong-Hsi Lee, 2018

            % Flip tangent to upper-half sphere
            tangent = tangent.*repmat(sign(tangent(:,3)),[1,3]);

            % Mean direction
            [U,~,~] = svd(tangent.'*tangent/size(tangent,1));

            % Align tangent to the main direction
            Lx = [0 0 0; 0 0 -1; 0 1 0]; Ly = [0 0 1; 0 0 0; -1 0 0]; Lz = [0 -1 0; 1 0 0; 0 0 0];
            uL = U(2,1)*Lx - U(1,1)*Ly; uL = uL/sqrt(U(1,1)^2+U(2,1)^2);
            thetar = acos(U(3,1));
            tangentr = expm(thetar*uL)*tangent.'; tangentr = tangentr.';

            % Calculate dispersion angle of all tangents
            w = tangentr(:,3);
            dispang = acosd(rms(w));
        end
        
        function pl = rotinv(fod,nord)
%ROTINV    Rotational invariants of the FOD.
%   ROTINV(fod,nord) returns rotaional invariants pl of fiber orientation
%   distribution fod up to the order nord. Only even orders are returned.
%   The first row of pl is p0, the second row of pl is p2, and so on. The
%   fod is the output of fodsht function, whose 'Density' must be set to
%   'low'.
%
% author: Hong-Hsi Lee, 2017
            str = which('analyzeseg.m');
            [root,~,~] = fileparts(str);
            load(fullfile(root,'TR_low.mat'))
            pts = TRCC.Points;
            wts = triarea/sum(triarea)*4*pi;
            [phi,~] = cart2pol(pts(:,1),pts(:,2));
            theta = acos(pts(:,3));
            dirs = [phi theta];
            
            fn = directSHT(nord,fod,dirs,'real',wts);
            pl = zeros(floor(nord/2)+1,1);
            plnorm = sqrt((2*(0:2:nord).'+1)/4/pi);
            for l = 0:2:nord
                list = sum(1:2:(2*l-1))+1:sum(1:2:(2*l+1));
                pl(floor(l/2)+1) = sqrt(sum(abs(fn(list)).^2));
            end
            pl = pl./plnorm;
        end
        
        function dispang = dispersionanglep2(p2)
%DISPERSIONANGLEP2    Dispersion angle based on p2.
%   DISPERSIONANGLEP2(p2) returns dispersion angle calculated from p2,
%   which is the second row of the output of rotinv function.
            dispang = acosd(sqrt(2*p2/3+1/3));
        end
        
        function [lambda,C] = rotinvpoisson(pl,l)
%ROTINVPOISSON    Fit rotational invariants to the Poisson model.
%   [lambda,C] = ROTINVPOISSON returns parameters of Poisson model for the 
%   rotational invariants pl (nl x nt) upto the order l (nl x 1).
%   Parameters include power exponent lambda and extrapolation C at l = 0.
%
% author: Hong-Hsi Lee, 2018
            A = [ones(numel(l),1),l(:)];
            X = A\log(pl);
            C = exp(X(1,:));
            lambda = exp(X(2,:));
        end
        
        function [kappa,odi] = bingham(tangent)
%BINGHAM    Fit FOD to Bingham distribution.
%   [kappa,odi] = BINGHAM(tangent) returns fitting parameter kappa, and
%   orientation distribution index odi, by fitting unit vector tangent
%   (N x 3) to the Bingham distribution.
%
% author: Hong-Hsi Lee, 2018
            tangent = [tangent;-tangent];
            B = bingham_fit(tangent);
            kappa = -B.Z;
            odi = 2/pi*atan(1./kappa);
        end
        
        function [diameter,zaxis] = smoothdiameter(diameter,zaxis,sigma,cutlength,vox)
            cutnum = ceil(cutlength/vox(3));
            zaxis = zaxis(cutnum+1:end-cutnum);
            zaxis = zaxis-zaxis(1);
            diameter = diameter(cutnum+1:end-cutnum);
            diameter = imgaussfilt(diameter,sigma/mean(diff(zaxis)));
        end
        
        function diameter = diameternuman(diameter)
            diameter = (mean(diameter.^6)/mean(diameter.^2)).^0.25;
        end
        
        function C = loglinearfit(diameter,outerdiameter)
            nl = (outerdiameter-diameter)/2;
            A = [ones(numel(diameter),1),diameter(:),log(diameter(:))];
            C = A\nl(:);
        end
        
        function gratio = loglinearmodel(diameter,C)
            gratio = diameter./(diameter + 2* (C(1) + C(2)*diameter + C(3)*log(diameter)) );
        end
        
        function [gmri,sigmagmri] = gratiomri(gratio,outerdiameter)
            gmri = sqrt(sum(gratio.^2.*outerdiameter.^2)/sum(outerdiameter.^2));
            sigmagmri = sqrt(std(gratio)^2/mean(gratio)^2 + 2*std(outerdiameter)^2/mean(outerdiameter)^2);
        end
    end
end
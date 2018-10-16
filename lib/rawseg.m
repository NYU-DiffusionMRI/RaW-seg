classdef rawseg < handle
    properties (Constant = true, Access = protected)
        
    end
    
    properties (GetAccess = public, SetAccess = protected)
    end
    
    properties (GetAccess = private, SetAccess = protected)
        
    end
    
    methods (Access = public)
        function this = rawseg(varargin)
            
        end
    end
    
    methods (Static) 
        function mask = magicwand(data,seed)
            [nx,ny,nz] = size(data);
            mask = false(nx,ny,nz);
            parfor i = 1:nz
                [xlist,ylist] = ind2sub([nx,ny],find(seed(:,:,i)));
                maski = magicwand2(squeeze(data(:,:,i)),1,xlist,ylist);
                maski = bwconvhull(~maski);
                xproj = sum(maski,2);
                yproj = sum(maski,1);
                Ix = find(xproj,1):find(xproj,1,'last');
                Iy = find(yproj,1):find(yproj,1,'last');
                maski = false(nx,ny); maski(Ix,Iy) = true;
                mask(:,:,i) = maski;
            end
        end
        
        function [datac,I,ccf,maskc] = distortioncorrect(data,mask)
            % detect distorted slices
            nk = 2;   % # neighboring slices for interp1
            th = 0.5; % threshold of corrcoef
            data = single(data);
            [nx,ny,nz] = size(data);
            I = [];
            ccf = [];
            for i = nk+1+1 : nz-nk-1
                % spline interpolation as a reference image
                datak = data(:,:,[i-nk:i-1,i+1:i+nk]);
                datak = reshape(datak,[nx*ny,2*nk]).';
                dataint = interp1([1:nk,nk+2:2*nk+1].',datak,nk+1,'spline');
                dataint = reshape(dataint.',[nx,ny]);

                % detect distorted slices with correlation coefficient
                cc = corrcoef(data(:,:,i),dataint);
                if cc(2) < th
                    I = cat(1,I,i);
                    ccf = cat(1,ccf,cc(2));
                end
            end

            % optical flow unwarping to correct distorted slices
            nk = 3;   % # neighboring slices for interp1
            [y,x] = meshgrid(1:nx,1:ny);
            datac = data;
            if nargin > 1, maskc = mask; end
            for i = 1:numel(I)
                % spline interpolation as a reference image
                Ik = I(i)-nk:I(i)+nk;
                [Ik,ia] = setdiff(Ik,I(i));
                datak = data(:,:,Ik);
                datak = reshape(datak,[nx*ny,2*nk]).';
                dataint = interp1(ia,datak,nk+1,'spline');
                dataint = reshape(dataint.',[nx,ny]);

                % Flow unwarping
                uv = estimate_flow_interface(double(dataint),double(squeeze(data(:,:,I(i)))),'classic+nl-fast');
                datac(:,:,I(i)) = interp2(squeeze(data(:,:,I(i))),x.'+squeeze(uv(:,:,1)),y.'+squeeze(uv(:,:,2)));
                if nargin ==2
                    maskc(:,:,I(i)) = interp2(squeeze(single(mask(:,:,I(i)))),x.'+squeeze(uv(:,:,1)),y.'+squeeze(uv(:,:,2))) > 0.5;
                end
            end
            datac = uint8(datac);
        end
        
        function [fiber,I] = randomhopping(medium,seed,Npar,Nstep)
            [nx,ny,nz] = size(medium);
            if isempty(Npar), Npar = 4e3; end
            if isempty(Nstep), Nstep = 16*nz^2; end
            nseed = size(seed,1);
            fiber = zeros(nx,ny,nz,'uint16');
            I = (1:nseed).';
            for i = 1:nseed
                fiberi = false(nx,ny,nz);
                xt = repmat(seed(i,:),[Npar,1]);
                for t = 1:Nstep
                    v = 6*rand(Npar,1);
                    v = [v<1, -1*(v>=1).*(v<2), (v>=2).*(v<3), -1*(v>=3).*(v<4),...
                        (v>=4).*(v<5), -1*(v>=5).*(v<6)];
                    v = [v(:,1)+v(:,2),v(:,3)+v(:,4),v(:,5)+v(:,6)];
                    xtmp = xt + v;

                    xtmp(:,1) = max(1,xtmp(:,1));
                    xtmp(:,1) = min(nx,xtmp(:,1));
                    xtmp(:,2) = max(1,xtmp(:,2));
                    xtmp(:,2) = min(ny,xtmp(:,2));
                    xtmp(:,3) = max(1,xtmp(:,3));
                    xtmp(:,3) = min(nz,xtmp(:,3));

                    list = sub2ind([nx,ny,nz],xtmp(:,1),xtmp(:,2),xtmp(:,3));
                    list1 = ~medium(list);
                    xt(list1,:) = xtmp(list1,:);
                    list2 = list(list1);
                    fiberi(list2) = true;
                end
                fiber = fiber + uint16(fiberi*i);
            end
        end
        
%         function fiberfill = fillhole(fiber,I)
%             if nargin < 2
%                 I = unique(fiber(:));
%                 if I(1) == 0, I = I(2:end); end
%             end
%             fiberfill = zeros(nx,ny,nz,'uint16');
%             se = strel('square',3);
%             for i = 1:numel(I)
%                 fiberi = fiber==I(i);
%                 for k = 1:nz
%                     fiberi(:,:,k) = imfill(imclose(fiberi(:,:,k),se),'holes');
%                 end
%                 fiberfill = fiberfill + fiberi*i;
%             end
%         end
        
%         function Iselect = select(fiber,I)
%             if nargin < 2
%                 I = unique(fiber(:));
%                 if I(1) == 0, I = I(2:end); end
%             end
%             [nx,ny,nz] = size(fiber);
%             Iselect = [];
%             for i = 1:numel(I)
%                 fiberi = fiber==I(i);
%                 fibere = imerode(fiberi,strel('cuboid',[5,5,1]));
%                 CC = bwconncomp(fibere);
%                 nVox = cellfun(@numel,CC.PixelIdxList);
%                 [~,Imax] = max(nVox);
%                 fibere(:) = false;
%                 if ~isempty(Imax)
%                     fibere(CC.PixelIdxList{Imax}) = true;
%                 end
%                 
%                 fiberd = imdilate(fibere,strel('cuboid',[7,7,3]));
%                 fiberd = fiberd.*fiberi;
%                 ratio = zeros(nz,1);
%                 for k = 1:nz
%                     fiberk = fiberi(:,:,k);
%                     fiberdk = fiberd(:,:,k);
%                     if nnz(fiberk) == 0, continue; end
%                     [Ix,Iy] = ind2sub([nx,ny],find(fiberk));
%                     if rank([Ix-Ix(1),Iy-Iy(1)])<2
%                         area_convhull = numel(Ix);
%                     else
%                         [~,area_convhull] = convhull(Ix,Iy);
%                     end
%                     ratio(k) = nnz(fiberdk)/area_convhull;
%                 end
%                 I1 = find(ratio,1);
%                 I2 = find(ratio,1,'last');
%                 ratio = ratio(I1+1:I2-1);
%                 if min(ratio) > 0.5
%                     Iselect = cat(1,Iselect,I(i));
%                 end
%             end
%         end
        
        function visualizefiber(mask,res)
            if nargin < 2
                res = [1 1 1];
            end
            [nx,ny,nz] = size(mask);
            [I,J,K] = ind2sub([nx,ny,nz],find(mask));
            I = I *res(1); J = J*res(2); K = K*res(3);
            scatter3(I,J,K,5,K,'.'); axis equal
%             xlim([1,nx*res(1)]); ylim([1,ny*res(2)]); zlim([1,nz*res(3)])
        end
        
        function animatefiber(filename,data,mask)
            rng(0); cmap = rand(128,3);
            % create 3D RGB data
            img = repmat(single(data(:)),[1,3])/single(max(data(:)));
            seg = zeros([numel(mask),3],'single');
            list = find(mask(:)>0);
            seg(list,:) = cmap(mod(mask(list),128)+1,:);
            % combine EM and segmentation using alpha compositing
            adata = 0.2; aseg = 0.8;
            imgc = img;
            imgc(list,:) = ( adata*img(list,:) + aseg*seg(list,:)*(1-adata) )/( adata + aseg*(1-adata));
            imgc = reshape(imgc,[size(mask),3]);
            imgc = permute(imgc,[1,2,4,3]);
            % write video
            v = VideoWriter(filename);
            v.FrameRate = 10;
            open(v);
            writeVideo(v,imgc);
            close(v);
        end
        
        function savefiber(filename,fiber)
            save(filename,'fiber')
        end
        
        function fiberiso = resize(fiber,vox)
            [nx,ny,nz] = size(fiber);
            vmax = max(vox);
            nxi = round(nx*vox(1)/vmax);
            nyi = round(ny*vox(2)/vmax);
            nzi = round(nz*vox(3)/vmax);
            fiberiso = imresize3(fiber,[nxi,nyi,nzi],'nearest');
        end
        
        function fiberfill = fillhole(fiber)
            fiberfill = false(size(fiber));
            zlist = find(sum(sum(fiber,1),2)).';
            for i = zlist
                fiberfill(:,:,i) = imfill(fiber(:,:,i),'holes');
            end
        end
        
        function L = watershed(fiber)
            D = bwdist(logical(fiber));
            L = watershed(D);
        end
        
        function myelin = myelinsheath(fiber,I,L,maskmy,maskfg,myelinmax,vox)
            kernelwidth = round(myelinmax/vox(1))*2 + 1; % dilation kernel width
            myelin = zeros(size(fiber),'uint16');
            parfor i = 1:numel(I)
                fiberi = fiber == I(i);
                IL = uint16(fiberi).*L; IL(IL==0) = []; IL = mode(IL);
                Li = imdilate(L==IL,strel('square',2));
                fiberdilate = imdilate(fiberi,strel('square',kernelwidth));
                fiberdilate = logical(fiberdilate) .* logical(Li) .* (~fiberi);
                myelini = fiberdilate .* maskmy .* maskfg;
                zlist = squeeze(find(sum(sum(myelini,1),2))).';
                for j = zlist
                    myelinj = squeeze(myelini(:,:,j));
                    CC = bwconncomp(myelinj);
                    numVoxels = cellfun(@numel,CC.PixelIdxList);
                    [~,Imax] = max(numVoxels);
                    myelinj = false(size(myelinj));
                    myelinj(CC.PixelIdxList{Imax}) = true;
                    myelini(:,:,j) = myelinj;
                end
                myelin = myelin + uint16(myelini)*I(i);
                i
            end
        end
        
        
        
    end
end
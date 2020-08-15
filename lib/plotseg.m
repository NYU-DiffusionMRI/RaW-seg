classdef plotseg < handle
    properties (Constant = true, Access = protected)
        
    end
    
    properties (GetAccess = public, SetAccess = protected)
    end
    
    properties (GetAccess = private, SetAccess = protected)
        
    end
    
    methods (Access = public)
        function this = plotseg(varargin)
            
        end
    end
    
    methods (Static) 
        function imgc = overlay(data,mask,alphadata,alphamask,cmap)
            ncmap = size(cmap,1);
            % create 2D RGB data
            img = repmat(single(data(:)),[1,3])/single(max(data(:)));
            seg = zeros([numel(mask),3],'single');
            list = find(mask(:)>0);
            listrem = mod(mask(list),ncmap);
            listrem(listrem==0) = ncmap;
            seg(list,:) = cmap(listrem,:);
            % combine EM and segmentation using alpha compositing
            adata = alphadata; amask = alphamask;% 0.2; 0.8
            imgc = img;
            imgc(list,:) = ( adata*img(list,:) + amask*seg(list,:)*(1-adata) )/( adata + amask*(1-adata));
            imgc = reshape(imgc,[size(mask),3]);
        end
        
        function plotbox(I,vox)
            imgu = squeeze(I(:,:,1,:));
            imgd = squeeze(I(:,:,end,:));

            imgl = squeeze(I(:,1,:,:));
            imgr = squeeze(I(:,end,:,:));

            imga = squeeze(I(1,:,:,:));
            imgp = squeeze(I(end,:,:,:));
            
            hold on
            [nx,ny,nz,~] = size(I); sx = ny/100*vox(1); sy = nx/100*vox(2); sz = nz/100*vox(3);

            xImage = [-0.5 0.5; -0.5 0.5];   % The x data for the image corners
            yImage = [-0.5 -0.5; 0.5 0.5];   % The y data for the image corners
            zImage = [0.5 0.5; 0.5 0.5];     % The z data for the image corners

            surf(xImage*sx,yImage*sy,zImage*sz,...
                 'CData',imgu,...
                 'FaceColor','texturemap','edgecolor','none');
            surf(xImage*sx,yImage*sy,-zImage*sz,...
                 'CData',imgd,...
                 'FaceColor','texturemap','edgecolor','none');
            surf(-zImage*sx,xImage*sy,yImage*sz,...
                 'CData',rot90(imgl),...
                 'FaceColor','texturemap','edgecolor','none');
            surf(zImage*sx,xImage*sy,yImage*sz,...
                 'CData',rot90(imgr),...
                 'FaceColor','texturemap','edgecolor','none');
            surf(xImage*sx,-zImage*sy,-yImage*sz,...
                 'CData',fliplr(rot90(imga,-1)),...
                 'FaceColor','texturemap','edgecolor','none');
            surf(xImage*sx,zImage*sy,-yImage*sz,...
                 'CData',fliplr(rot90(imgp,-1)),...
                 'FaceColor','texturemap','edgecolor','none');
            axis equal off
            view([37.5 35]); grid on
            light('Position',[-1,-1,1],'Style','infinite')
            light('Position',[-1,-1,0.7],'Style','infinite')
            light('Position',[0.5,-1,0],'Style','infinite')
            material dull
            camva;
        end
        
        function TR = bw2triangulation(BW,varargin)
            sf = 0.95;
            if nargin >1
                sf = varargin{1};
            end
            [I,J,K] = ind2sub(size(BW),find(BW));
            [k,~] = boundary(I,J,K,sf);
            TR = triangulation(k,I,J,K);
        end
    end
end
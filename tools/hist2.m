function h = hist2(x, y, Nbins, xlimits, ylimits, col)

if nargin > 3, 
    hh = histogram2(x, y, 'NumBins', Nbins, 'xbinlimits', xlimits, 'ybinlimits', ylimits, 'visible', 'off', 'normalization', 'probability'); 
else    
    hh = histogram2(x, y, 'NumBins', Nbins, 'visible', 'off', 'normalization', 'probability'); 
end

if isempty(col), col = 'blue'; end 

% [Y,X] = meshgrid(hh.YBinEdges(1:end-1)+0.5*hh.BinWidth(1), hh.XBinEdges(1:end-1)+0.5*hh.BinWidth(2)); Z = hh.Values; 
[Y,X] = meshgrid(hh.YBinEdges(1:end-1)*0.5 + hh.YBinEdges(2:end)*0.5, hh.XBinEdges(1:end-1)*0.5 + hh.XBinEdges(2:end)*0.5); Z = hh.Values;

surf(X, Y, Z, 'AlphaData', Z, 'FaceAlpha', 'Interp', 'FaceColor', col, 'EdgeColor', 'none'); view(0,90); axis equal square
    
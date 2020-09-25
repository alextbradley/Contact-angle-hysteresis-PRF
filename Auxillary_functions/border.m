function [X, Y, xgrid, ymin, ymax] = border(xx, yy, npts)
%Returns a vector of border points surrounding points xx, yy for use in
%fill.
    xgrid = linspace(min(xx),max(xx), npts);
    dx = diff(xgrid);
    dx = dx(1);
    ymin = nan(1,length(xgrid));
    ymax = nan(1,length(xgrid));
    
    for i = 1:length(xgrid)
        idx = (abs(xgrid(i) - xx) < dx/2); %points within dx of xgrid(i)
        ymin(i) = min(yy(idx));
        ymax(i) = max(yy(idx));
    end
    
    X = [xgrid, flip(xgrid)];
    Y = [ymin, flip(ymax)];
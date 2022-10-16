function figpref(varargin)
% This selects default font sizes and line widths to make figures
% more suitable for specific destinations (for example, larger
% fonts for small bitmaps files intended for presentation on the
% Web)
%
% Usage:
%   o  place in your matlab path
%   o  run without parameters for a menu
%   o  you can skip the menu by supplying the desired choice as the
%      first parameter
%
% Version: 0.0.1
%
% Author: Colin Macdonald (http://www.math.sfu.ca/~cbm)
% license: GPL
%
% Please send me improvements, etc at cbm[at]sfu[dot]ca

% Option 4 added by Dominic Vella
  
  if (nargin == 0) 
    disp(' Setup figure preferences');
    disp('   0) on screen viewing');
    disp('   1) small bitmap creation (for Web)');
    disp('   2) small .eps file creation (for side-by-side on page)');
    disp('   3) large .eps file creation (full width figures)');
    disp('   4) Full page width in LaTeX');
    disp('   9) other');
   
    wh = input('choose #:  ');
  else
    wh = varargin{1};
  end
  
  % First, reset everything to the way it was when matlab started
  % I found these defaults by starting matlab and using
  % get(0, 'whatever')
  % IMPORTANT: This should reset everything changed by the other
  % modes!
  set(0, 'DefaultAxesFontSize', 10)
  set(0, 'DefaultAxesLineWidth', .5);
  set(0, 'DefaultLineLineWidth', .5);
  set(0, 'DefaultPatchLineWidth', .5);
  set(0, 'DefaultLineMarkerSize', 6);
  
  switch wh
   case 0
    disp('setting figure preferences for on screen viewing');
    % this is done by the defaults above
    
   case 1
    disp('setting figure preferences for small bitmap creation');
    disp('** This is intended to create small (-r40) bitmaps');
    disp('   for example for Web presentation.');
    set(0, 'DefaultAxesFontSize', 20)
    set(0, 'DefaultLineMarkerSize', 12);
    
   case 2
    disp('setting figure preferences for small .eps file creation');
    disp('** This is intended to create figures which will be');
    disp('   displayed side-by-side with another figure in LaTeX.');
    % this is what Trefethen used to create the figures for
    % "Spectral Methods in Matlab"
    set(0, 'defaultaxesfontsize', 12);
    set(0, 'defaultaxeslinewidth', .7);
    set(0, 'defaultlinelinewidth', .8);
    set(0, 'defaultpatchlinewidth', .7);

   case 3
    disp('setting figure preferences for standard figures file creation');
    set(0, 'defaultaxesfontsize', 20);
    set(0, 'defaultaxeslinewidth', 1.5);
    set(0, 'defaultlinelinewidth', 1.4);
    set(0, 'defaultpatchlinewidth', .7);
        
   case 4
    disp('** This is intended to create figures which will use');
    disp('   the full width of a page in LaTeX.');
    set(0, 'defaultaxesfontsize', 14);
    set(0, 'defaultaxeslinewidth', 2);
    set(0, 'defaultlinelinewidth', 1.5);
    set(0, 'defaultpatchlinewidth', 0.7); 
 
   case 9
    disp('Please modify the source code to add new modes');
    disp('or modify existing ones.');
    disp('Please contribute changes back to cbm[at]sfu[dot]ca');
    
   otherwise
    warning('invalid figures preference call');
  end
function [fitresult,gof] = LinearFit1(x,y)
%LINEARFIT1 Secondary fitting program used for testing and forcing bad data
%fits
%   A more sturctured linear fit that is better at forceing a startpoint to
%   ensure good drift calcs on bad data. 
%% Fit:
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'm*x+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.406924655543266 0.278498218867048];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end


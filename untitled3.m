load("Geigerdata.mat");
nuclearCounts = VarName1;
numPts = numel(nuclearCounts);

nBins = round(sqrt(numPts)); % define the number of bins for histogram
h = histogram(nuclearCounts,nBins)% Create a handle (h) for the histogram
xData = h.BinEdges(1:end-1)'+ h.BinWidth/2; % convert a row vector into a
%column vector & assign xData to midpoints of bins
yData = h.Values';% 

% Assign initial parameter values - best guesses.

amplitudeStart = 100;
sigmaStart = 15;
meanCountsStart = 100;
11
%Define options for fit - parameter ranges, method.

fitOpts = fitoptions('Method','NonlinearLeastSquares',...
 'Lower',[-10*amplitudeStart -10*meanCountsStart -10*sigmaStart],...
 'Upper',[10*amplitudeStart 10*meanCountsStart 10*sigmaStart],...
 'StartPoint', [amplitudeStart meanCountsStart sigmaStart]);

% Define model function

ft = fittype( 'ampl*exp(-((x-meanX).^2/(2*sigma*sigma)))',...
'coefficients',{'ampl','meanX','sigma'}, 'independent', {'x'},...
'dependent', 'y');

% Peform fit using 'fit' function in Matlab

[fitresult, gof, fout] = fit(xData, yData, ft, fitOpts);

%Print results of fitting the model to the data
amplitude = fitresult.ampl
meanCounts = fitresult.meanX
sigma = fitresult.sigma
gof
% Compute results of fitting the model to the data
yFit = amplitude*exp(-((xData-meanCounts).^2/(2*sigma*sigma)));

dof = nBins-3; %number of degrees of freedom (nPoints-no. of fit
%parameters)
reducedChiSquare = (1/dof)*sum(((yFit-yData).^2)./(yData+2))

hfigure = figure('Color', 'w');
%plot bar graph and fitted curve
hAx(1) = axes();
h_bar = bar(xData, yData, 'b', 'Parent',hAx(1));% bar graph experimental
%data with blue = 'b' bars
set(hAx(1), 'Box','off')
hlegend1 = legend(h_bar, 'Exp Data', 'Location','NorthEast'); %display
%legend 1
hold on
hAx(2) = copyobj(hAx(1),gcf); %# copy first axes to second plot
delete( get(hAx(2),'Children') ) %# delete its children
hfit = plot(fitresult);
set(hfit,'LineWidth',3, 'Color', 'r');
title('1D Gaussian fit to Decay Counts Frequency Distribution');
% Label axes
xlabel( 'counts' );
ylabel( 'frequency' );
set(gcf, 'Color', 'w');
set(hAx(2), 'Color','none', 'XTick',[], ...
 'YAxisLocation','right', 'Box','off') %# make it transparent
grid on


%om = sigma/sqrt(nuclearCounts)
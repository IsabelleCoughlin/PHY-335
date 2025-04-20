
N = %  total number of samples
p = % probability of decay
n = p*N %the average number of counts in time interval
o = sqrt(n) %standard deviation
om = o/sqrt(N) %standard error

histogram(%tablename)


filename = filename = '/Users/isabe/Downloads/Lab1.xlsx';
nuclearDecayData =importdata(filename);
nuclearCounts = nuclearDecayData.data(:,2);
numPts = numel(nuclearCounts);% 
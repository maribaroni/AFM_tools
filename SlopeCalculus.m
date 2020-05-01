%%
% Script responsible for getting the trace and retrace curves from a
% PFC file (peak force - Nanoscope file) and calculating the slope.
%
% NSMatlabUtilities toolbox (from Brucker corporation) is needed to get the curves.
%
% @author: Mariana P. M. A Baroni, PhD
% @last access: April 30, 2020
%
% Don't forget to cite it properly!


%****Sanitizing
clear all %clear variables
close all %clear figures
clc %clear command window
%*****

NSMU = NSMatlabUtilities();%invoking the toolbox

%******* Current folder needed to get the PFC file.
currentFolder = pwd;
NSMU.Open(strcat(currentFolder, '\TesteLuciana.pfc')); %change the file name here (Windows based system)

%******* Getting the number of curves from the PFC file using the toolbox
[NumberOfChannels] = GetNumberOfChannels(NSMU);
NumberOfCurves = NSMU.GetNumberOfForceCurves();
NumberOfCurves = NumberOfCurves - NumberOfCurves*.1; %dispose the last 10% of curves

%******* Saving memory to store all curves on different files (Trace and retrace)
%If you know the dimension number
n = 128;%dimension number (taked from AFM file)
matrixTrace = zeros(n,n);
matrixRetrace = zeros(n,n);

%******* Saving memory to store the slope of all curves
%if you know how many slopes that you will have
%vector = zeros(n,1);

%******* Loop to get each curve
interval = 1; % interval is a parameter to set an interval to get the curves
for i = 1 : interval : NumberOfCurves
    
    [trace, retrace, scaleUnit, dataTypeDesc] = GetForceVolumeForceCurveData(NSMU, NumberOfChannels, i, NSMU.METRIC);
    %******* Plot curves
    % figure(1)
    % plot(trace)
    % hold on
    % title('Using GetForceVolumeCurveData - TRACE')
    % figure(2)
    % plot(retrace)
    % hold on
    % title('Using GetForceVolumeCurveData - RETRACE')
    
    %Set the curve
    %f = flipud(trace); %flip array in up/down direction %TRACE
    f = flipud(retrace); %flip array in up/down direction %RETRACE
    
    %******* Calculate the slope
    %All points above x-axis
    minimum = min(f);
    f = f + abs(minimum);
    
    %Interval to calculate the slope
    % n0 = ceil(size(f,1)*.10); %dispose the 10% first data (without flip)
    % nf = ceil(size(f,1)*.80); %working with the 80% first data (without flip)
    n0 = ceil(size(f,1)*.85); %dispose the 85% first data (with flip)
    nf = ceil(size(f,1) - size(f,1)*.01); %dispose the 1% last data (with flip)
    
    h = 1; %Assumption: distance 1 between the data points
    
    interval = n0+1 : h : nf-1;
    for j = interval
        
        %calculating the derivatives
        K1 = firstDerivative(f(j-1), f(j), f(j+1) ,h);
        K2 = secondDerivative(f(j-1), f(j), f(j+1) ,h);
        
        error = 0.01; %plausible error
        
        %we also can use the minimum value of function
        %    if (f(j) <= error)
        %        index = j;
        %        break;
        %    end
        
        % Using derivatives
        
        % if (K2<0) %%in case of second derivative (assumption negative)
        %  slope) or (K2<0) assumption positive slope
        if( (abs(K1) <= error) || (K2 < 0) ) %using both first or second derivative
            index = j;
            break;
        end
        
    end
    
    %z0 = index - 1; %store z0 considering considering the value before the index founded (without flip)
    z0 = index + 2; %flip
    
    %Linear polyfit
    [slope] = polyfit(n0:z0, f(n0:z0)',1);
    vector(i) = slope(1); %store the angular coefficient
    
    %%Plot
    % figure(3)
    % plot(f); %curva
    % hold on
    % y = polyval(slope,n0:z0);
    % plot(n0:z0,y, '-*r');
    % pause
    
    %*******Storing the curves in differente matrices
    matrixTrace(:,i) = trace;
    matrixRetrace(:,i) = retrace;
    
    i; %print the number of curve in working (just to follow the execution)
    
end

%Histogram of slope
figure(4)
histogram(vector, 50) %50 bins

figure(5)
plot(vector)

%Saving the matrices on files (take several minutes to finish)
%save('TesteTrace.txt','matrixTrace', '-ascii');
%save('TesteRetrace.txt','matrixRetrace', '-ascii');

%%%%%%%%%%%%%%%%%%%%%%FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******Function to calculate the first derivative
%Assumption: h constant
function [K1] = firstDerivative(fx0, fx1, fx2, h)

%centered differences - ordem error 4
K1 = (fx2 - fx0) / (2 * h);

end

%*******Function to calculate the second derivative
%Assumption: h constant
function [K2] = secondDerivative(fx0, fx1, fx2, h)

%centered differences - ordem error 4
K2 = ( fx2 - (2 * fx1) + fx0 ) / (h * h);

end

%******Closing toolbox
%NSMU.close();

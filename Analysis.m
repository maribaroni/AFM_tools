%%
% Script responsible for getting the trace and retrace curves from a
% PFC file (peak force - Nanoscope file) and calculating the slope
% and Young's Modulus.
%
% NSMatlabUtilities toolbox (from Brucker corporation) is needed
% to get the curves.
%
% @author: Mariana P. M. A Baroni, PhD
% @last access: February 01, 2021
%
% Don't forget to cite it properly!

%****Sanitizing
clear all %clear variables
close all %clear figures
clc %clear command window
%*****

NSMU = NSMatlabUtilities();%invoking the toolbox

% GUI element to allow selection of a directory where the data is
fileLocationsName = uigetdir();

% Get only .pfc files from this directory
fileLocations = dir(fullfile(fileLocationsName, '*.pfc'));

%Store all curves
matrixRetrace = [];

%Store all Modulus Young
matrixYoung = [];

%Store all Slope
matrixSlope = [];

nameFiles = string.empty;

count2 = 1;
% Loop to go through each file
for currentFile = 1:length(fileLocations)
    
    %Retrieve the file name
    fileName=fileLocations(currentFile).name;
    
    disp('===============================================================================================================');
    disp(' Analysing New file!');
    disp([' ' fileName]);
    disp('===============================================================================================================');
    
    directoryFileName=char(strcat(fileLocationsName, '\', fileName));
    
    %Open up the file
    NSMU.Open(directoryFileName);
    
    %Return image pixels and number of force curves
    %in each scan line of the specific force volume/peak force
    %capture file.
    [sizeOfImage,forVolPixel] = NSMU.GetForceVolumeScanLinePixels();
    
    %Returns the number of force curves in the file
    NumberOfCurves = NSMU.GetNumberOfForceCurves();
    
    %Returns the Poisson ratio from header file
    PR = NSMU.GetPoissonRatio();
    
    %Returns the Tip radius from the header file
    TR = NSMU.GetTipRadius();
    
    %Store a vector with each value of Youg's Modulus
    modulusYoung = double.empty;
    
    %Store a vector with each value of Youg's Modulus
    slopesPFC = double.empty;
    
    %Store error curves
    errorCurves = double.empty;
    
    count = 1;
    count1 = 0;
    
    %Loop to calcule the slope and Young's Modulus for each curve
    for i = 1 : NumberOfCurves
        
        count1 = count1 + 1;

        [chanDesc] = NSMU.GetDataTypeDesc(2); %Chan 1 is Image channel, so grab force curve data type from chan 2

        %get F vs tip-sample separation plot of curve i
        [xTrace, xRetrace, yTrace, yRetrace, xLabel, yLabel] = NSMU.CreatePeakForceForceCurveZplot(i, NSMU.METRIC, 1); %isSeparation? yes = 1, no = 0
        minimum = yRetrace(size(yRetrace,1));
        yRetrace = yRetrace + abs(minimum);
                
        %Store retrace
        matrixRetrace(:,count1) = yRetrace;%flipud(yRetrace);
        
        %Get the index of contact point
        contactPointIndex = FindContactPoint(xRetrace, yRetrace);
        
        %Compute the adhesion region (adhesion force)
        [regionBottom, regionTop] = ComputeContactRegionBounds(xRetrace, yRetrace, contactPointIndex, 10, 90);
        
        %Indexes of adhesion region
        [index1, index2] = ComputeMarkers(xRetrace, yRetrace, regionTop, regionBottom);
        
        %Curves with problems
        if ( index1 <= 0 || index2 <= 0)
            errorCurves(count) = i;
            count = count + 1;
        else
            
            %Calculate the force/deflection
            K = ExponentialFit(xRetrace, yRetrace, index1, index2, contactPointIndex);
            
            %Calculate the Young's modulus
            modulus = GetYoungsModulus(K, PR, TR);
            
            %Store the Young's modulus
            modulusYoung(count1)=modulus;
            
            %Linear polyfit
            %   n0 = index2;
            %  nf = index1;
            n0 = index1;
            nf = index2;
            
            [slope] = polyfit(xRetrace(n0:nf), yRetrace(n0:nf)',1);
            slopesPFC(count1) = slope(1); %-1*slope(1); %store the angular coefficient
            
            %Evaluate the Polinomio
            y = polyval(slope,xRetrace(n0:nf));
            
            %Plot each %5000 curves
            %         %Plot all curves
            if mod(i,5000) == 0
                fig = figure(1);
                plot (xRetrace, yRetrace, 'b-','LineWidth',1);
                xlabel(xLabel);
                ylabel(yLabel);
                hold on;
                plot(xRetrace(n0:nf),y, '-*r');
                %  set(gca,'xdir','reverse','FontSize',14);
            end
            
        end
        
    end
    
    fig1 = figure(2);
    histogram(modulusYoung,100);
    xlabel(" Young's Modulus (MPa) ");
    set(gca,'FontSize',14);
    
    fig2 = figure(3);
    histogram(slopesPFC,100);
    xlabel(" Slope (nm/V)");
    set(gca,'FontSize',14);
    
    %Save files
    newFileName = extractBefore(fileName,".pfc");
    str = strcat(fileLocationsName, "\", newFileName,"_Young.txt");
    str1 = strcat(fileLocationsName, "\",newFileName,"_Slopes.txt");
    str2 = strcat(fileLocationsName, "\",newFileName,"_ErrorCurves.txt");
    str3 = strcat(fileLocationsName, "\",newFileName,"_YoungHistogram.jpg");
    str4 = strcat(fileLocationsName, "\",newFileName,"_SlopeHistogram.jpg");
    str5 = strcat(fileLocationsName, "\",newFileName,"_PFCCurves.jpg");
    
    save(str, 'modulusYoung', '-ascii');
    save(str1, 'slopesPFC', '-ascii');
    save(str2, 'errorCurves', '-ascii');
    
    saveas(fig, str5);
    saveas(fig1, str3)
    saveas(fig2, str4)
    
    matrixYoung(:,currentFile) = modulusYoung;
    matrixSlope(:,currentFile) = slopesPFC;
    nameFiles(currentFile) = newFileName;
    
    disp('***Finish***')
    
    close all
    
end

fig3 = figure(4);
boxplot(matrixYoung,nameFiles) %whisker ommitts outliers
xlabel('Samples')
ylabel("Young's Modulus (MPa)")
set(gca,'FontSize',14);

fig4 = figure(5);
boxplot(matrixSlope,nameFiles, 'whisker',10^8);
xlabel('Samples')
ylabel("Slope (nm/V)")
set(gca,'FontSize',14);

%Save files for all curves
str6 = strcat(fileLocationsName, "\","AllPFCCurves.txt");
save(str6,'matrixRetrace', '-ascii');

str7 = strcat(fileLocationsName, "\","AllModulusYoung.txt");
save(str7,'matrixYoung', '-ascii');

str8 = strcat(fileLocationsName, "\","AllSlopes.txt");
save(str8,'matrixSlope', '-ascii');

str9 = strcat(fileLocationsName, "\","BoxPlotYoung.jpg");
saveas(fig3, str9)

str10 = strcat(fileLocationsName, "\","BoxPlotSlope.jpg");
saveas(fig4, str10)

disp('===============================================================================================================');
disp('Starting PCA analysis...')

%+++++++++++++++++++++++++++++++PCA Analysis of PFC curves: pca(X)
%X is a Matrix, Rows of X correspond to observations and columns correspond to variables.

%------------- PCA - Matrix Retrace
%Each column of matrixRetrace is a observation, so, it is necessary
%transpose it
[coeff,score,latent,tsquared,explained,mu] = pca(matrixRetrace);
fig5 = figure(6);
gscatter(score(:,1),score(:,2))
xlabel(strcat('1st Principal Component ',num2str(explained(1))))
ylabel(strcat('2nd Principal Component ',num2str(explained(2))))
set(gca,'FontSize',14);
str11 = strcat(fileLocationsName, "\","PCA2_RetraceCurves.jpg");
saveas(fig5, str11)
% %Biplot
% clear figure(6) fig5
% fig5 = figure(6);
% vlabs = nameFiles;
% biplot(coeff(:,1:2),'Scores',score(:,1:2),'Color','b','Marker','o','VarLabels',vlabs);
% set(gca,'FontSize',14);
% str11 = strcat(fileLocationsName, "\","BiPlot2_RetraceCurves.jpg");
% saveas(fig5, str11)
% clear figure(6) fig5
% fig5 = figure(6);
% vlabs = nameFiles;
% biplot(coeff(:,1:3),'Scores',score(:,1:3),'Color','b','Marker','o','VarLabels',vlabs);
% set(gca,'FontSize',14);
% str11 = strcat(fileLocationsName, "\","BiPlot3_RetraceCurves.jpg");
% saveas(fig5, str11)


fig6 = figure(7);
scatter3(score(:,1),score(:,2),score(:,3))
xlabel(strcat('1st Principal Component ',num2str(explained(1))))
ylabel(strcat('2nd Principal Component ',num2str(explained(2))))
zlabel(strcat('3rd Principal Component ',num2str(explained(3))))
set(gca,'FontSize',14);
str12 = strcat(fileLocationsName, "\","PCA3_RetraceCurves.jpg");
saveas(fig5, str12)
clear coeff score

%------------- PCA - Matrix Slope
%Each column of matrixSlope is a observation, so, it is necessary
%transpose it
[coeff,score,latent,tsquared,explained,mu] = pca(matrixSlope');
fig7 = figure(8);
gscatter(score(:,1),score(:,2),nameFiles')
xlabel(strcat('1st Principal Component ',num2str(explained(1))))
ylabel(strcat('2nd Principal Component ',num2str(explained(2))))
set(gca,'FontSize',14);
str11 = strcat(fileLocationsName, "\","PCA2_Slope.jpg");
saveas(fig7, str11)
% %Biplot
% clear figure(8) fig7
% fig7 = figure(8);
% vlabs = num2cell(nameFiles);
% biplot(coeff(:,1:2),'Scores',score(:,1:2),'Color','b','Marker','o','VarLabels',vlabs);
% set(gca,'FontSize',14);
% str11 = strcat(fileLocationsName, "\","BiPlot2_Slope.jpg");
% saveas(fig7, str11)
% clear figure(8) fig7
% fig7 = figure(8);
% vlabs = nameFiles;
% biplot(coeff(:,1:3),'Scores',score(:,1:3),'Color','b','Marker','o','VarLabels',vlabs);
% set(gca,'FontSize',14);
% str11 = strcat(fileLocationsName, "\","BiPlot3_Slope.jpg");
% saveas(fig7, str11)

fig8 = figure(9);
scatter3(score(:,1),score(:,2),score(:,3));%,nameFiles')
xlabel(strcat('1st Principal Component ',num2str(explained(1))))
ylabel(strcat('2nd Principal Component ',num2str(explained(2))))
zlabel(strcat('3rd Principal Component ',num2str(explained(3))))
set(gca,'FontSize',14);
str12 = strcat(fileLocationsName, "\","PCA3_Slope.jpg");
saveas(fig8, str12)

clear coeff score

%------------- PCA - Matrix Young
%Each column of matrixYoung is a observation, so, it is necessary
%transpose it
[coeff,score,latent,tsquared,explained,mu] = pca(matrixYoung');
fig9 = figure(10);
gscatter(score(:,1),score(:,2),nameFiles')
xlabel(strcat('1st Principal Component ',num2str(explained(1))))
ylabel(strcat('2nd Principal Component ',num2str(explained(2))))
set(gca,'FontSize',14);
str11 = strcat(fileLocationsName, "\","PCA2_Young.jpg");
saveas(fig9, str11)
% %Biplot
% clear figure(10) fig9
% fig9 = figure(10);
% vlabs = nameFiles;
% biplot(coeff(:,1:2),'Scores',score(:,1:2),'Color','b','Marker','o','VarLabels',vlabs);
% set(gca,'FontSize',14);
% str11 = strcat(fileLocationsName, "\","BiPlot2_Young.jpg");
% saveas(fig9, str11)
% clear figure(10) fig9
% fig9 = figure(10);
% vlabs = nameFiles;
% biplot(coeff(:,1:3),'Scores',score(:,1:3),'Color','b','Marker','o','VarLabels',vlabs);
% set(gca,'FontSize',14);
% str11 = strcat(fileLocationsName, "\","BiPlot3_Young.jpg");
% saveas(fig9, str11)

fig10 = figure(11);
scatter3(score(:,1),score(:,2),score(:,3));%,nameFiles')
xlabel(strcat('1st Principal Component ',num2str(explained(1))))
ylabel(strcat('2nd Principal Component ',num2str(explained(2))))
zlabel(strcat('3rd Principal Component ',num2str(explained(3))))
set(gca,'FontSize',14);
str12 = strcat(fileLocationsName, "\","PCA3_Young.jpg");
saveas(fig10, str12)

clear coeff score


function [contactPointIndex] = FindContactPoint(xData, yData)
% Finds intersection point of the curve and line connected the first and last point.
% curveStartIndex, contactPointIndex andcurveEndIndex are 0 if no intersection exists.

curveStartIndex = 0;
contactPointIndex = 0;
curveEndIndex = 0;
xSize = max(size(xData));
ySize = max(size(yData));
%if(xSize == ySize)

%Method for contact point is subtract the line connected the first and last point from all points in the curve, and use the lowest point
slope = (yData(ySize) - yData(1)) / (xData(xSize) - xData(1));
minY = 0;
for i = 1:ySize
    yVal = yData(i) - slope * xData(i);
    if (i == 1 || yVal < minY)
        minY = yVal;
        contactPointIndex = i;
    end
end
%Determine direction of curve. Method assumes end of the contact region is higher that the non-contact region
if (yData(1)> yData(ySize))
    curveStartIndex = ySize;
    curveEndIndex = 1;
else
    curveStartIndex = 1;
    curveEndIndex = ySize;
end;
%Method is find highest point after contact point
if curveEndIndex > contactPointIndex
    increment = 1;
else
    increment = -1;
end
maxY = 0;
i = contactPointIndex + increment;
while i * increment <= curveEndIndex * increment
    if yData(i) > maxY
        maxY = yData(i);
    end
    i = i + increment;
end
targetMaxY = (maxY - yData(contactPointIndex)) + yData(contactPointIndex);
maxIndex = contactPointIndex;
i= contactPointIndex + increment;
while i * increment <= curveEndIndex * increment
    if (yData(i) > targetMaxY)
        break;
    end
    maxIndex = i;
    i = i + increment;
end
curveEndIndex = maxIndex;
end

function [regionBottom, regionTop] = ComputeContactRegionBounds(xData, yData, contactPointIndex, minForcePercent, maxForcePercent)

[curveStartIndex, curveEndIndex] = deal(0);
xSize = max(size(xData));
ySize = max(size(yData));

if xData(1) < xData(xSize)
    curveStartIndex = xSize;
    curveEndIndex = 1;
else
    curveStartIndex = 1;
    curveEndIndex = xSize;
end
maxRegion = yData(curveEndIndex);
minRegion = yData(contactPointIndex);
regionBottom = (maxRegion - minRegion) * minForcePercent / 100 + minRegion;
regionTop = (maxRegion - minRegion) * maxForcePercent / 100 + minRegion;
end


function [index1, index2] = ComputeMarkers(xData, yData, regionTop, regionBottom)
[index1, index2] = deal(-1);
nPts = max(size(yData));
step = 1;
startIndex = 1;
endIndex = nPts;
%if it's reverse index
%aqui    if xData(1) > xData(nPts)
%aqui        startIndex = nPts;
%aqui        endIndex = 1;
%aqui       step = -1;
%aqui    end
i = startIndex;
while i~=endIndex
    if (index1 == -1)
        if (yData(i) <= regionTop)
            index1 = i;
        end
    elseif (yData(i) <= regionBottom)
        index2 = i - step;
        break;
    end
    i = i+ step;
end
end


function E = GetYoungsModulus(K, PoissonRatio, TipRadius)
%Returns YoungsModulus, E in MPa
%Use Hertz sphere model to calculate Young's Modulus.
PoissonRationSqrd = PoissonRatio^2;
sqrtTipRadius = sqrt(TipRadius);
%Hertz model: K = 4/3 * E/(1-PR^2) * R ^ 1/2
E = (1 - PoissonRationSqrd) * .75 * K/sqrtTipRadius;
end

function [K,scale] = ExponentialFit (xData, yData, index1, index2, contactPointIndex)
K = 0;
xSize = max(size(xData));
ySize = max(size(yData));
if xSize ~= ySize
    error('xData must be the same size as yData.')
end
nBins = abs(index2 - index1) + 1;
startIndex = min(index1, index2);
endIndex = max(index1, index2);
%add the contact regions to the vectors
fitXData = xData(startIndex:(endIndex - 1));
fitYData = yData(startIndex:(endIndex - 1));
%add the contact point to the front of the vectors
fitXData = [xData(contactPointIndex); fitXData];
fitYData = [yData(contactPointIndex); fitYData];
%transfer the data so x=0 is the first data point
x0 = fitXData(1);
x1 = fitXData(nBins);
scale = abs(x1 - x0)/(x1 - x0);
yMin = fitYData(1);
%not flipped
if fitYData(1) <= fitYData(nBins)
    fitXData = (fitXData - x0) * scale;
    %flipped
else
    yMin = fitYData(nBins);
    fitXData = (x1 - fitXData)*scale;
end
%get y in PN
fitYData = (fitYData - yMin) * 1000;
%get the first guess, K = mean(F(x)/x^E)
n = 0;
for i=1:nBins
    if fitXData(i)>0 && fitYData(i)>=0
        K = K + fitYData(i)/fitXData(i)^1.5;
        n = n+1;
    end
end
if n > 0
    K = K /n;
end
end

%%
% Script responsible for ..... force volume
%
% NSMatlabUtilities toolbox (from Brucker corporation) is needed to get the curves.
%
% @author: Mariana P. M. A Baroni
% @last access: July 15, 2020
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
name='30hz'; %can be 1hz, 5hz, 15hz 
NSMU.Open(strcat(currentFolder, '\novos testes', '\NE_L_', name)); %change the file name here (Windows based system)

%[deflSens, deflSensUnits] = GetDeflectionSensitivity(NSMU, true);

%get number of image pixels & number of force curves in each scan line
[imagePixel, forVolPixel] = NSMU.GetForceVolumeScanLinePixels();
NumberOfCurves = NSMU.GetNumberOfForceCurves();

%---------------------------------
%Display image
f = figure();
movegui(f, 'northwest');

%Get FV Image Data
[data, scaleUnit, dataTypeDesc] = NSMU.GetForceVolumeImageData(NSMU.METRIC);

image(flipud(data), 'CDataMapping', 'scaled');
set(gca, 'YDir', 'normal');
axis('tight', 'square');
colormap('gray')%('Copper');
hold on;
plot(1, 1, 's', 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b')
plot(imagePixel, 1, 's', 'MarkerSize', 10, 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c')
plot(imagePixel, imagePixel, 's', 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g')
%colorbar();

title(strcat('Force Volume', {' '}, dataTypeDesc, ' Image'));
xLabel = NSMU.GetScanSizeLabel();
xlabel(xLabel);

newdata = flipud(data);
save(strcat('FV_',name,'.txt'),'newdata', '-ascii');

NSMU.Close();

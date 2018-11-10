%% SCRABBLE DEMO: 
% A representative imputation result using simulated data based on 
% down-sampling of real bulk RNA-Seq data. The drop-out rate parameter \labmda 
% equals to 0.1 that corresponds to 82% zeros in the data matrix.

%% Clear all variables
clear all
clc
close all

%% Load the data
% There are three datasets in the .mat file. There are the true data set,
% Drop-out data set, and the imputed data set by SCRABBLE.
load('demo_data_HF.mat')

%% Prepare the data
% We construct the data structure which is taken as one of the input of
% *scrabble* https://geo2.ggpht.com/cbk?panoid=Tj72QX_B4DCrqVcMhsCUsA&output=thumbnail&cb_client=search.TACTILE.gps&thumb=2&w=408&h=200&yaw=238.76057&pitch=0&thumbfov=100function.
data.data_sc = data_sc;
data.data_bulk = data_bulk;

%% Prepare the parameter for SCRABBLE
% set up the parameters used in example
parameter = [100,2e-7];
nIter = 100;

%% Run SCRABBLE
dataRecovered = scrabble(data,parameter,nIter);

%% Plot the results
gcf = figure(1);
set(gcf, 'Position', [100, 500, 1200, 300])
subplot(1,3,1)
imagesc(log10(data_true+1))
title('True Data')
axis off
subplot(1,3,2)
imagesc(log10(data_sc+1))
title('Drop-out Data')
axis off
subplot(1,3,3)
imagesc(log10(dataRecovered+1))
title('Imputed Data by SCRABBLE')
axis off
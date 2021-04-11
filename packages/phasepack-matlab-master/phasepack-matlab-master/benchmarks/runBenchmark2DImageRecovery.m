%                           runBenchmark2DImageRecovery.m
% 
% This example benchmarks algorithms based on their ability to reconstruct
% an image using synthetic measurements (random Gaussian).  The benchmark 
% shows how the different methods behave as the number of Fourier masks
% increases.
% 
% 
% This script does the following: 
% 1. Set up parameters and create a list of
% algorithm structs. 
% 
% 2. Invoke the general benchmark function benchmarkPR.  This function will
% reconstruct the image using different numbers of masks, and produce a
% curve for each algorithm.
%
%
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017

%% -----------------------------START----------------------------------

function [rec,algorithms_return] = runBenchmark2DImageRecovery(opt)

%% 1.Set up parameters
% Choose x label, values shown on the x axis.  For 2D images, we use masks
% to acquire Fourier measurements.  More masks means more measurements.
xitem = 'masks';            % Vary the number of masks on the x-axis
xvalues = [opt.numMasks];  % Choose these values for the number of masks
yitem = 'reconerror';       % Plot the reconstruction error on the y-axis


% Tell the benchmark function to use a 2D image
dataSet = '2DImage';

% Set up general parameters
params.verbose = false;
params.imagePath = 'data/shapes.png';   % Where is the image located?
params.numTrials = 1;                   % How many random trials to run for each algorithm

% Create a list of algorithms structs
wf = struct('initMethod','spectral','algorithm','wirtflow');
twf = struct('algorithm','twf'); 
rwf = struct('algorithm','rwf');
ampflow = struct('algorithm','amplitudeflow');
taf = struct('algorithm','taf');
raf = struct('algorithm','raf');
fienup = struct('algorithm','fienup');
gs = struct('algorithm','gerchbergsaxton');
cd = struct('algorithm','coordinatedescent','maxIters',60);
kac = struct('algorithm','kaczmarz','maxIters',60);
pmax = struct( 'algorithm','phasemax', 'maxIters',113);
plamp = struct('algorithm', 'phaselamp');  %maybe out of memory
plift = struct('algorithm','phaselift','maxIters',1000);  %maybe out of memory                                            


% Grab your pick of algorithms.
algorithms = {gs,fienup,wf,twf,rwf,ampflow,taf,raf,pmax,cd,kac};
algorithms_return = {'gs','fienup','wf','twf','rwf','ampflow','taf','raf','pmax','cd','kac'};
% Run benchmark
rec = benchmarkSynthetic(xitem, xvalues, yitem, algorithms, dataSet, params,opt);

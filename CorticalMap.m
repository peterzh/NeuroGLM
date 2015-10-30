%This script concerns creating maps of the parameters across the cortex

%% Define initialisations
model = 'C^N subset';
expRefs = {};

%% Load laser inactivation information (coord, on/off flag)
dat.expFilePath('Hopkins',
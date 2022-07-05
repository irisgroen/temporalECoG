% Figure 1-1 :
% Individual patient electrode coverages on brain mesh
%
% This script will generate one or more separate figures for each patient,
% showing different views of the implanted hemisphere(s)
%
% WARNING: running this script all at once with create a lot of figures (27
% in total). It is advised to run cell-by-cell for each patient in turn.
%
% 2022 Iris Groen

%close all;
atlasName           =  {'wang15_fplbl_norm'};
projectDir          = bidsRootPath; 

specs = [];
specs.plotelecs     = 'yes';
specs.plotlabel     = 'no';

%% Extended Data Figure 1_1: all patients with electrodes in final dataset:

%% patient 2
subject             = 'p02'; 
specs.plotmesh      = 'right';
specs.plotelecrad   = 2.3;
specs.plotlabel = 'yes';
% lateral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(-60,-10);

% medial view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(50,0);

% ventral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(180,-90);

%% patient 3
subject             = 'p03'; 
specs.plotmesh      = 'left';
specs.plotelecrad   = [];

% lateral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(-35,0);

% dorsal view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(20,30);

% ventral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(0,-90);

%% patient 4
subject             = 'p04'; 
specs.plotmesh      = 'right';
specs.plotelecrad   = [];

% lateral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(35,0);

% ventral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(0,-90);

%% patient 5
subject             = 'p05'; 
specs.plotmesh      = 'right';
specs.plotelecrad   = 2.3; % looks like sizes for second clinical grid are not correct?

% lateral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(30,0);

% dorsal view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(0,45);

% ventral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(0,-90);

%% patient 6
subject             = 'p06'; 
specs.plotmesh      = 'left';
specs.plotelecrad   = []; 

% lateral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(-50,5);

% medial view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(50,-10);

% ventral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(0,-90);

%% patient 7
subject             = 'p07'; 
specs.plotmesh      = 'left';
specs.plotelecrad   = []; 

% lateral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(-45,10);

% medial view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(35,0);

% ventral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(0,-90);

%% patient 10 (HD patient)
subject             = 'p10'; 
specs.plotmesh      = 'right';
specs.plotelecrad   = []; 

% lateral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(45,10);

% medial view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(-30,0);

% ventral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(0,-90);

% dorsal view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(0,70);

%% patient 11 (HD patient)
subject             = 'p11'; 
specs.plotmesh      = 'left';
specs.plotelecrad   = []; 

% lateral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(-45,-5);

% medial view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(50,-10);

% ventral view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(0,-90);

% dorsal view
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
view(0,45);

%% End of Figure 1_1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

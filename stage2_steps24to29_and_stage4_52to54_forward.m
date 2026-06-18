% Author: Paolo Belardinelli

% MRI Loading
mripath = '/Users/paolobelardinelli/Downloads/Scripts (1)/data_freesurfer_and_mri/S1or';
cd ./data_freesurfer_and_mri/S1or;
subjectname = 'S1ctf';

mri = ft_read_mri(fullfile(mripath,sprintf('%s.mgz',subjectname)));
mri = ft_determine_coordsys(mri, 'interactive', 'yes');

%% MRI check
cfg = [];

ft_sourceplot(cfg, mri);

print -dpng original_mri.png
%% Mri info
disp(mri)


%% Load Electrodes File
% filepath = '/Users/paolobelardinelli/Desktop/';
% filename = 'EEGMarkers20180615181822532.xml';
% elec = ft_read_sens(fullfile(filepath,filename));
 filepath = '/Users/paolobelardinelli/Downloads/';
filename = 'SPACETIME_020_20241017_prepost.txt';
 elec = ft_read_sens(fullfile(filepath,filename));
elec.elecpos(:,1)= -elec.elecpos(:,1);
elec.chanpos(:,1)= -elec.chanpos(:,1);
%% Segmentation: Labeling of Volumes 

cfg           = [];
cfg.coordsys  = 'neuromag';
cfg.output    = {'brain','skull','scalp'};
segmentedmri  = ft_volumesegment(cfg, mri);

save segmentedmri segmentedmri

disp(segmentedmri)

%% Segmentation: Boundary Mesh Creation 

cfg=[];
cfg.tissue={'brain','skull','scalp'};
cfg.numvertices = [2500 2000 1000]; %this can be changed 
bnd=ft_prepare_mesh(cfg,segmentedmri);

save bnd bnd

%% check meshes

figure;
ft_plot_mesh(bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(bnd(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(bnd(1),'edgecolor','none','facecolor',[0.4 0.6 0.4]);
%% check electrodes on scalp
disp(elec)
figure;

ft_plot_mesh(bnd(3), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
hold on;
ft_plot_sens(elec, 'coilshape', 'circle')

%% Just in case: realignment
cfg           = [];
cfg.method    = 'interactive';
cfg.elec      = elec;
cfg.headshape = bnd(3);
elec_realigned  = ft_electroderealign(cfg);
%% add toolboxes paths

TOOLBOXPATH = '/Users/paolobelardinelli/Desktop/toolboxes';
addpath(genpath(fullfile(TOOLBOXPATH, 'hbf_distribution_open_v170624')))
addpath(genpath(fullfile(TOOLBOXPATH, 'hbf_distribution_open_v170624')))
addpath(fullfile(TOOLBOXPATH, 'plotroutines_v170706'))
addpath(fullfile(TOOLBOXPATH, 'stenroos_functions'))
addpath(fullfile(TOOLBOXPATH, 'phastimate'))
addpath(fullfile(TOOLBOXPATH, 'freezeColors'))

%% load relevant data

subject = 'S1';

load(fullfile(mripath,sprintf('%s_sourcemodel_15684.mat',subject)));

headmodel = [];
headmodel.bmeshes(1).p = bnd(1).pos;% boundary meshes
headmodel.bmeshes(1).e = bnd(1).tri;
headmodel.bmeshes(2).p = bnd(2).pos;% boundary meshes
headmodel.bmeshes(2).e = bnd(2).tri;
headmodel.bmeshes(3).p = bnd(3).pos;% boundary meshes
headmodel.bmeshes(3).e = bnd(3).tri;
headmodel.smesh.p = sourcemodel.pos;% cortical surface mesh
headmodel.smesh.e = sourcemodel.tri;
headmodel.smesh.nn = CalcNodeNormals(headmodel.smesh); % direction of dipoles normal to surface mesh
headmodel.elec_ft = elec; % this satisfies http://www.fieldtriptoolbox.org/reference/ft_datatype_sens/
headmodel.elec_ft = ft_convert_units(headmodel.elec_ft, 'm');
headmodel.elec_ft.label(strcmp(headmodel.elec_ft.label, 'REF')) = {'FCz'};
headmodel.elec_ft.chanpos = headmodel.elec_ft.chanpos(4:end,:);
headmodel.elec_ft.elecpos = headmodel.elec_ft.elecpos(4:end,:);
headmodel.elec_ft.chantype = headmodel.elec_ft.chantype(4:end,:);
headmodel.elec_ft.chanunit = headmodel.elec_ft.chanunit(4:end,:);
headmodel.elec_ft.label = headmodel.elec_ft.label(4:end);
headmodel.refchannel = 'FCz';

figure
PlotMesh(headmodel.smesh);
hold on
PlotMesh(headmodel.bmeshes(1), 'facealpha', 1,'edgecolor', 'none', 'facecolor', [0.65 0.65 0.65])

% depending on the plotting results, the smesh and bmesh(1) could need a slight adjustment;

headmodel.bmeshes(1).p = 1.01*headmodel.bmeshes(1).p;% 


%% preprocessing

H{1} = headmodel.bmeshes(1);
H{2} = headmodel.bmeshes(2);
H{3} = headmodel.bmeshes(3);

headmodel.elec = hbf_ProjectElectrodesToScalp(elec_realigned.chanpos, H);
bmeshes3 = H; %3-shell model
ci = [1 1/80 1] *.33; %conductivities
co = [1/80 1 0] *.33;
D = hbf_BEMOperatorsPhi_LC(bmeshes3);
Tphi_full = hbf_TM_Phi_LC_ISA2(D, ci, co, 1);
Tphi_elecs = hbf_InterpolateTfullToElectrodes(Tphi_full, bmeshes3, headmodel.elec);
LFMphi_dir = hbf_LFM_LC(bmeshes3, Tphi_elecs, headmodel.smesh.p, headmodel.smesh.nn); % leadfield for normally oriented sources
LFMphi_xyz = hbf_LFM_LC(bmeshes3, Tphi_elecs, headmodel.smesh.p);

headmodel.leadfield = LFMphi_dir;
headmodel.label = elec_realigned.label;

clear( 'ci', 'co', 'D', 'Tphi_full', 'Tphi_elecs',  'LFMphi_xyz');

%% prepare layout

layout = ft_prepare_layout([], struct('elec', elec_realigned));

%% Sensitivity profiles

L = LFMphi_dir;
elec_aligned=headmodel.elec;
Eind=find(strcmp(headmodel.label,'C3'));


%%% Plot sensitivity profile
va=[-90 45];
cmap=jet(20);
cmap([1 2 19 20],:)=1;
set(figure(2),'position',[0 100 800 500]);clf;hold on
ltest=L(Eind,:);
ltest=ltest/max(abs(ltest));
PlotDataOnMesh(headmodel.smesh,ltest,'view',va,'colormap',cmap);

PlotPoints(elec_aligned.pproj,'k.',10);
PlotPoints(elec_aligned.pproj(Eind,:),'b.',20);
camzoom(1.3);

%% Visualize


fig = figure;
PlotMesh(bmeshes3{end-2}, 'figure', fig); % csf to bone
PlotMesh(bmeshes3{end-1}, 'facecolor', [248 232 192]./255, 'figure', fig); %bone to scalp
PlotMesh(bmeshes3{end}, 'figure', fig); % scalp to air
PlotMesh(headmodel.smesh, 'figure', fig);

save headmodel headmodel
%% Prepare Spatial FIlters
spf = [];
spf.c3_lap = zeros(1, length(headmodel.label));
spf.c3_lap(strcmp(headmodel.label, 'C3')) = 1;
spf.c3_lap(strcmp(headmodel.label, 'FC1')) = -0.25;
spf.c3_lap(strcmp(headmodel.label, 'FC5')) = -0.25;
spf.c3_lap(strcmp(headmodel.label, 'CP1')) = -0.25;
spf.c3_lap(strcmp(headmodel.label, 'CP5')) = -0.25;



%% Visualize Sensitivity Profiles and sensors

    load('mymap.mat')
    weights = spf.c3_lap;
    sens_prof = weights * headmodel.leadfield;
    sens_prof = sens_prof / max(abs(sens_prof)); %normalize
    
    SensProfPlot = PlotDataOnMesh(headmodel.smesh, sens_prof, 'view', [-90 45]);
    colormap(mymap), colorbar('off')
    
    freezeColors
    hold on
    elec_pos = headmodel.elec.pproj(weights ~= 0, :);
    elec_xyz = mat2cell(elec_pos, length(elec_pos), [1 1 1]); % actual positions
    plot3(elec_xyz{:}, 'Color', 'b', 'Marker', 'o', 'LineStyle', 'none', 'MarkerSize', 8)




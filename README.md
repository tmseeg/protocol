# Description of MATLAB Scripts

## stage2_steps12to18and 22_sourcespace.m

Derives a cortical source mesh by reading the NIfTI file into MATLAB with FieldTrip73. Dipoles can optionally be oriented normal to the surface of the mesh. Fiducials are interactively marked on the MR image, usually nasion, left and right preauricular points, as well as a generic positive z point.
Reads the NIfTI MR volume into MATLAB (```ft_read_mri```). Checks the coordinate system of the image (```ft_determine_coordsys```).
Defines fiducial reference points (i.e., nasion, left preauricular point, right preauricular point and a generic z>0 point) in the sensor coordinate system (```ft_volumerealign```).
Reslices MR volume to obtain cubic voxels (```ft_volumereslice```) and check that the whole head is still contained in the 3D MR image (```ft_determine_coordsys```). If not, repeats operations from Step 13. This step can be problematic but it is necessary since later FreeSurfer will make the voxels cubic anyway but we need to derive the transformation matrix from the original NIfTI coordinate system to sensor space and to Anterior-Commissure-Posterior Commissure (acpc) space and, therefore, need consistency from the beginning with the image dimension.
1.	Writes the first ```.mgz``` MR volume in sensor coordinates for further Boundary Element Method (BEM) segmentation consistent with the source mesh we are creating here. The sensor space-based head model is where we will source-extract cortical activity (```ft_volumewrite```).
1.	Defines the acpc allowing FreeSurfer to determine an origin point from which to start the segmentation procedure, together with a generic ‘right hemisphere point’ to inform FreeSurfer (```ft_volumerealign```).
1.	Writes the second ```.mgz``` volume in acpc coordinates for FreeSurfer white and grey matter segmentation (```ft_volumewrite```).

### After Freesurfer recon-all execution from Bash

Merges left and right hemisphere meshes into one combined mesh using FieldTrip and co-register the resulting source-space mesh to the sensor coordinate system: the two left and right decimated cortical meshes are GIfTI (```.gii```) files, one for each hemisphere. The coordinates of the vertices are still defined in the acpc coordinate system. The two meshes can be merged by means of ```ft_read_headshape``` (geometrically, this is a juxtaposition of the two meshes, which are already spatially contiguous). Finally, the coordinates of the mesh points can be transformed into the sensor space leveraging the two space-transformation matrices generated in Steps 16 and 17. Moreover, a now useless field is removed and the mesh is saved to disk. The mesh is then plotted to check for anomalies and inaccuracies.

## stage2steps24to29_and_stage4_52to54_forward.m

Extracts head anatomical boundary. Loads the ```.mgz``` volume file generated with sourcespace.m in sensor coordinates. 
The quality of this image can be further visually inspected through ```ft_sourceplot```.
 3-shells (scalp/outer skull, internal skull, brain envelope) are segmented out of the MR volume. Volumes external to the source space are first labeled (as concentric volumes: scalp, inner skull, brain envelope) and then segmented as triangulated meshes with SPM through FieldTrip functions (```ft_volumesegment```) to label each voxel belonging to the three volumes, and ```ft_prepare_mesh``` to generate the mesh comprising of vertices and edges (or triangles)). Finally, plot the meshes to check their correctness, e.g., exclude the presence of holes in the meshes (```ft_plot_mesh```).
It aligns source mesh within head meshes. In case the source mesh exceeds the internal BEM mesh (brain envelope) in one or more gyri, it may be necessary either to slightly translate or adjust either the brain envelope or the source mesh itself. 
Validation and plotting: all meshes need to be concentric, otherwise the leadfield function will not generate a correct leadfield.
Finally, it saves the head model: source meshes, boundary meshes and leadfield are saved (the meshes won’t change throughout different sessions but the EEG electrode positions and the leadfields will).
It is also possible to check consistency of sensors with the sensitivity profile at source, to check that coordinates of sensor and source space are aligned correctly.

## stage3_step048_create_epochs.m

This files loads raw EEG data in Bittium NeurOne format from the ```data_neurone``` folder. EEG channels and EMG channels are stored separately. Data is epoched around the TMS events (identified by an ```{'A - Out'}``` trigger marker). MATLAB is required to run this script. The required EEGLAB toolbox with the Bittium NeurOne file format loading extension is provided in the toolboxes folder.

Two variants of epoched data are generated for use in subsequent steps:

-	Pre-stimulus data ([-1.5 -0.005]) in fieldtrip format, downsampled to 1kHz, containing all channels, saved in the ```data_epoched``` folder as ```OPENLOOP_ftdata.mat```.
-	Peri-stimulus data ([-1.5 1.5])  in EEGLAB format, containing only channels that are part of the spatial filter (```{'C3', 'FC1', 'FC5', 'CP1', 'CP5'}```), saved in the ```data_epoched``` folder as OPENLOOP_eeglab.mat.

## stage3_step049_estimate_spectrum.m

A power spectrum is estimated from the epoched data using fieldtrip toolbox. MATLAB is required to run this script. The required fieldtrip toolbox is provided in the toolboxes folder.
The script illustrates how to apply a specific single-channel spatial filter montage and alternatively, how to transform all channels using a Surface Laplacian. The aperiodic 1/f fractal background part of the spectrum is estimated using IRASA (using the fieldtrip toolbox).

## stage4_53to57_sourceVSsensor.m

The Leadfield matrix (dimensions: number of dipoles x number of sensors) is calculated as the transfer matrix between source dipoles and sensors. Since electrical conductivities in the different compartments of the head are very different (for example, the skull is a bad conductor), a three-compartment volume conductor model is applied using the three meshes generated in Stage 2 and the solution Boundary Element Method is provided considering three homogeneous and concentric layers: intracranial space (conductivity 0.33 S/m), skull (0.0041 S/m), and scalp (0.33 S/m).
Then data covariance is estimated from EEG prior to TMS pulse.
Inverse solution in time domain from forward model and covariance matrix is obtained using LCMV beamforming.
One can interactively select the cortical area from which to reconstruct the source activity through the ```mf_selectdipole``` function.
Finally, power spectra for sensor Hjorth, M1 and S1 handknob can be calculated.

## stage5_step075_simulate_phase_accuracy.m

This script loads epoched data and applies a spatial filter. The phase at the centre of each data segment is determined using standard acausal (i.e., considering data before and after the time point of interest) signal processing methods (band-pass filter followed by Hilbert transform). This is possible only in the absence of a stimulus artifact, which is why this script uses pre-stimulus or resting-state data. For each segment, the script also estimates the phase at the centre of the epoch using the same causal method (i.e., using the preceding data only) that is implemented by the real-time algorithm. The difference between the two estimates is considered the phase error and is plotted as a circular histogram.

## stage6_step099_quality_control.m

This script visualizes the data recorded in an EEG-triggered TMS measurement. It loads epoched peri-stimulus data in EEGLAB format. A spatial filter is applied. A spectrum is estimated from the pre-stimulus data.
The following panels are generated:
-	Timing of triggers (stimulus number vs. time) with a stimulus interval distribution
-	Average TMS-evoked EEG potential in original reference format
-	A selection of 10 single trials are plotted
-	The average pre-stimulus spatially filtered signal is plotted (this should show in oscillation in phase-triggered experiments and no oscillation in a random-phase or open-loop control experiment.)
-	The full and 1/f-corrected power spectral density is shown
-	Average TMS-pulse related EMG data (showing MEP in supra-threshold TMS)

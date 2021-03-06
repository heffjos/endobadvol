clear
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The folder that contains your subject folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Exp = '/nfs/turbo/berent-lab/metabolic/Endopoid/Data';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Path where your logfiles will be stored
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LogTemplate = './Logs';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Path where your images are located
%%%
%%%  Variables you can use in your template are:
%%%       Exp      = path to your experiment directory
%%%       Subject  = name of subject from SubjDir (using iSubject as index of row)
%%%       Run      = name of run from RunDir (using iRun as index of row)
%%% Examples:
%%% ImageTemplate = '[Exp]/Subjects/[Subject]/func/[Run]/';
%%% ImageTemplate = '[Exp]/Subjects/[Subject]/TASK/func/[Run]/'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImagePathTemplate = '[Exp]/[Subject]/func/visual/[Run]/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% spm_select filter used to collect subject images
%%% Type help regexp and spm_filter for description how to create filters
%%% generally it will always be something like '^sw3mm_ra_spm8_run.*nii'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BaseFileSpmFilter = '^s888_w2mm_rtrun.*nii';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A list of run folders where the script can find the images to use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RunDir = {
'run_01';
'run_02';
'run_03';
'run_04';
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The list of subjects to process
%%% The format is 'subjectfolder',[runs to include]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SubjDir = {
% 'abs13end00001_00568', [1 2 3];
% % 'abs13end00002_00565', [NA NA NA];
% 'abs13end00053_04228', [1 2 3];
% % 'abs13end00056_04375', [NA NA NA];
% 'abs13end00246_04277', [1 2 3];
% 'abs13ins00003_01346', [1 2 3];
% 'abs13ins00004_01411', [1 2 3];
% 'abs13ins00005_01533', [1 2 3];
% 'abs13ins00006_01591', [1 2 3];
% 'abs13ins00007_01590', [1 2 3];
% 'abs13ins00008_01659', [1 2 3];
% 'abs13ins00009_01770', [1 2 3];
% 'abs13ins00010_01687', [1 2 3];
% 'abs13ins00011_01690', [1 2 3];
% 'abs13ins00013_01841', [1 2 3];
% 'abs13ins00014_02465', [1 2 3];
% 'abs13ins00015_02327', [1 2 3];
% % 'abs13ins00016_02304', [NA NA NA]; % missing t1spgr
% 'abs13ins00017_02371', [1 2 3];
% 'abs13ins00020_02581', [1 2 3];
% 'abs13ins00021_02495', [1 2 3];
% 'abs13ins00024_03288', [1 2 3];
% 'abs13ins00025_02644', [1 2 3];
% 'abs13ins00026_02659', [1 2 3];
% 'abs13ins00028_02674', [1 2 3];
% 'abs13ins00035_03073', [1 2 3];
% 'abs13ins00036_02976', [1 2 3];
% 'abs13ins00037_02988', [1 2 3];
% 'abs13ins00038_03123', [1 2 3];
% 'abs13ins00039_03046', [1 2 3];
% 'abs13ins00040_03279', [1 2 3];
% 'abs13ins00042_03392', [1 2 3];
% 'abs13ins00046_03742', [1 2 3];
% 'abs13ins00049_03933', [1 2 3];
% 'abs13ins20017_03308', [1 2 3];
'abs13ins20024_04204', [1 2 3];
% 'abs13ins20036_03982', [1 2 3];
% 'abs13ins20042_04134', [1 2 3];
};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Location of masterdata CSV file
%%%
%%%  Variables you can use in your template are:
%%%       Exp            = path to your experiment directory
%%%        *             = wildcard (can only be placed in final part of template)
%%% Examples:
%%% MasterTemplate='[Exp]/Scripts/MasterData/EurekaDM_Master.csv';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MasterDataFilePath = '/nfs/turbo/berent-lab/metabolic/Endopoid/endoeprime/MasterDataFiles/MDF_Visual.csv';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Column number in the MasterData file where subject numbers are located
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SubjColumn = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Column number in the MasterData file where run numbers are located
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RunColumn = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Column number(s) in the MasterData file where conditions numbers are located
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CondColumn = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Column number in the MasterData file where your Onset times are located
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TimeColumn = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Column number in the MasterData file where your Durations are located
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DurationColumn = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% List of conditions in your model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ConditionName = {
'Match';
'Delay1';
'Delay4';
};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% If you are including any Parametric regressors in your model
%%% syntax: 'Parameter Name', Master Data File column, Condition Number as listed
%%%         in ConditionName for parametric regressor, polynomial order
%%%
%%% If using multiple condition columns above, the parametric regressor
%%% will be used with every condition number as given in the third column.
%%%
%%% NOTE: If you want to include parametric regressors for a condition, the 
%%% values of your regressor MUST change over trials.  SPM can not include
%%% a constant parametric regressor and it will cause problems with contrasts.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParametricList = {

};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Location of user specified regressor files. To use these make sure 
%%%  column 2 of RegOp is set to 1 above.
%%%
%%%  syntax: File location, number of regressors to use OR or a vector specifying
%%%          columns to include, derivative order to include for regressor, 
%%%          polynomial term (applies to derivatives as well if being used)
%%%
%%%  Variables you can use in your template are:
%%%       Exp         = path to your experiment directory
%%%       Subject     = folder name of current subject
%%%       Run         = folder name of current run
%%%        *          = wildcard (can only be placed in final part of template)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RegFilesTemplate = {
    '[Exp]/[Subject]/func/visual/[Run]/realign.txt',Inf,0,1;
    '/nfs/turbo/berent-lab/Researchers/heffjos/BadVolumeComparison/[Subject]/func/visual/[Run]/CensorVectors.csv', Inf, 0, 1;
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% List of contrasts to add to the estimated model
%%% Format is 'Name of contrast' [Cond1 Param1...N]...[CondN Param1...N] [Reg1...RegN]
%%% You need to properly balance/weight your contrasts below as if it was just one run/session
%%% The script will handle balancing it across runs. To mark all values for [Reg1...RegN] as 0,
%%% set the last field as '[]' without quotes.
%%% If using fir basis functions, condition vectors must be at least as
%%% long as the number of bins being used.
%%%
%%% If you are using FIR basis functions, you can automatically generate
%%% contrasts for each regressor using the FirDoContrasts flag  in the advanced
%%% section.  Contrasts listed here are still generated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ContrastList = {
    'Match-Delay1' 1 -1 0 [];
    'Match-Delay4' 1 0 -1 [];
    'Delay1-Match' -1 1 0 [];
    'Delay4-Match' -1 0 1 [];
    'Delay1-Delay4' 0 1 -1 [];
    'Delay4-Delay1' 0 -1 1 [];
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Path for output images
%%%
%%%  Variables you can use in your template are:
%%%       Exp        = path to your experiment directory
%%%       Subject    = name of subject from SubjDir (using iSubject as index of row)
%%%       Run        = name of run from RunDir (using iRun as index of row)
%%% Examples:
%%% OutputTemplate = '[Exp]/Subjects/[Subject]/func/[Run]/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OutputDir = fullfile(pwd, './FL_CensorOnly/[Subject]/Visual');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The TR your data was collected at
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TR = 2;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~ Advanced ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run mode
%%% 1 = regular
%%% 2 = Contrast add on
%%% 3 = test without running anything
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mode = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set the contrast start point, used only if Mode = 2
%%% 1 = Overwrite Previous Contrasts
%%% 2 = Append new contrasts to previous ones 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StartOp=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analyze data on a local drive on the machine you run the script on
%%% If your data is located on /net/data4 or a similar network drive, using
%%% this option will greatly reduce the required processing time.
%%% IMPORTANT NOTE: Due to the method of sandboxing, using this WILL
%%% OVERWRITE existing results without prompting you, so please be sure
%%% your paths are all correct before running.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UseSandbox = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ref point for data out of 16, use same fraction as ref slice for slice timing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fMRI_T0 = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CondThreshold
%%% 0 = only remove empty conditions
%%% 1 = remove singleton conditions (useful b/c SPM won't estimate a beta for
%%% parameters that modulate a singleton condition)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ConditionThreshold = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This sets whether the models are the same for everyone or subject-specific
%%%   0 - each person has different models, grab data from section of MasterData based on subject index
%%%   1 - each person has identical models, grab all from first block of MasterData
%%%   NOTE: Regressors still use subject index so are not identical across subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IdenticalModels = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Total number of trials in your experiment for a subject (only used if IdenticalModels = 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TotalTrials = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Basis function to use in model.  Valid options are 'hrf' or 'fir'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Basis = 'hrf';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Only used if Basis = 'hrf'
%%% Indicates wheter to model the derivates of the canonical hrf function.
%%% 0 - none
%%% 1 - derivative
%%% 2 = derivative and dispersion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HrfDerivative = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Only used if Basis = 'fir'
%%% Indicates the poststimulus duration for the hemodynamic response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FirDuration = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Only used if Basis = 'fir'
%%% Indicates the bins to use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FirBins = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Only used if Basis = 'fir'
%%% Flag used to automatically generate contrasts when FIR basis functions
%%% are used.  The contrasts generated have a 1 at each FIR regressor
%%% and zero everywher else.
%%% Contrasts listed in ContrastList WILL still be generated
%%% 0 - create no contrasts
%%% 1 - create a contrast for each bin for all conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FirDoContrasts = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run-specific contrasts (NumContrasts rows each of NumRuns elements)
%%% This allows you to set up a list of weights for each run for each
%%% contrast. This combines with the ContrastList above. Leaving a
%%% particular row empty ([]) is equivalent to setting the contrast weights
%%% as 1 for each run (i.e. no change from the standard contrast method).
%%% Example:
%%% If you have 4 runs and want to compare a condition in run 1 against a
%%% condition in run 4 you would set the contrast weight for that contrast
%%% as [1 0 0 -1].  If you want to only look at the contrast averaged
%%% across runs 1 and 2 you would set the weights as [1 1 0 0]. If you just
%%% want the standard contrast (the average across all runs) you can set it
%%% either to [] or [1 1 1 1].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ContrastRunWeights = {
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Volumespecifier specifies the number of volumes to use for each run.
%%% The number of columns must equal the number of runs specified in 
%%% RunDir.  The top rows indicate from which volume to start and the 
%%% bottom row indicates which volume to stop.  For instance,
%%% VolumeSpecifier = [10 20; 130 150] will start from the 10th volume
%%% and stop at the 130th volume for the first run in RunDir.  The second
%%% run in RunDir will start at the 20th volume and stop at the 150th
%%% volume.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VolumeSpecifier = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%ScaleOp - 'Scaling' = do proportonal scaling
%%%          'None' = do standard grand mean scaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ScaleOp = 'None';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Path and name of explicit mask to use at first level.
%%% Leave this blank ('') to turn off explicit masking
%%%
%%%  Variables you can use in your template are:
%%%       Exp      = path to your experiment directory
%%%       Subject  = name of subject from SubjDir (using iSubject as index of row)
%%% Examples:
%%% ExplicitMask = '[Exp]/Subjects/[Subject]/mask.nii';
%%% ExplicitMask = '/MethodsCore/ConnTool/Templates/3mm_mask.nii'
%%% ExplicitMask = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExplicitMask = '/nfs/turbo/berent-lab/metabolic/Endopoid/Data/[Subject]/func/visual/VisualMask.nii';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Use AR(1) auto-regression correction or not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
usear1 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% High pass filter cutoff value in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hpf = 128;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The main default that impacts first level analysis is the implicit 
%%% masking threshold. The default is 0.8 which means that voxels that have
%%% a value of less than 80% of the grand mean will be masked out of the
%%% analysis.  This default value can be problematic in some susceptibility
%%% prone areas like OFC.  A more liberal value like 0.5 can help to keep
%%% these regions in the analysis.  If you set this value very low, you'll
%%% want to use an explicit mask to exclude non-brain regions from
%%% analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spmdefaults = {
    'mask.thresh'   0.0001;
};

global mcRoot;

mcRoot = '/usr/local/MethodsCore'

addpath(fullfile(mcRoot,'matlabScripts'))
addpath(fullfile(mcRoot,'FirstLevel'))
addpath(fullfile(mcRoot,'FirstLevel','functions'));
addpath(fullfile(mcRoot,'SPM','SPM8','spm8_with_R6313'))

FirstLevel_mc_central


Original = '/nfs/turbo/berent-lab/metabolic/Endopoid/Analyses/FirstLevel/abs13ins20024_04204/';
CensorOnly = fullfile(pwd, 'FL_CensorOnly/abs13ins20024_04204');

VerbalMask = '/nfs/turbo/berent-lab/metabolic/Endopoid/Data/abs13ins20024_04204/func/verbal/VerbalMask.nii';
VisualMask = '/nfs/turbo/berent-lab/metabolic/Endopoid/Data/abs13ins20024_04204/func/visual/VisualMask.nii';

% let's do verbal first
nii = nifti(VerbalMask);
VerbalMask = nii.dat(:,:,:,:) == 1;
[X Y Slices] = size(VerbalMask);
VerbalMask = reshape(VerbalMask, X*Y, Slices);

VerbalSim = zeros(Slices, 2); % columns = number of contrasts
for i = 1:2
    OrigNii = nifti(fullfile(Original, sprintf('Verbal/con_%04d.hdr', i)));
    OrigData = OrigNii.dat(:,:,:,:);
    OrigData = reshape(OrigData, X*Y, Slices);

    CensorNii = nifti(fullfile(CensorOnly, sprintf('Verbal/con_%04d.hdr', i)));
    CensorData = CensorNii.dat(:,:,:,:);
    CensorData = reshape(CensorData, X*Y, Slices);

    for k = 1:Slices
        VerbalSim(k, i) = corr(OrigData(VerbalMask(:, k), k), CensorData(VerbalMask(:, k), k));
    end
end
    
    
% now let's do visual
nii = nifti(VisualMask);
VisualMask = nii.dat(:,:,:,:) == 1;
[X Y Slices] = size(VisualMask);
VisualMask = reshape(VisualMask, X*Y, Slices);

VisualSim = zeros(Slices, 2); % columns = number of contrasts
for i = 1:6
    OrigNii = nifti(fullfile(Original, sprintf('Visual/con_%04d.hdr', i)));
    OrigData = OrigNii.dat(:,:,:,:);
    OrigData = reshape(OrigData, X*Y, Slices);

    CensorNii = nifti(fullfile(CensorOnly, sprintf('Visual/con_%04d.hdr', i)));
    CensorData = CensorNii.dat(:,:,:,:);
    CensorData = reshape(CensorData, X*Y, Slices);

    for k = 1:Slices
        VisualSim(k, i) = corr(OrigData(VisualMask(:, k), k), CensorData(VisualMask(:, k), k));
    end
end


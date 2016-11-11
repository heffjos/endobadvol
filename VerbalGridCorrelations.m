clear all;

Original = '/nfs/turbo/berent-lab/metabolic/Endopoid/Analyses/FirstLevel/abs13ins20024_04204/';
CensorOnly = fullfile(pwd, 'FL_CensorOnly/abs13ins20024_04204');

VerbalMask = '/nfs/turbo/berent-lab/metabolic/Endopoid/Data/abs13ins20024_04204/func/verbal/VerbalMask.nii';

% calculate verbal now
nii = nifti(VerbalMask);
Mask = nii.dat(:,:,:,:);
[X Y Slices] = size(Mask);
Vol = zeros(X, Y, Slices);

XIndex = [2:3:X]';
YIndex = 2:3:Y;
ZIndex = zeros(1, 1, length([2:3:Slices]));
ZIndex(:) = 2:3:Slices;

XCoord = repmat(XIndex, 1, length(YIndex), length(ZIndex));
YCoord = repmat(YIndex, length(XIndex), 1, length(ZIndex));
ZCoord = repmat(ZIndex, length(XIndex), length(YIndex), 1);
Idx = sub2ind(size(Vol), XCoord, YCoord, ZCoord);

% remove centers outside mask
Vol(Idx) = 1;
Idx = find((Vol .* Mask) == 1);
[XCoord YCoord ZCoord] = ind2sub(size(Vol), Idx);
Vol = Vol .* 0;

% create 3x3x3 touching boxes
TmpXGrid = repmat([1:X]', 1, Y, Slices);
TmpYGrid = repmat(1:Y, X, 1, Slices);
TmpZGrid = zeros(1, 1, Slices);
TmpZGrid(:) = 1:Slices;
TmpZGrid = repmat(TmpZGrid, X, Y, 1);

XIndex = (TmpXGrid > 0) .* (TmpXGrid < 4);
YIndex = (TmpYGrid > 0) .* (TmpYGrid < 4);
ZIndex = (TmpZGrid > 0) .* (TmpZGrid < 4);
TmpIndex = XIndex .* YIndex .* ZIndex;
ReferenceBox = find(TmpIndex);
[RX RY RZ] = ind2sub(size(Vol), ReferenceBox);
IndOffsets = [RX-2 RY-2 RZ-2];

% create volume of ROIs for posterity
fprintf(1, 'Creating ROI volumes.\n');

for i = 1:length(XCoord)
    XIdx = XCoord(i) + IndOffsets(:, 1);
    YIdx = YCoord(i) + IndOffsets(:, 2);
    ZIdx = ZCoord(i) + IndOffsets(:, 3);

    ToKeep = logical( (XIdx <= X) .* (YIdx <= Y) .* (ZIdx <= Slices) );
    Idx = sub2ind(size(Vol), XIdx(ToKeep), YIdx(ToKeep), ZIdx(ToKeep));
    Idx( Mask(Idx) == 0 ) = [];
    
    if any(Vol(Idx))
        fprintf(1, 'Overlap Verbal index: (%d, %d, %d)\n', ...
            XCoord(i), YCoord(i), ZCoord(i));
        break;
    end
    Vol(Idx) = i;
end
V = spm_vol(VerbalMask);
V.fname = 'Output/Masks/VerbalRoi3.nii';
V.dt(1) = 16; % float32
V.pinfo(3) = 352;
spm_write_vol(V, Vol);

% calculate verbal correlations now
fprintf(1, 'Calculate verbal correlations now.\n');

UniqueCon = [1];
VerbalSim3 = zeros(length(XCoord), length(UniqueCon)); % columns = number of contrasts
for iCon = 1:length(UniqueCon)
    OrigNii = nifti(fullfile(Original, sprintf('Verbal/con_%04d.hdr', UniqueCon(iCon))));
    OrigData = OrigNii.dat(:,:,:,:);
    
    CensorNii = nifti(fullfile(CensorOnly, sprintf('Verbal/con_%04d.hdr', UniqueCon(iCon))));
    CensorData = CensorNii.dat(:,:,:,:);

    for iCoord = 1:length(XCoord)
        XIdx = XCoord(iCoord) + IndOffsets(:, 1);
        YIdx = YCoord(iCoord) + IndOffsets(:, 2);
        ZIdx = ZCoord(iCoord) + IndOffsets(:, 3);

        ToKeep = logical( (XIdx <= X) .* (YIdx <= Y) .* (ZIdx <= Slices) );
        Idx = sub2ind(size(Vol), XIdx(ToKeep), YIdx(ToKeep), ZIdx(ToKeep));
        Idx( Mask(Idx) == 0 ) = [];
    
        VerbalSim3(iCoord, iCon) = corr(OrigData(Idx), CensorData(Idx));
    end
end

% create low correlation masks
fprintf(1, 'Creating low correlation masks.\n');

Vol = Vol .* 0;
Vol1 = Vol;
Vol2 = Vol;
Vol3 = Vol;
for i = 1:size(VerbalSim3, 1)
    if VerbalSim3(i, 1) < 0.9
        XIdx = XCoord(i) + IndOffsets(:, 1);
        YIdx = YCoord(i) + IndOffsets(:, 2);
        ZIdx = ZCoord(i) + IndOffsets(:, 3);

        ToKeep = logical( (XIdx <= X) .* (YIdx <= Y) .* (ZIdx <= Slices) );
        Idx = sub2ind(size(Vol), XIdx(ToKeep), YIdx(ToKeep), ZIdx(ToKeep));
        Idx( Mask(Idx) == 0 ) = [];

        Vol1(Idx) = 1;
    end

    if VerbalSim3(i, 1) < 0.8
        XIdx = XCoord(i) + IndOffsets(:, 1);
        YIdx = YCoord(i) + IndOffsets(:, 2);
        ZIdx = ZCoord(i) + IndOffsets(:, 3);

        ToKeep = logical( (XIdx <= X) .* (YIdx <= Y) .* (ZIdx <= Slices) );
        Idx = sub2ind(size(Vol), XIdx(ToKeep), YIdx(ToKeep), ZIdx(ToKeep));
        Idx( Mask(Idx) == 0 ) = [];

        Vol2(Idx) = 1;
    end

    if VerbalSim3(i, 1) < 0.7
        XIdx = XCoord(i) + IndOffsets(:, 1);
        YIdx = YCoord(i) + IndOffsets(:, 2);
        ZIdx = ZCoord(i) + IndOffsets(:, 3);

        ToKeep = logical( (XIdx <= X) .* (YIdx <= Y) .* (ZIdx <= Slices) );
        Idx = sub2ind(size(Vol), XIdx(ToKeep), YIdx(ToKeep), ZIdx(ToKeep));
        Idx( Mask(Idx) == 0 ) = [];

        Vol3(Idx) = 1;
    end
end
V = spm_vol(VerbalMask);
V.fname = 'VerbalCon1Thr0.9.nii';
V.dt(1) = 16; % float32
V.pinfo(3) = 352;
spm_write_vol(V, Vol1);

V.fname = 'VerbalCon1Thr0.8.nii';
spm_write_vol(V, Vol2);

V.fname = 'VerbalCon1Thr0.7.nii';
spm_write_vol(V, Vol3);

% now write out correlation values
fprintf(1, 'Writing correlation values.\n');

fid = fopen('VerbalCorr.csv', 'w');
fprintf(fid, ['Participant,' ...
    'Task,' ...
    'ContrastNum,' ...
    'X,' ...
    'Y,' ...
    'Z,' ...
    'Correlation\n']);
for iCon = 1:size(VerbalSim3, 2)
    for iRow = 1:size(VerbalSim3, 1)
        fprintf(fid, 'abs13ins20024_04204,');
        fprintf(fid, 'Verbal,');
        fprintf(fid, 'Con%d,', UniqueCon(iCon));
        fprintf(fid, '%d,', XCoord(iRow));
        fprintf(fid, '%d,', YCoord(iRow));
        fprintf(fid, '%d,', ZCoord(iRow));
        fprintf(fid, '%0.2f\n', VerbalSim3(iRow, iCon));
    end
end
fclose(fid);

fprintf(1, 'All Done!\n');
    

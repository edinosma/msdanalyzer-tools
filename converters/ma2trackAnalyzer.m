% Convert msdanalyzer tracks to a format Track Analysis understands. Written by Edin Osmanbasic on 13.6.23 under GPLv3 license.
%% Declare variables
ma = importdata(""); % Path to msdanalyzer variable
SAVE_FILE = false;

%% Calc variables
lTracks = cellfun('size', ma{1}.tracks, 1);
offset  = cumsum(lTracks);

%% Move data
% Init
tr = nan(sum(lTracks), 4);
R = ma{1}.tracks{1}(:,2:end);
frames = 1:lTracks(1);
tr(1:offset(1),:) = [R, frames', ones(lTracks(1), 1)];

% Loop
for i = 2:length(lTracks)
    R = ma{1}.tracks{i}(:,2:end);
    frames = 1:lTracks(i);
    colStart = offset(i-1) + 1;
    tr(colStart:offset(i), :) = [R, frames', repmat(i, lTracks(i), 1)];
end

%% Save
if SAVE_FILE
    save("convertedTracks.mat", "tr", "-v7.3")
end

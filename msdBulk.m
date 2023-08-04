% Calculate MSDs of Troika trjR variables in bulk. Written by Edin Osmanbasic on 4/3/23.
% DOI for citing msdanalyzer: 10.1083/jcb.201307172
clearvars -except trjR;

%% Declare variables
[filename,filepath] = uigetfile('*.mat','SPT files?','MultiSelect','on');

TIME_LAG = 1;      % In units of seconds.
PX_SIZE  = 0.1067; % In units of spatial unit (nm) per px.
                   % For simulated data, disable this by setting to 1.

    ALPHA_CLIP = 0.5; % What percentage of data do we want to find alpha from?
DIFFUSION_CLIP =   1; % What percentage of data do we want to find diffusions from?

SAVE_DIR = []; % Where to save msdanalyzer vars? Use [] to not save. 

%% Convert lone file names from strings to cells
if (ischar(filename))
    filename = {filename};
end

%% Allocate memory
N_FILES = 1:length(filename);
[ma2{N_FILES}] = deal(msdanalyzer(2, "nm", "s"));
trjRs = cell(1,N_FILES(end));

parfor ifile = N_FILES
    disp(string(ifile) + '/' + string(max(N_FILES)) + ' trjRs loaded.');
    thisfile = filename{ifile};
    trjRs{ifile} = importdata([filepath thisfile]);
end

N_FRAMES = cellfun('size',trjRs,1);
N_PARTICLES = cellfun('size',trjRs,3);
msd_trjR = cell(N_FILES(end), max(N_PARTICLES));

%% Loop: Convert Troika trjR to msdanalyzer trjR
% trjR is stored as [y x w]. We want [t x y].
for iConvert = N_FILES
    n_frames = (1:N_FRAMES(iConvert))';
    
    for i = 1:N_PARTICLES(iConvert)
        msd_trjR{iConvert,i} = [(n_frames .* TIME_LAG) ...          % t
                            (trjRs{iConvert}(:,2,i) .* PX_SIZE) ... % x
                            (trjRs{iConvert}(:,1,i) .* PX_SIZE)];   % y, as [t x y]

        % Remove (0,0)s from trjR
        msd_zeros = (msd_trjR{iConvert,i}(:,2) == 0); % x indicies (0 in x is 0 in y)    
        msd_trjR{iConvert,i}([msd_zeros msd_zeros msd_zeros]) = [];
        msd_trjR{iConvert,i} = reshape(msd_trjR{iConvert,i},[],3); % reshape array to X by 3
    end
end

%% Calc MSD
ticBytes(gcp); tic;
parfor iMSD = N_FILES % This may take a little bit!
    disp('Calculating MSD (and fits) ' + string(iMSD) + ' out of ' + string(max(N_FILES)));

    % Only take in cells with data (remove padding)
    ma2{iMSD} = ma2{iMSD}.addAll(msd_trjR(iMSD, 1:N_PARTICLES(iMSD)));
    
    ma2{iMSD} = ma2{iMSD}.computeMSD; % MSD
    ma2{iMSD} = ma2{iMSD}.fitMSD; % fit MSD for diffusion
    ma2{iMSD} = ma2{iMSD}.fitLogLogMSD(ALPHA_CLIP); % fit log MSD for alpha
end
toc; tocBytes(gcp);

%% Save data to file "ma_filename.mat"
if ~isempty(SAVE_DIR)
    for iSaving = N_FILES
        ma = ma2{iSaving};
        save(SAVE_DIR + 'ma_' + filename{iSaving},'ma','-v7.3');
    end
end
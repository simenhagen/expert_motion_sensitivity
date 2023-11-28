%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Motion study
%           
%  Hagen, Vuong, Scott, Curran, Tanaka, 2021, JOV
%
%  Script by Hagen and Vuong 



clear all; clc;

addpath(fullfile(pwd,'stair')); %%% QCV %%%


% Debug flag. 1 = debugging
DEBUG = 0;


if DEBUG
    Screen('Preference',    'SkipSyncTests', 1);
end
%Screen('Preference', 'ConserveVRAM', 64);
%[keyboardIndices, productNames, allInfos] = GetKeyboardIndices('Dell USB Entry Keyboard');

%set up variables for staircase  (1-up/3-down --> ~79% acc)  (1-up/2-down --> ~70%)
startndots = [ 4+randi(5,1)  10+randi(5,1) ];
initial_stepsize = [ 4 3 ];   
nreps_init_stepsize = [ 2 2 ]; % four initial reversals not counted
final_stepsize = 3;
rules = [ 3 1 12 ]; % 3-correct -> up, 1-wrong -> down, n reversals to stop
stair{1} = STAIR_init(startndots(1),[ initial_stepsize; nreps_init_stepsize ],final_stepsize,rules); % birds upright
stair{2} = STAIR_init(startndots(2),[ initial_stepsize; nreps_init_stepsize ],final_stepsize,rules); % birds upright
stair{3} = STAIR_init(startndots(1),[ initial_stepsize; nreps_init_stepsize ],final_stepsize,rules); % birds inverted
stair{4} = STAIR_init(startndots(2),[ initial_stepsize; nreps_init_stepsize ],final_stepsize,rules); % birds inverted
stair{5} = STAIR_init(startndots(1),[ initial_stepsize; nreps_init_stepsize ],final_stepsize,rules); % human upright
stair{6} = STAIR_init(startndots(2),[ initial_stepsize; nreps_init_stepsize ],final_stepsize,rules); % human upright
stair{7} = STAIR_init(startndots(1),[ initial_stepsize; nreps_init_stepsize ],final_stepsize,rules); % human inverted
stair{8} = STAIR_init(startndots(2),[ initial_stepsize; nreps_init_stepsize ],final_stepsize,rules); % human inverted
nstairs = 8;

nconditions = nstairs/2;  % (birds, humans; upright, inverted)

trial_conds_prac = [];
for h = 1:2 % stimuli class (bird, human)
    for i = 1:2 % orientation (upright, inverted)
        for j = 1:2 % 1 = intact, 2 = scrambled
            for k = 1:7 % stimulus index
                 for r = 1:1 % number of reps
                    trial_conds_prac = [ trial_conds_prac; h i 1 j k ];
                end
            end
        end
    end
end

totaltrials_prac_humans = length(trial_conds_prac);
ord = randperm(size(trial_conds_prac,1));
trial_conds_prac = trial_conds_prac(ord,:);

trial_conds_exp = [];
for i = 1:nstairs % for each stair
    for j = 1:2 % 1 = intact, 2 = scrambled
        for k = 1:7 % for each stimulus original = 1:7
            for r = 1:20 % number of reps
                trial_conds_exp = [ trial_conds_exp; i j k ];
            end
        end
    end
end

totaltrials = length(trial_conds_exp);
% add object and orientation column to trial_conds_exp
objcond = [ones(totaltrials/2,1); ones(totaltrials/2,1)+1]; % add bird/human column
%objcond = [ones(totaltrials/2,1)+1; ones(totaltrials/2,1)+1]; % add bird/human column  % temp fix, human/bird only
orientationcol  = [ones(totaltrials/4,1); ones(totaltrials/4,1)+1; ones(totaltrials/4,1); ones(totaltrials/4,1)+1];
%orientationcol  = [ones(totaltrials/4,1)+1; ones(totaltrials/4,1)+1; ones(totaltrials/4,1)+1; ones(totaltrials/4,1)+1]; % temp fix, up/inv only

trial_conds_exp = [objcond orientationcol trial_conds_exp]; % 1 = birds; 2 = humans
% randomly shuffle trial list
ord = randperm(size(trial_conds_exp,1));
trial_conds_exp = trial_conds_exp(ord,:);



% check for Opengl compatibility, abort otherwise:
AssertOpenGL;

% Reseed the random-number generator for each expt.
rand('state',sum(100*clock));

% Make sure keyboard mapping is same across OS
KbName('UnifyKeyNames');

% Define response keys
signalresp = KbName('f');
noiseresp = KbName('j');

% Initiate variable
trial = [];

% Input files
% stimdir_birds = 'C:\Users\Vizcoglab\Desktop\Simen\motion\version_3.1\bird_stims';
% stimdir_humans = 'C:\Users\Vizcoglab\Desktop\Simen\motion\version_3.1\human_stims';
% stimdir_birds_demo = '/home/simen/Dropbox/Research/Victoria/Experiments/Motion/sensivity_experiment/Task/bird_stims/unused_subset';
% stimdir_birds = '/home/simen/Dropbox/Research/Victoria/Experiments/Motion/sensivity_experiment/Task/bird_stims/used_subset';
%stimdir_humans = '/home/simen/Dropbox/Research/Victoria/Experiments/Motion/sensivity_experiment/Task/human_stims';
stimdir_birds_demo = '/Users/u297172/Dropbox/Research/Victoria/Experiments/Motion/sensivity_experiment/Task/bird_stims/unused_subset';
stimdir_birds = '/Users/u297172/Dropbox/Research/Victoria/Experiments/Motion/sensivity_experiment/Task/bird_stims/used_subset';
stimdir_humans = '/Users/u297172/Dropbox/Research/Victoria/Experiments/Motion/sensivity_experiment/Task/human_stims';


%stimdir_birds_demo = '/home/simen/Dropbox/Research/Experiments/Motion/version_3.3/bird_stims/unused_subset';
%stimdir_birds = '/home/simen/Dropbox/Research/Experiments/Motion/version_3.3/bird_stims/used_subset';
%stimdir_humans = '/home/simen/Dropbox/Research/Experiments/Motion/version_3.3/human_stims';
fnames_birds_demo = dir(fullfile(stimdir_birds_demo, '*.mat'));
fnames_birds = dir(fullfile(stimdir_birds,'*.mat'));
fnames_humans = dir(fullfile(stimdir_humans, '*.mat'));

% Variables
rotVect_test = round(150:10:210);         % 7 Human-view points ranging from 150-210  
rotVect_prac = rotVect_test+10;
subjgroup_label = {'pilot', 'control', 'expert'};
objtype_label = {'bird', 'human'};
orientation_label = {'upright', 'inverted'};
intactScramble_label = {'intact', 'scrambled'};


% ------------------------
% set dot field parameters
% ------------------------
%mon_width = 47; % Horizonal dimension of display screen
mon_width = 31; % Horizonal dimension of display screen (Lenovo laptop)
v_dist = 70;    % Viewing distance
dot_w = 0.1;   % Dot size (deg)
fix_r = 0.15;    % Radius of fixation point (deg)
waitframes = 1; % Number of video refresh intervals to show each image before
        
%----------------------
% file/subject handling
%----------------------

% Poll subject information and create output file
repeat=1;
while (repeat)
    prompt= {'Participant number (format: 001)', 'Age', 'Gender', 'Group (1==Pilot; 2==Control; 3==Expert)', 'Block (1 or 2'};
    defaultAnswer={'','','','', ''};
    options.Resize = 'on';
    answer=inputdlg(prompt,'Participant Information',1, defaultAnswer, options);
    [subNo, subjAge, subjGender, subjGroup, blockOrder]=deal(answer{:});
    if isempty(str2num(subNo)) || ~isreal(str2num(subNo)) || isempty(str2num(subjAge)) || ~isreal(str2num(subjAge)) || isempty(num2str(subjGender) || isempty(str2num(subjGroup)) || ~isreal(str2num(subjGroup) || isempty(str2num(blockOrder)) || ~isreal(str2num(blockOrder))))
        h=errordlg('Please fill in all information correctly','Input Error');
        repeat=1;
        uiwait(h);
    else
        if str2num(subjGroup) == 1
            OutputFile=['data/subj_' subNo '_pilot.txt'];	% set data file directory and name of the file
        elseif str2num(subjGroup) == 2
            OutputFile=['data/subj_' subNo '_control.txt'];
        elseif str2num(subjGroup) == 3
            OutputFile=['data/subj_' subNo '_expert.txt'];
        end
        if exist(OutputFile,'file')~=0
            button=questdlg(['Overwrite subj' subNo '.txt?']);
            if strcmp(button,'Yes'); repeat=0; end
        else
            repeat=0;
        end
        % create .txt file to store data
        OutputTextFile = fopen(OutputFile, 'wt');
    end
end

%-------------
% Experiment
%-------------
try
    
    % Retrieve screen number of display screen (max screen number = external
    % screen
    screens = Screen('Screens');
    screenNumber = max(screens);
    
    % Hide the mouse cursos
    HideCursor;
   
    % Open dobbel buffered screen on display screen
    %screenRes = [ 0 0 1024 768 ];
    %Screen('Resolution', screenNumber, screenRes(1), screenRes(2));
    [w, wRect] = Screen('OpenWindow', screenNumber, 0); % wRect(3)=horizontal dim; wRect(4)=vertical dim;
    %[w, wRect] = Screen('OpenWindow', screenNumber, 0, [0 0 720 450]);
    %[w, wRect] = Screen('OpenWindow', screenNumber, 0, [0 0 1024 768]);

        
    % define size of rectangle used during debug. 
    cRect = CenterRect([0 0 500 500], wRect);
    
    % Enable alpha blending with proper blend-function. We need it
    % for drawing of smoothed points:
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    [centerdim(1), centerdim(2)] = RectCenter(wRect);
    center = [0, 0];
    fps = Screen('FrameRate', w);  % frames per second
    ifi = Screen('GetFlipInterval', w);
    if fps==0
        fps=1/ifi;
    end;
    white = WhiteIndex(w);
    
    % Load GetSecs, KbCheck, and WaitSecs to avoid delay during critical
    % times of experiment
    GetSecs;
    KbCheck;
    WaitSecs;
    
    % Priority for script execution = realtime priority
    Priority(MaxPriority(w));
    
    % Run through practice and test phase
    if DEBUG
        startphase = 5; % skip practice
    else
        startphase = 1;
    end
    
    % Select what object category to present in which block
    if str2num(blockOrder) == 1 
        block2ind = 1; 
        block3ind = 2; 
        block2_label = 'bird';
        block3_label = 'human';
        demo1_label = 'Bird';
        demo2_label = 'Human';
    else                
        block2ind = 2;
        block3ind = 1;
        block2_label = 'human';
        block3_label = 'bird';
        demo1_label = 'Human';
        demo2_label = 'Bird';
    end
    
    % Initial instruction screen
    str1 = 'In a given trial you will see dynamic point-light displays (white dots) representing either a flying bird or a human walker.\n\n';
    str2 = 'However, the white-dots are sometimes spatially scrambled (re-arranged) \n so that the GLOBAL pattern of movement is disrupted (LOCAL movement of dots are still present, but spatial configuration is disrupted).\n\n';
    str3 = 'Your task is to discriminate COHERENT (global) from SCRAMBLED (local only) motion of birds and humans. \n\n';
    str4 = 'For example, if it is a bird trial, and you notice that the white dots make up a bird in flight, then PRESS _ f _ for COHERENT motion.\n';
    str5 = 'In contrast, if the dynamic dots DO NOT make up a bird in flight (only scrambled bird motion), then PRESS _ j _ for SCRAMBLED motion. \n\n';
    str6 = 'Finally, to make the task a bit more tricky, the COHEREN and SCRAMBLED displays will be embedded in additional SCRAMBLED dots.\n';
    str7 = 'For example, if COHERENT bird motion is displayed, it will be embedded in additional scrambled dots\n\n';
    str8 = 'Note that the COHERENT and SCRAMBLED motion will always be centered in the middle of the screen\n';
    str9 = 'Once the point-light display is terminated, please respond as accurately and quickly as you can\n\n';
    str10 = 'Click _ mouse _ to see a demo of the point-lights.';
    message = [str1, str2, str3, str4, str5, str6, str7, str8, str9, str10];

    Screen('TextSize', w, 15);
    DrawFormattedText(w, message, 'center', 'center', WhiteIndex(w));
    Screen('Flip', w);
    GetClicks(w);
    %Screen('Flip', w);
    
    startphase = 1;
    %startphase = 4;
    for phase = startphase:8
               
        % Set parameters depending on phase
        if phase == 1
            phaselabel1 = 'Demo Intact';
            phaselabel2 = demo1_label;
            idx = find((trial_conds_prac(:,1)==block2ind & trial_conds_prac(:,2)==1 & trial_conds_prac(:,4)==1 & trial_conds_prac(:,5)==1) | (trial_conds_prac(:,1)==block2ind & trial_conds_prac(:,2)==2 & trial_conds_prac(:,4)==1 & trial_conds_prac(:,5)==1));
            trial_conds = trial_conds_prac(idx,:);
        elseif phase == 2
            phaselabel1 = 'Demo Scrambled';
            phaselabel2 = demo1_label;
            idx = find((trial_conds_prac(:,1)==block2ind & trial_conds_prac(:,2)==1 & trial_conds_prac(:,4)==2 & trial_conds_prac(:,5)==1) | (trial_conds_prac(:,1)==block2ind & trial_conds_prac(:,2)==2 & trial_conds_prac(:,4)==2 & trial_conds_prac(:,5)==1));
            trial_conds = trial_conds_prac(idx,:);
        elseif phase == 3
            phaselabel1 = 'Practice';
            phaselabel2 = demo1_label;
            idx = find(trial_conds_prac(:,1)==block2ind);
            trial_conds = trial_conds_prac(idx,:);
        elseif phase == 4
            phaselabel1 = 'Test';
            phaselabel2 = block2_label;
            idx = find(trial_conds_exp(:,1)==block2ind);   % index rows by block2ind
            trial_conds = trial_conds_exp(idx,:);
        elseif phase == 5 
            phaselabel1 = 'Demo Intact';
            phaselabel2 = demo2_label;
            idx = find((trial_conds_prac(:,1)==block3ind & trial_conds_prac(:,2)==1 & trial_conds_prac(:,4)==1 & trial_conds_prac(:,5)==1) | (trial_conds_prac(:,1)==block3ind & trial_conds_prac(:,2)==2 & trial_conds_prac(:,4)==1 & trial_conds_prac(:,5)==1));
            trial_conds = trial_conds_prac(idx,:);
        elseif phase == 6
            phaselabel1 = 'Demo Scrambled';
            phaselabel2 = demo2_label;
            idx = find((trial_conds_prac(:,1)==block3ind & trial_conds_prac(:,2)==1 & trial_conds_prac(:,4)==2 & trial_conds_prac(:,5)==1) | (trial_conds_prac(:,1)==block3ind & trial_conds_prac(:,2)==2 & trial_conds_prac(:,4)==2 & trial_conds_prac(:,5)==1));
            trial_conds = trial_conds_prac(idx,:);
        elseif phase == 7  
            phaselabel1 = 'Practice';
            phaselabel2 = demo2_label;
            idx = find(trial_conds_prac(:,1)==block3ind);
            trial_conds = trial_conds_prac(idx,:);    
        elseif phase == 8  
            phaselabel1 = 'Test'; 
            phaselabel2 = block3_label;
            idx = find(trial_conds_exp(:,1)==block3ind);
            trial_conds = trial_conds_exp(idx,:);
        end
        
        % Human angles to be displayed
        rotVect = rotVect_prac;
        if phase == 4 || phase == 8
            rotVect = rotVect_test;
        end
            
        
        % total trials if not terminated by stairs
        if ~DEBUG
            ntrials = size(trial_conds,1); 
        else
            ntrials = 7; 
        end

        % Create vectors specifying Break-trials
        breakVect = 75:75:ntrials;    
        
        
        Screen('TextSize', w, 20);
        
        % insert instructions about birds versus humans block. 
        str1 = sprintf('%s Session - Block %i\n', phaselabel1, phase);
        str2 = sprintf('\n\n This block contains upright and inverted %s point-lights \n', phaselabel2);
        str3 = sprintf('Press _ %s _ for COHERENT %s motion \n and _ %s _ for SCRAMBLED %s motion \n', KbName(signalresp), phaselabel2,  KbName(noiseresp), phaselabel2);
        %str4 = '\nRemember that the coherent/scrambled motion will be embedded in additional scrambled dots (i.e., noise).\n';
        message = [str1 str2 ' \n\n' str3 '\n Click _ mouse _ to begin'];
        
        % Draw and display introduction screen
        Screen('TextSize', w, 25);
        DrawFormattedText(w, message, 'center', 'center', WhiteIndex(w));
        Screen('Flip', w);
        
        % Wait for mouse click
        GetClicks(w);
        
        % Flip screen back to background color
        vbl = Screen('Flip', w);
        
        % Wait some time before starting trial
        WaitSecs(1.000);
        
        % ---------------------------------------
        % initialize dot positions and velocities
        % ---------------------------------------
        ppd = pi * (wRect(3)-wRect(1)) / atan(mon_width/v_dist/2) / 360;    % pixels per degree
        s = dot_w * ppd;
        fix_cord = [centerdim-fix_r*ppd centerdim+fix_r*ppd];
        colvect = white;
        
        % Do initial flip...
        vbl=Screen('Flip', w);
        
        trial = 1;
        while ~isempty(trial_conds) && trial <= ntrials
            
            % Take a break every 75 trials. 
            if find(breakVect == trial) >= 1
                str1 = 'Please take a break \n\n Click _ mouse _ when ready to continue';
                DrawFormattedText(w, str1, 'center', 'center', WhiteIndex(w));
                Screen('Flip', w);
                % Flip after mouse click 
                GetClicks(w);
                Screen('Flip', w);
                
            end
            
            % clear data containers
            clear dat
            clear mat
            
           % get conditions from row 1, then delete %%% QCV
            objtypeindex = trial_conds(1,1);     % 1=birds,  2=humans
            orientationindex = trial_conds(1,2); % 1=upright, 2=inverted
            stairindex = trial_conds(1,3);
            intactScrambleindex = trial_conds(1,4);
            stimuliindex = trial_conds(1,5);
            trial_conds(1,:) = [];
            
            if phase == 4 || phase == 8
                
                % check for end of stair
                if stair{stairindex}.reversal.count==stair{stairindex}.reversal.stop+1
                    idx = find(trial_conds(:,3)==stairindex);
                    trial_conds(idx,:) = []; % delete remaining stair
                end
                
                % update tTest %%% QCV
                if isempty(stair{stairindex}.intensity.level)
                    tTest = stair{stairindex}.startintensity;
                else
                    tTest = stair{stairindex}.intensity.level(end);
                end
                
            end
            Screen('Flip', w);
            
            % ITI
            WaitSecs(0.500);
            
            % initialize KbCheck and variables to make sure they're
            % procperly initialized/allocted by Matlab
            [KeyIsDown, endrt, KeyCode] = KbCheck;
            %[KeyIsDown, endrt, KeyCode] = PsychHID('KbCheck', keyboardIndices);

            
            % Load stimulus            
                if objtypeindex == 1                    % 1 = birds
                    if phase == 1 || phase == 2 || phase == 3 || phase == 5 || phase == 6 || phase == 7 % demo and practice-stage: select different stimuli than experiment set
                        stimdir_objects = stimdir_birds_demo;
                        fnames_objects = fnames_birds_demo;
                    else
                        stimdir_objects = stimdir_birds;
                        fnames_objects = fnames_birds;
                    end

                    load(fullfile(stimdir_objects,fnames_objects(stimuliindex).name),'dat', 'mov');
                    nframes = size(dat,3);              
                else
                    stimdir_objects = stimdir_humans;    % 2 = humans
                    fnames_objects = fnames_humans;
                    load(fullfile(stimdir_objects,fnames_objects(1).name)); % only one walker: thus index = 1 
                    nframes = size(mat,3);
                end
              
            %-----------------------
            % if bird, process
            %-----------------------
            if exist('dat')
                % remove occluded dots (replace with NaN)
                szx = size(mov,2);
                for i = 1:size(mov,4)
                    idx = find(dat(:,1,i) < 20 & dat(:,2,i) > (szx-20));
                    if ~isempty(idx)
                        dat(idx,:,i) = NaN;
                    end
                end
                span = 5;
                dats = smooth_points(dat,span); % smooth dots (moving average)
                datt = remove_translation(dats,1); % remove translation of smoothed dots
                mat = datt;
                rotflag = 0;
            else
                rotflag = 1;        % rotate or not
                rotdeg = rotVect(stimuliindex);  
                whichaxis = 1;  % y-axis of rotation: upright
                plotflag = 0;   % show dots
                mat = rotxyz(mat,rotdeg,whichaxis,plotflag);
            end
            
            % create noise
            if phase == 1 || phase == 2 || phase == 5 || phase == 6
                ndots = 0;
            elseif phase == 3 || phase == 7
                ndots = 8;
            elseif phase == 4 || phase == 8
                ndots = tTest; % Original
                %ndots = 8; % temp fix, rec screen; ndots used in mov
                %ndots = 0; % temp fix for recording screen

            end
            maxdots = 3000;
            %%%%%%%%%%%%%%%%%%%%%%%%
            % make sure matrix is N x 2 x FRAMES
            %%%%%%%%%%%%%%%%%%%%%%%%
            dim = size(mat);
            mattmp = zeros(dim(1),2,dim(3));
            for n = 1:dim(3)
                mattmp(:,:,n) = mat(:,1:2,n);
            end
            mat = mattmp;
            dim = size(mat);

             % scale size of birds to approximately the same size
             if objtypeindex == 1 
                if stimuliindex == 1
                    scalefactor = 2.3;
                elseif stimuliindex == 2
                    scalefactor = 1.9;
                elseif stimuliindex == 4
                    scalefactor = 1.9;
                elseif stimuliindex == 5
                    scalefactor = 1.3;
                elseif stimuliindex == 6
                    scalefactor = 1.8;
                elseif stimuliindex == 3 || stimuliindex == 7
                    scalefactor = 1.7;
                end
                mat = mat * scalefactor; 
             end
                       
            % Center bird stimulus
            if objtypeindex == 1
                centeringvect = centerdim - mat(6,:,1);
            else
                centeringvect = centerdim - mat(7,:,1);
            end
            centeringvect = repmat(centeringvect,[size(mat,1),1,nframes]);
            mat = mat + centeringvect;
            
            % randomly flip stimuli along horizontal axis
            hor_flip_idx = randi(2);
            if hor_flip_idx == 2
                hor_flip_vect = repmat([wRect(3),0],[size(mat,1),1,nframes]);
                mat_temp = hor_flip_vect - mat;
                mat(:,1,:) = mat_temp(:,1,:);
            end
            
            % if inversion trial, then invert the stimuli.   
            if orientationindex == 2   % 2 == inversion
%                invertvect = repmat(wRect(3:4),[size(mat,1),1,nframes]);
                 invertvect = repmat([0,wRect(4)],[size(mat,1),1,nframes]); 
                 mat_temp = invertvect - mat;
                 mat(:,2,:) = mat_temp(:,2,:);
%                mat = invertvect - mat;
            end
                             
            %%%%%%%%%%%%%%%%%%%%%%%%
            % distribute evenly as possible in an arbitrary 3 x 3 grid
            %%%%%%%%%%%%%%%%%%%%%%%%
            xx = reshape(squeeze((mat(:,1,:))),dim(1)*dim(3),1);
            yy = reshape(squeeze((mat(:,2,:))),dim(1)*dim(3),1);
            xmn = floor(min(xx)); xmx = floor(max(xx)); ymn = floor(min(yy)); ymx = floor(max(yy));
            xs = floor((xmx-xmn)/3); ys = floor((ymx-ymn)/3);
            recti = [];
            for i = 1:3
                for j = 1:3
                    recti = [ recti; xmn+(i-1)*xs+1 xmn+(i*xs) ymn+(j-1)*ys+1 ymn+(j*ys) ];
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % create an arbitrarily large set of dots (maxdots)
            %%%%%%%%%%%%%%%%%%%%%%%%
            scramblevect_signal = [];
            grid_index = [];
            g = 1; % keep track of grid position
            for n = 1:maxdots
                scramblevect_signal = [ scramblevect_signal; randi([recti(g,1),recti(g,2)],[1,1]) randi([recti(g,3),recti(g,4)],[1,1]) ];
                grid_index = [ grid_index; g ];
                g = g + 1;
                if g > 9, g = 1; end
            end
            temp_index = randi(dim(1),maxdots,1);
            temp_noise = mat(temp_index,:,:);
            temp_scramble = [];
            for n = 1:maxdots
                temp_scramble = [ temp_scramble; scramblevect_signal(n,:) ];
            end
            temp_scramble = repmat(temp_scramble,[1,1,dim(3)]);
            tmp = temp_noise - repmat(temp_noise(:,:,1),[1,1,dim(3)]);
            matx = tmp + temp_scramble;
            
            % throw out dots that goes beyond border by jit amount, can change jit
            % amount
            checkdots = ones(maxdots,1);
            jit = 25;
            for i = 1:maxdots
                if any(matx(i,1,:)<(xmn-jit)) | any(matx(i,1,:)>(xmx+jit)) | any(matx(i,2,:)<(ymn-jit)) | any(matx(i,2,:)>(ymx+jit))
                    checkdots(i) = 0;
                end
            end
            idxo = find(checkdots==1);
            
            % for checking final distribution
            check_params = [
                temp_index(idxo(ndots+1:ndots+dim(1))) grid_index(idxo(ndots+1:ndots+dim(1))); ...
                temp_index(idxo(1:ndots)) grid_index(idxo(1:ndots)) ];
            
            mats = matx(idxo(ndots+1:ndots+dim(1)),:,:);    %%% FINAL SCR MATRIX
            matx = matx(idxo(1:ndots),:,:);                 %%% FINAL NOISE MATRIX
            
            if intactScrambleindex == 1
                mat_ind = mat;
                intactScramble_name = 'intact';
            else
                mat_ind = mats;
                intactScramble_name = 'scrambled';
            end
            
            % Draw initial fixation point
            Screen('FillOval', w, uint8(white), fix_cord); % fixation point until next flip
            Screen('TextSize', w, 24);                     % set txt size for msg screen
            Screen('Flip', w);
            WaitSecs(0.500);
            Screen('Flip', w);
            
            % Birds: randomly select a starting frame from the first 4 seconds (4sec*33frames/sec), and
            % play for 75 frames (2.5 seconds)
            if objtypeindex == 1 && stimuliindex == 6
                nframes_birds = 87+30;
                startframe = randi([30,nframes_birds-75]); % skip beginning frames of this move due to weird dots. 
                endframe = startframe + 75;
            elseif objtypeindex == 1 && (stimuliindex == 1 || stimuliindex == 2 || stimuliindex == 3 || stimuliindex == 4 || stimuliindex == 5 || stimuliindex == 7)
                nframes_birds = 87;         % 87 frames = max(size(duck2). Select 87 first frames of all the movies, then select a random subset of 75 frames (later in the script).  
                startframe = randi(nframes_birds-75);
               endframe = startframe+75;
            elseif objtypeindex == 2
                startframe = 1;
                endframe = nframes;
            end
             
            stim_onset = GetSecs; % used to log time of stimuli duration
            for frame = startframe:endframe
                
                if ~DEBUG
                    % combine stimuli and noise
                    xymatrix = vertcat(mat_ind(:,:,frame), matx(:,:,frame));
                    Screen('DrawDots', w, xymatrix', s, colvect, center, 1);  % change 1 to 0 or 4 to draw square dots
                else
                    Screen('DrawDots', w, mat_ind(:,:,frame)', s, colvect, center, 1);  % change 1 to 0 or 4 to draw square dots
                    Screen('DrawDots', w, matx(:,:,frame)', s, [255 0 0], center, 1);  % change 1 to 0 or 4 to draw square dots
                    DrawFormattedText(w, sprintf('%s, ndots = %i, %s',fnames_objects(1).name,ndots,intactScramble_label{intactScrambleindex}), 'center', 700, WhiteIndex(w));
                    Screen('FrameRect', w, [255 255 255], cRect);
                end
                Screen('DrawingFinished', w); %Tell PTB that no further drawing commands will follow before Screen('Flip')
                                
                %vbl = Screen('Flip', w, vbl + 60/25*((waitframes - 0.5)*ifi));  % Presents at 25Hz?
                %vbl = Screen('Flip', w, vbl + 2*((waitframes - 0.5)*ifi)); % Presents at 30Hz (60/30=2) %ORIGINAL
                vbl = Screen('Flip', w, vbl + 4*((waitframes - 0.5)*ifi)); % Presents at 30Hz (120/30=4) % for lenovo (refresh rate=120)

                   
            end
            
             % record stimuli duration for birds.
             stim_offset = GetSecs;
             stim_dur = stim_offset-stim_onset;

            
            WaitSecs(0.100);
            
            % Create response screen.
            DrawFormattedText(w, str3, 'center', 'center', WhiteIndex(w));
            [VBLTimestamp startrt] = Screen('Flip', w);
            
            repeat = 1;
            while repeat
                % poll for a resp
                % during test phase, subjects can response
                % before stimulus terminates
                
                if ( KeyCode(signalresp)==1 || KeyCode(noiseresp)==1 )
                    repeat = 0;
                    break;
                end
                [KeyIsDown, endrt, KeyCode]=KbCheck;
                %[KeyIsDown, endrt, KeyCode] = PsychHID('KbCheck', keyboardIndices);

                
                % Wait 1 ms before checking the keyboard again to prevent
                % overload of the machine at elevated Priority():
                WaitSecs(0.001);
            end
            
            % Compute response time
            rt = round(1000*(endrt-startrt));
            
            % Compute acc
            if ((intactScrambleindex == 1 && KeyCode(signalresp) == 1) || (intactScrambleindex == 2 && KeyCode(noiseresp) == 1))
                acc = 1;
                acc_text = 'Correct';
            else
                acc = 0;
                acc_text = 'Incorrect';
            end
            
            % Reverse acc to create response variable used for QUEST.
            if acc == 0
                response = 1;
            elseif acc  == 1
                response = 0;
            end
            
            % Record key pressed
            keyresp = KbName(find(KeyCode,1,'first'));
            %keyresp = KbName(KeyCode);
            
            % Log correct response
            if intactScrambleindex == 1
                corrResp = KbName(signalresp);
            elseif intactScrambleindex == 2
                corrResp = KbName(noiseresp);
            end
            
            % log stimuli name
            if objtypeindex == 1
                [stimpath, stimname, ext] = fileparts(fnames_objects(stimuliindex).name);
            else
                stimname = strcat('humanWalk_', num2str(rotVect(stimuliindex)), 'deg');
            end
            
            % update stair %%% QCV
            if phase == 4 || phase == 8
                stair{stairindex} = STAIR_update(stair{stairindex},response);
            end
            % Provide feedback during practice.
            if phase == 1 || phase == 2 || phase == 3 || phase == 5 || phase == 6 || phase == 7
                tTest = ndots;
                str_temp = sprintf('\n\n Practice trial %i out of %i', trial, ntrials);
                text = [acc_text str_temp];
                DrawFormattedText(w, text, 'center', 'center', WhiteIndex(w));
                Screen('Flip', w);
                WaitSecs(2.0);
            end
            
            % Log results
            fprintf(OutputTextFile, '\n %s\t %s\t %s\t %s\t %s\t %s\t %i\t %i\t %i\t %s\t %s\t %i\t %s\t %i\t %s\t %i\t %i\t %s\t %s\t %i\t %i\t %i\t', ...
                subNo, ...
                subjAge, ...
                subjGender, ...
                subjGroup, ...
                subjgroup_label{str2num(subjGroup)}, ...
                blockOrder, ...
                phase, ...
                trial, ...
                objtypeindex, ...
                objtype_label{objtypeindex}, ...
                stimname, ...
                orientationindex, ...
                orientation_label{orientationindex}, ...
                intactScrambleindex, ...
                intactScramble_label{intactScrambleindex}, ... 
                stairindex, ...
                stim_dur, ...
                keyresp, ...
                corrResp, ...
                rt , ...
                acc , ...
                tTest);
            
            % update trial
            trial = trial + 1;
            
        end
        
    end
    
    % save mat file
    output_name = sprintf('stair_subj%s_%s', subNo, subjgroup_label{str2num(subjGroup)});
    save(output_name, 'stair');
    
    % Ending screen
    finishtxt = 'Thank you for your participation \n\n Click _ mouse _ to exit. ';
    DrawFormattedText(w, finishtxt, 'center', 'center', WhiteIndex(w));
    Screen('Flip',w);
    
    GetClicks(w);
    
    % Cleanup at end of experiment
    sca; % close window
    ShowCursor; % show mouse cursor
    fclose('all');  % close result file
    Priority(0); % switch back to normal priority
catch
    sca;
    ShowCursor;
    fclose('all');
    Priority(0);
    
    % output last error msg
    psychrethrow(psychlasterror)
end % try ... catch


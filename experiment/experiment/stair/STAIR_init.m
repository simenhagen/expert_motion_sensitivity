function s = STAIR_init(startintensity,initsteps,finalstep,rules)

%
% startintensity: integer value for start intensity e.g., 10
% initialsteps: matrix with:
%        row 1: desired initial step sizes 
%        row 2: number of changes for each step size in row 1
%        e.g., [ 5 2; 2 2 ];
% finalstep: integer value for final step size, e.g., 1
% rules: row array with the up rule, down rule, and stop rule e.g., [ 1 3 8 ]
%
% NOTE: no checks so make sure input is correct!
%
% 08/12/14, qcv
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set initial intensity level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s.startintensity = startintensity;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set step size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% can have some initial large steps if desired
s.step.initialsizes = initsteps(1,:);   % step size
s.step.initialns = initsteps(2,:);      % how many times at each size
s.step.finalsize = finalstep;           % set final step size, also start tracking reversals from here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set rules:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 up/1 down = 50.0%
% 1 up/2 down = 70.7%
% 1 up/3 down = 79.4%
% 1 up/4 down = 84.1%
s.rule.up.val = rules(1);       % response = 0
s.rule.up.c = 0;                % count the ups
s.rule.down.val = rules(2);     % response = 1    
s.rule.down.c = 0;              % count the downs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set reversals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s.reversal.seq = [];            % store rev sequence (1 = rev, 0 = no rev)
s.reversal.count = 1;
s.reversal.stop = rules(3);     % stopping rule, counting starts from final step size
s.intensity.changecount = 1;    % count changes (irrespective of direction)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% store intensity and responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s.intensity.level = [];
s.intensity.dir = [];           % -1 = decrease, 1 = increase, 0 = no change
s.response = [];

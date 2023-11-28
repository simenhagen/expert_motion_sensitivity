function s = STAIR_update(s,response)

% update response
s.response = [ s.response response ];

% get current intensity
if isempty(s.intensity.level)
    intensity = s.startintensity;
else
    intensity = s.intensity.level(end);
end

% generate step sequence for initial changes
steps = [];
for ii = 1:size(s.step.initialsizes,2)
    steps = [ steps s.step.initialsizes(ii)*ones(1,s.step.initialns(ii)) ];
end

% get current step size
if s.intensity.changecount <= size(steps,2)
    step = steps(s.intensity.changecount);
else
    step = s.step.finalsize;
end

reversed = 0;

% get current rule (up or down)
if response
    %curr_rule = s.rule.down;
    dir = -1;
    s.rule.down.c = s.rule.down.c + 1;
    if s.rule.down.c == s.rule.down.val
        intensity = intensity + dir*step;
%         if intensity < 10 % HACK to prevent 0 noisedots
%             intensity = s.startintensity;
%         end
        if intensity < 1 % HACK to prevent 0 noisedots
            intensity = 0;
        end
        s.intensity.dir = [ s.intensity.dir dir ];  % store decrease
        s.rule.down.c = 0;                          % reset counter
        s.rule.up.c = 0;                            % need to also reset other counter
        s.intensity.changecount = s.intensity.changecount + 1;
        % check for reversal
        tmp = fliplr(s.intensity.dir);
        if size(tmp,2) > s.rule.down.val
            if tmp(s.rule.down.val+1)~=dir
                reversed = 1;
            end
        end
    else
        s.intensity.dir = [ s.intensity.dir 0 ];    % store no change
    end
else
    %curr_rule = s.rule.up;
    dir = 1;
    s.rule.up.c = s.rule.up.c + 1;
    if s.rule.up.c == s.rule.up.val
        intensity = intensity + dir*step;
        s.intensity.dir = [ s.intensity.dir dir ];  % store decrease
        s.rule.up.c = 0;                            % reset counter
        s.rule.down.c = 0;                          % need to also reset other counter
        s.intensity.changecount = s.intensity.changecount + 1;
        % check for reversal
        tmp = fliplr(s.intensity.dir);
        if size(tmp,2) > s.rule.up.val
            if tmp(s.rule.up.val+1)~=dir
                reversed = 1;
            end
        end
    else
        s.intensity.dir = [ s.intensity.dir 0 ];    % store no change
    end
end

%update intensity
s.intensity.level = [ s.intensity.level intensity ];

% update reversal
if step == s.step.finalsize
    s.reversal.count = s.reversal.count + reversed;
end
s.reversal.seq = [ s.reversal.seq reversed ];


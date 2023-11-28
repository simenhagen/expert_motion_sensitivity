function STAIR_plot(stairmat)

close all

% load mat file
% fname = sprintf('ADAPT_stair_%03d.mat',subno);
% load(fname);

st = stairmat;

nstairs = size(st,2);

for idx = 1:nstairs
    
    intensity = [ st{idx}.startintensity st{idx}.intensity.level(1:end-1) ];
    response = st{idx}.response;
    
    if st{idx}.rule.up.val > st{idx}.rule.down.val
        response = 1 - response;
    end
    
    % compute overall ac
    ac = mean(response)*100;
    
    figure(idx);
    hold on
    for rr = 1:length(response)
        if response(rr)==1
            plot(rr,intensity(rr),'ko');
        else
            plot(rr,intensity(rr),'ro');
        end
    end
    xlabel('trial');
    ylabel('stimulus intensity');
    axis([0.5 length(response)+0.5 min(intensity)-0.5 max(intensity)+0.5 ]);
    plot(intensity);
    text(length(response)-4,min(intensity)+5,sprintf('Accuracy\n%.1f%%',ac));
    
end

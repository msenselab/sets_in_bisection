
%% auditory presentation
function presentAudio(a, tone_s, dur)
    fixdur = 0.3+0.3*rand;%random 0.3-0.6s fixation duration
    tone_start_idx = round((fixdur)*a.freq);  %28800- 57600
    tone_end_idx = dur*a.freq + tone_start_idx-1;
    mix_range = tone_start_idx:tone_end_idx;
    tone_mix(:,mix_range)=tone_s(:,mix_range);
    a.prepare(tone_mix);
    clearvars  tone_mix;
    a.present;
    WaitSecs(fixdur + dur); % keep same for all trials
    a.stop; %stop when the sound is finished
    %WaitSecs(0.3);
end


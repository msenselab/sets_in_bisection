  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name:          bisection_exp2.m
%%
%% Author:         Fiona Zhu
%% Description:   Main function for running the bisection study to
%% investigate the effect of mean of the distribution while keeping the spacing the same.
%%
%% note: 
%%      seven intervals (400, 600, 800, 1000, 1200, 1400, and 1600 ms)
%%      two types of sampled distributions will be tested
%%          Condition A has a descending frequency distribution of durations and
%%          condition B has an ascending frequency distribution of durations.
%% Date:          17-11-2018
%% Modified by Xiuna Zhu (Fiona)
%% Contact:   fiona.zhu1230@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bisection_exp2
% -------------------------------------
% This is an example code for bisection task to study frequency distribution
% Participants have to identify the presented visual stimulus was close to the short or to the long
% note:
%--------------------------------------

try
    num_dur = 7;
    num_dur_all = 1/2*num_dur*(num_dur+1);
    blockTrials =num_dur_all*2; % number of trials per block
  
    % block factor : frequency distribution A or B
    % A: descending frequency distribution
    % B: ascending frequency distribution
    exp = CExp(2,[1 num_dur_all],'blockRepetition', 6, 'blockFactors', 2);%1.modality(only audio) 2. duration  3. frequency distribution A or B
    
    %parameters
    frequency = 1000;
    durations = [400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0]/1000;
    %frequency distribution
    durList = 1: num_dur_all;
    durList_A = 1: num_dur_all;
    durList_B = 1: num_dur_all;
    Freq = 1: num_dur;
    for i = 1: num_dur
        j = 1/2*i*(i-1)+1;
        durList_B(j:1/2*i*(i+1)) = durations(i);  %B ascending frequency distribution
        durList_A(j:1/2*i*(i+1)) = durations(1+num_dur-i); % A: descending frequency distribution
        %definition of Fi = i/(n+1); 
        Freq(i) = 2*i/(num_dur*(num_dur+1));
    end

    % enquire subject information
    exp.subInfo('session: AF(1) or DF(2)?)', '1');
    if(strcmp(exp.sName, 'cancelled'))
        disp('This experiment has been cancelled!')
        return;
    end
    
    %set the order
    if(exp.sPara == 1)
        exp.seq = sortrows(exp.seq, 3);  %A: ascending frequency condition firstly
    else
        exp.seq = sortrows(exp.seq, -3); %B: descending frequency condition firstly
    end
    
    % new two practice blocks
    prac_seq = exp.genTrials(2,[1 num_dur_all], 2);
    
    % add 2 practice blocks' trials at beginning of part 1
    exp.seq = [prac_seq(1:blockTrials*2,:); exp.seq];
    exp.maxTrls = exp.maxTrls + blockTrials*2;
    
    
    % input device
    kb=CInput('k',[0 1],{'leftArrow','rightArrow'});
    
    % audio device
    %devices = PsychPortAudio('GetDevices' );
    %Option1: Using the computer's default audio device
    a = CAudio;
    %Option2: Using the Motu's audio device
    %a = CAudiom('hardware', 'Motu', 'channels',2); % 2 channels
       
    
    % visual display
    % v = CDisplay('bgColor',48,'fontSize',18, 'skipSync',1); % for test
    v = CDisplay('bgColor',48,'fontSize',18,'monitorSize',22, 'fullWindow', 1);
    
    HideCursor;
    
    infoText = init_text;
    
    %        ltone = a.genTone(500, durList(end)); % long tone
    noisedur = 2*durations(end);
    noise = a.genPinkNoise(frequency/2,noisedur,-50)/noisedur; % 60 dBA
    tone =  a.genTone(frequency*4, noisedur)*1.6; % 60 dBA paras: obj, freq, duration
    tone_s = tone+noise; %mixing noise and tone in stumuli
    tone_mix = noise;  % 2*144000
    warning = a.genTone(4000,0.02); % a warning signal
    
    v.dispText(infoText.instruction);
    kb.wait;
    for iTrl=1:exp.maxTrls
        if (mod(iTrl, blockTrials) == 1) && (iTrl == 0.5 * exp.maxTrls + 1) % session info
            v.dispText(infoText.startSession);
            kb.wait;
        elseif mod(iTrl, blockTrials) == 1  % block info
            txtBlockInfo = ['Block ' num2str(floor(iTrl/blockTrials)+1) '\n' infoText.startBlock];
            v.dispText(txtBlockInfo);
            kb.wait;
        end
        
        cond = exp.getCondition;
        if(cond(3) == 1)
            durList = durList_A;
        else
            durList = durList_B;
        end
        
        % fixtion with a warning tone
        v.dispFixation(1);
        a.prepare(warning);
        a.present;
        WaitSecs(0.5);
        a.stop;
        
        %present audio
        %presentAudio(a, tone_s, durList(cond(2)));

        fixdur = 0.3+0.3*rand;%random 0.3-0.6s fixation duration
        tone_start_idx = round((fixdur)*a.freq);  %28800- 57600
        tone_end_idx = durList(cond(2))*a.freq + tone_start_idx-1;
        mix_range = tone_start_idx:tone_end_idx;
        tone_mix(:,mix_range)=tone_s(:,mix_range);
        a.prepare(tone_mix);
        a.present;
        WaitSecs(0.3 + durList(cond(2))); % keep same for all trials
        a.stop; %stop when the sound is finished
        %WaitSecs(0.3);
        v.flip(1);
        %acquire response
        v.dispText(infoText.question);
        key = kb.response;
        if iTrl <= blockTrials*2
            % feedback
            feedbacktext = infoText.short;
            if durList(cond(2)) == durations((num_dur+1)/2)
                if(rand() >= 0.5)
                    feedbacktext = infoText.long;
                end
            elseif durList(cond(2)) > durations((num_dur+1)/2)
                feedbacktext = infoText.long;
            end
            v.dispText(feedbacktext);
            WaitSecs(1);
        end
        
        % store duration and responses
        exp.setResp([durList(cond(2)) key]);
        v.flip(1); %clear screen
        % inter trial interval (ITI)
        WaitSecs(0.9+0.2*rand);
        
        if  mod(iTrl,blockTrials) == 0 && iTrl < exp.maxTrls
            v.dispText(infoText.blockInfo);
            kb.wait;
        end
        if kb.wantStop
            v.dispText(infoText.stopInfo);
            disp(infoText.stopInfo);
            break;
        end
    end
    
    
    %closing the experiment
    exp.saveData;   %save data
    
    v.dispText(infoText.thankyou);
    kb.wait;
    v.close;
    a.close;
    
catch ME
    disp(ME.message);
    disp(ME.stack);
    for i=1:length(ME.stack)
        disp(ME.stack(i).name);
        disp(ME.stack(i).line);
    end
    v.close;
    a.close;
end
end

function infoText = init_text
% specify experimental text
infoText.instruction = ['Duration Bisection Experiment \n Please follow the instruction the experimenter gave.\n',...
    ' In this experiment, you will hear a tone or see a square. \nYour task is to judge if the stimuli presented is a short or long stimuli.\n', ...
    ' You will learn what long or short is after few trials.'];
infoText.practice = 'This is the practice block \n Press a key to begin';
infoText.formalInfo = 'The formal experiment will start soon, \n please press a key to start\n';
infoText.blockInfo = ['Please take a rest, and when you are ready, please press a key to start a new block. \n ', ...
    'Your task is to judge if the given stimulus is short (left arrow key) or long (right arrow key)'];
%infoText.question = '? \n\n\n\n\n Short (<-) or Long (->)?';
infoText.question = '?';
infoText.stopInfo = 'Stop in runing experiment...';
infoText.startBlock = 'Please press any key to start this block';
infoText.thankyou = 'The experiment is finished! \nThank you very much!';
infoText.short = 'The given stimulus is short!';
infoText.long = 'The given stimulus is long!';
infoText.startSession = ['The first part was finished. \n '...
    'When you are ready, please press a key to start the second part. \n ', ...
    'Your task is to judge if the given stimulus is short (left arrow key) or long (right arrow key)'];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name:          bisection_exp1.m
%%
%% Author:         Fiona Zhu
%% Description:   Main function for running the bisection study to
%% investigate the spacing effect.
%%
%% note:
%%      positive skew condition (400, 504, 636, 800, 1008, 1270, and 1600 ms)
%%      negative skew condition (400, 730, 992, 1200, 1366, 1496, and 1600 ms)
%%      two types of condition will be tested
%% Date:          17-11-2018
%% Modified by Xiuna Zhu (Fiona)
%% Contact:   fiona.zhu1230@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bisection_exp1
% -------------------------------------
% This is an example code for bisection task to study frequency distribution
% Participants have to identify the presented visual stimulus was close to the short or to the long
%--------------------------------------

try
    num_dur = 7;
    pracBlockNum = 2;
    num_dur_all = 1/2*num_dur*(num_dur+1);
    blockTrials =num_dur_all*2; % number of trials per block (56 trials)
    
    % block factor : frequency distribution A or B
    % A: positive skew condition
    % B: negative skew condition
    exp = CExp(2,[num_dur_all],'blockRepetition', 6, 'blockFactors', 2);  %672 trials
    %1.modality(only audio) 2. duration  3. frequency distribution A or B
    
    %parameters
    frequency = 1000;
    dur =  [400.0, 504.0, 636.0, 800.0, 1008.0, 1270.0, 1600.0;
        400.0, 730.0, 992.0, 1200.0, 1366.0, 1496.0, 1600.0]/1000;
    
    
    weigth = [1   1     1     1   1   1   1];
    currmean1 = sum(dur(1,:).*weigth/7); %code to check if mean == 0.8883
    currmean2 = sum(dur(2,:).*weigth/7); %code to check if mean == 1.1120
    x_mean1_2 = (dur(1,:)-currmean1).^2;
    x_mean2_2 = (dur(2,:)-currmean2).^2;
    sigma1 = sqrt(sum(x_mean1_2.*weigth/7)); %sigma1 =  0.4007
    sigma2 = sqrt(sum(x_mean2_2.*weigth/7)); %sigma2 =  0.4009
    
    
    
    meanDur = 1;
    durList = repmat(dur,1,4);
  
    
    % mean of first set mean(dur(1,:)) = 0.8883, and mean of second set mean(dur(2,:)) = 1.1120
    
    % enquire subject information
    exp.subInfo('Positive: 1(Y) or 2(N)?)', '1');
    if(strcmp(exp.sName, 'cancelled'))
        disp('This experiment has been cancelled!')
        return;
    end
    
    %set the order
    if(exp.sPara == 1)
        exp.seq = sortrows(exp.seq, 2);  %A: positive skew condition firstly
    else
        exp.seq = sortrows(exp.seq, -2); %B: negative skew condition firstly
    end
    
    % new two practice blocks
    if pracBlockNum > 0
        prac_seq = exp.genTrials(2,[num_dur_all], pracBlockNum);
        %add 2 practice blocks' trials at beginning of part 1
        exp.seq = [prac_seq(1:blockTrials*pracBlockNum,:); exp.seq];
        exp.maxTrls = exp.maxTrls + blockTrials*pracBlockNum;
    end
    
    
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
    
    
    infoText = init_text;
    noisedur = 2*dur(end);
    noise = a.genPinkNoise(frequency/2,noisedur,-50)/noisedur; % 60 dBA
    tone =  a.genTone(frequency*4, noisedur)*1.6; % 60 dBA paras: obj, freq, duration
    tone_s = tone+noise; %mixing noise and tone in stumuli
    tone_mix = noise;  % 2*144000
    warning = a.genTone(4000,0.02); % a warning signal
    
    % display instructions and wait for keypress
    pracTex = '\n There are ';
    if pracBlockNum == 1
        pracTex =[pracTex  'one practice block,'];
    elseif pracBlockNum > 1
        pracTex = [pracTex int2str(pracBlockNum) ' practice blocks, and '];
    end
    
    blockText = [pracTex int2str(exp.maxTrls/blockTrials) ' blocks in total.\n\n'];
    v.dispText([infoText.instruction, blockText, infoText.startBlock]);
    
    kb.wait;
    for iTrl=1:exp.maxTrls
        if pracBlockNum >= 1 && mod(iTrl, blockTrials) == 1 && floor(iTrl/blockTrials) < pracBlockNum
            v.dispText(['Practice Block ' num2str(floor(iTrl/blockTrials)+1) ' \n Please press any key to start']);
            kb.wait;
        elseif (mod(iTrl, blockTrials) == 1) && (floor(iTrl/blockTrials) == (exp.maxTrls/(2*blockTrials) + pracBlockNum))
            % show session info at block 16
            v.dispText(infoText.startSession);
            kb.wait;
        elseif mod((iTrl- pracBlockNum*blockTrials-1), blockTrials) == 0
            v.dispText(['Block ' num2str(floor((iTrl-pracBlockNum*blockTrials)/blockTrials)+1) '\nPlease press any key to start this block']);
            kb.wait;
        end
        
        
        cond = exp.getCondition;
        targetDur = durList(cond(2), cond(1));
        
        % fixtion with a warning tone
        v.dispFixation(1);
        a.prepare(warning);
        a.present;
        WaitSecs(0.5);
        a.stop;
        
        %present audio
        
        fixdur = 0.3+0.3*rand;%random 0.3-0.6s fixation duration
        tone_start_idx = round((fixdur)*a.freq);  %28800- 57600
        tone_end_idx = targetDur*a.freq + tone_start_idx-1;
        mix_range = tone_start_idx:tone_end_idx;
        tone_mix(:,mix_range)=tone_s(:,mix_range);
        a.prepare(tone_mix);
        a.present;
        WaitSecs(0.3 + targetDur); % keep same for all trials
        a.stop; %stop when the sound is finished
        %WaitSecs(0.3);
        v.flip(1);
        %acquire response
        v.dispText(infoText.question);
        key = kb.response;
        if iTrl <= blockTrials*2
            % feedback
            feedbacktext = infoText.short;
            if targetDur  >  meanDur
                feedbacktext = infoText.long;
            end
            v.dispText(feedbacktext);
            WaitSecs(1);
        end
        
        % store duration and responses
        exp.setResp([targetDur key]);
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
    ' In this experiment, you will hear a tone then a noise. \nYour task is to judge if the noise presented is a short or long stimuli.\n', ...
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
infoText.startSession = ['The first part was finished. Please take a longer break. \n '...
    'When you are ready, please press a key to start the second part. \n ', ...
    'Your task is to judge if the given stimulus is short (left arrow key) or long (right arrow key)'];

end

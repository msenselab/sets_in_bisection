%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name:          bisection_S3.m 
%%
%% Author:         Fiona Zhu
%% Description:   Main function for running the bisection study 2 to
%% investigate the influence of ensemble variance on temporal bisection 
%%
%% note: 
%% two types of sampled distributions will be tested: 
%%       a normal, and a U-shaped distribution 
%% 8 intervals (400,  550, 700, 850, 1000, 1150, 1300, and 1450 ms)
%% Two practice blocks : uniform distribution 
%% Author: Xiuna Zhu (Fiona)
%% Date: 22-1-2018
%% Contact: fiona.zhu1230@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bisection_S3
% -------------------------------------
% This is the code for bisection task to study2 with a U-shaped
% distribution and normal distribution
% weight = [1/8      1/8    1/8      1/8     1/8     1/8     1/8      1/8;  %uniformal
%          30/72    2/72    2/72     2/72    2/72    2/72    2/72    30/72;  %U-shaped
%          2/72    2/72     2/72    30/72    30/72   2/72    2/72    2/72 ]; %normal
%Participants have to identify the presented visual stimulus was close to the short or to the long
% note:
%--------------------------------------

try
    close all;  clc;
    num_dur = 8;
    pracBlockNum = 2;
    num_dur_all = 1/2*num_dur*(num_dur+1); %36
    blockTrials =num_dur_all*2; % number of trials per block 72
    exp = CExp(2,[1 blockTrials], 'blockRepetition', 2, 'blockFactors', 2);
    %1.    
    %2. duration 
    %3. distribution: 1 = uniform; 2 = U-shaped; 3 = normal. 
    
    
    idx = find(exp.seq(:,3)==1);
    exp.seq(idx,3) = 3;
    
    frequency = 1000;
    durations = [400.0, 550.0, 700.0, 850.0, 1000.0, 1150.0, 1300.0, 1450.0]/1000;
    durList = ones(3, blockTrials);

    % parameters for uniform distribution
    for i = 1: num_dur
        for g = 1 : blockTrials/num_dur
            j = i+(g-1)*num_dur;
            durList(1, j) = durations(i);  
        end 
    end
    
    % parameters for U-shaped distribution
    weight_u = [30 2 2 2 2 2 2 30];
    index = 0;
    for i = 1: num_dur
        for g = 1 : weight_u(i)
            durList(2,index+g) = durations(i);  
        end 
        index = index + weight_u(i);
    end
    %currmean = sum(durations.*weight_u/blockTrials); %code to check if mean == 0.9250
    %[mu,s,muci,sci] = normfit(durList(2,:)) %sigma = 0.4940

    % parameters for normal distribution
    weight_nor = [2 2 2 30 30 2 2 2];
    index = 0;
    for i = 1: num_dur
        for g = 1 : weight_nor(i)
            durList(3, index+g) = durations(i);  
        end 
        index = index + weight_nor(i);
    end
       %currmean = sum(durations.*weight_nor/blockTrials); %code to check if mean == 0.9250
     %[mu,s,muci,sci] = normfit(durList(3,:)) %sigma = 0.1762
    
    
        
   % enquire subject information
    exp.subInfo('U: 1(y) or 2(N)?)', '1');
    if(strcmp(exp.sName, 'cancelled'))
        disp('This experiment has been cancelled!')
        return;
    end
    
    %set the order
    if(exp.sPara == 1)
        exp.seq = sortrows(exp.seq, 3);  %first U-shaped distribution
    else
        exp.seq = sortrows(exp.seq, -3);  %first normal distribution
    end
   
    
    % new two practice blocks (uniform distribution)
    prac_seq = exp.genTrials(2,[1 blockTrials], 1);
    
    % add 2 practice blocks' trials at beginning of part 1
    exp.seq = [prac_seq(1:blockTrials*pracBlockNum,:); exp.seq];
    exp.maxTrls = exp.maxTrls + blockTrials*pracBlockNum;
    
    % input device
    kb=CInput('k',[0 1],{'leftArrow','rightArrow'});
    
    % audio device
    a = CAudio; 
    
    % visual display
    % v = CDisplay('bgColor',48,'fontSize',18, 'skipSync',1); % for test
    v = CDisplay('bgColor',48,'fontSize',18,'monitorSize',22, 'fullWindow', 0);
    
    HideCursor;
    
    infoText = init_text;
    
    noisedur = 2*durations(end);
    noise = a.genPinkNoise(frequency, noisedur,-50)/noisedur; % 60 dBA
    tone =  a.genTone(frequency*4, noisedur)*1.6; % 60 dBA 
    tone_s = tone+noise; %mixing noise and tone in stumuli %307200
    tone_mix = noise;  % 2*144000
    warning = a.genTone(4000,0.02); % a warning signal
    
    blockText = ['\n There are ' int2str(pracBlockNum) ' practice blocks, and '...
    int2str(exp.maxTrls/blockTrials) ' blocks in total.\n\n'];
    v.dispText([infoText.instruction, blockText, infoText.startBlock]);
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
        tone_end_idx = durList(cond(3), cond(2))*a.freq + tone_start_idx-1;
        mix_range = tone_start_idx:tone_end_idx; % 38400
        tone_mix(:,mix_range)= tone_s(:,mix_range); %2*307200
        a.prepare(tone_mix);
        a.present;
        WaitSecs(0.3 + durList(cond(3), cond(2))); % keep same for all trials
        a.stop; %stop when the sound is finished
        %WaitSecs(0.3);
        v.flip(1);
        %acquire response
        v.dispText(infoText.question);
        key= kb.response;
        if iTrl <= blockTrials*pracBlockNum
            % feedback
            feedbacktext = infoText.short;
            if durList(cond(3), cond(2)) >= durations(floor(num_dur/2) +1)
                feedbacktext = infoText.long;
            end
            v.dispText(feedbacktext);
            WaitSecs(1);
        end
        
        % store duration and responses
        exp.setResp([durList(cond(3), cond(2)) key]);
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
    ' In this experiment, you will hear a noise after a tone . \nYour task is to judge if the noise is presented for a short or long time.\n', ...
    ' During practice blocks, there is feedback telling you the right answer after each trial.'];
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

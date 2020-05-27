%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name:          bisection_exp3.m
%%
%% Author:         Fiona Zhu
%% Description:   Main function for running the bisection study 2 to
%% investigate the influence of ensemble variance on temporal bisection 
%%
%% note: 
%% three types of sampled distributions will be tested: 
%%       a uniform, a normal, and a U-shaped distribution 
%% seven intervals: (400, 600, 800, 1000, 1200, 1400, and 1600 ms)
%% Author: Xiuna Zhu (Fiona)
%% Date: 27-12-2018
%% Contact: fiona.zhu1230@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bisection_exp3
% -------------------------------------
% This is the code for bisection task to study2 with a U-shaped distribution
%The U-shaped distribution has effectively parameters:
% a= 0.4, b=1.6, 
%gravitational balance center or offset :beta = (b+a)/2 =1
%vertical scale alpha =12 /(b-a)^{3} =  6.94444...
%	PDF =  alpha* (x-beta)^{2}
% durations = [400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0]/1000;
% weight = [1/7      1/7      1/7    1/7     1/7     1/7    1/7;  %uniformal
%          4/56     7/56     10/56  14/56   10/56   7/56    4/56;  %normal
%          12/56    9/56     6/56   2/56    6/56    9/56    12/56]; %U-shaped
%Participants have to identify the presented visual stimulus was close to the short or to the long
% note:
%--------------------------------------

try
    num_dur = 7;
    num_dur_all = 1/2*num_dur*(num_dur+1); %28
    blockTrials =num_dur_all*2; % number of trials per block 56
    exp = CExp(2,[1 blockTrials],'blockRepetition', 2, 'blockFactors', 3);
    %1.    
    %2. duration 
    %3. distribution: 1 = uniform; 2 = U-shaped; 3 = normal
    
    frequency = 1000;
    durations = [400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0]/1000;
    durList = ones(3, blockTrials);

    % parameters for uniform distribution
    for i = 1: num_dur
        for g = 1 : 8
            j = i+(g-1)*num_dur;
            durList(1, j) = durations(i);  
        end 
    end
    
    % parameters for U-shaped distribution
    weight_u = [12 9 6 2 6 9 12];
    index = 0;
    for i = 1: num_dur
        for g = 1 : weight_u(i)
            durList(2,index+g) = durations(i);  
        end 
        index = index + weight_u(i);
    end
    %currmean = sum(durations.*weight_u/56); %code to check if mean == 1
     %[mu,s,muci,sci] = normfit(durList2,:)  %sigma = 0.4671

    
    
    % parameters for normal distribution
    weight_nor = [4 7 10 14 10 7 4];
%     currmean = sum(durations.*weight_nor/56); %code to check if mean == 1
%     x_mean2 = (durations-currmean).^2;
%     sigma = sqrt(sum(x_mean2.*weight_nor/56)); %sigma = 0.3251
    index = 0;
    for i = 1: num_dur
        for g = 1 : weight_nor(i)
            durList(3, index+g) = durations(i);  
        end 
        index = index + weight_nor(i);
    end
    %[mu,s,muci,sci] = normfit(durList3,:)  %sigma = 0.3251
    

    %set the order
    exp.seq = sortrows(exp.seq,-3); 
   
    
%     % add uniform distribution at the end 
%     exp.seq = [exp.seq; exp.seq(1:244,:)];
%     exp.maxTrls = exp.maxTrls + 244;
    
    
    % new two practice blocks
    prac_seq = exp.genTrials(2,[1 blockTrials], 1);
    
    % add 2 practice blocks' trials at beginning of part 1
    exp.seq = [prac_seq(1:blockTrials*2,:); exp.seq];
    exp.maxTrls = exp.maxTrls + blockTrials*2;
    
    
    % enquire subject information
    exp.subInfo();
    if(strcmp(exp.sName, 'cancelled'))
        disp('This experiment has been cancelled!')
        return;
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
    
    HideCursor;
    
    infoText = init_text;
    
    noisedur = 2*durations(end);
    noise = a.genPinkNoise(frequency, noisedur,-50)/noisedur; % 60 dBA
    tone =  a.genTone(frequency*4, noisedur)*1.6; % 60 dBA 
    tone_s = tone+noise; %mixing noise and tone in stumuli %307200
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
        initime = GetSecs;
        [key, restime]= kb.response;
        rt = restime-initime;
        if iTrl <= blockTrials*2
            % feedback
            feedbacktext = infoText.short;
            if durList(cond(3), cond(2)) == durations((num_dur+1)/2)
                if(rand() >= 0.5)
                    feedbacktext = infoText.long;
                end
            elseif durList(cond(3), cond(2)) > durations((num_dur+1)/2)
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

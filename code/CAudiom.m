classdef CAudiom < handle
    % a class managing audio interface
    % last modified by Strongway: genNoise with randn, normalize noise
    % 2016.02.02 Add Motu device (Strongway)
    properties
        freq;
        channels;
        pins; 
        latency; 
        pHandle;
        hwName; 
        hwType;
    end
    
    methods
        function obj = CAudiom(varargin)
            %constructor: frequency, channels, latency
            p = inputParser;
            p.addParamValue('freq',96000,@isnumeric);
            p.addParamValue('channels',2,@isnumeric);
            p.addParamValue('latency', 0.001, @(x) x>0);
            p.addParamValue('pins',[1 2],@isnumeric); %define the pin number of the channels
            p.addParamValue('hardware','sndCard', @(x) any(strcmpi(x,{'sndCard','Aout','Ain','Motu','Soundflower'})));  % Bing add 'Soundflower', 18.Sep.2018
            p.parse(varargin{:});
            
            obj.freq = p.Results.freq;
            obj.channels = p.Results.channels;
            obj.latency = p.Results.latency;
            obj.hwName = p.Results.hardware; 
            obj.pins = p.Results.pins;
            switch obj.hwName
                case {'sndCard'}
                    obj.hwType = 1;
                    try
                        InitializePsychSound(1);
                        obj.pHandle = PsychPortAudio('Open', [],[],2,obj.freq,obj.channels);
                        PsychPortAudio('LatencyBias', obj.pHandle, obj.latency);
                        obj.latency = PsychPortAudio('LatencyBias', obj.pHandle);
                    catch ME
                        disp(ME.message);
                    end
                case {'Aout','Ain'}
                    obj.hwType = 2;
                    obj.pHandle = analogoutput('nidaq',obj.hwName);
                    achan = addchannel(obj.pHandle,obj.pins);
                    obj.channels = size(achan,1);
                    set(obj.pHandle,'SampleRate',obj.freq);
                    set(obj.pHandle,'TriggerType','manual');
                case {'Motu'}
                    obj.hwType = 1;
                    try
                        InitializePsychSound(1);
                        dv = PsychPortAudio('GetDevices');
                        for id = 1:length(dv)
                           if (strncmp('MOTU Audio ASIO',dv(id).DeviceName,11) &&  dv(id).NrOutputChannels > 0)
                                break;
                            end
                        end
                        obj.pHandle = PsychPortAudio('Open', dv(id).DeviceIndex,[],2,obj.freq,obj.channels);
                        PsychPortAudio('LatencyBias', obj.pHandle, obj.latency);
                        obj.latency = PsychPortAudio('LatencyBias', obj.pHandle);
                    catch ME
                        disp(ME.message);
                    end
                case {'Soundflower'}
                    obj.hwType = 1;
                    try
                        InitializePsychSound(1);
                        dv = PsychPortAudio('GetDevices');
                        for id = 1:length(dv)
                            if strcmp('Soundflower',dv(id).DeviceName)
                                break;
                            end
                        end
                        
                        obj.pHandle = PsychPortAudio('Open', dv(id).DeviceIndex,[],2,obj.freq,obj.channels);
                        PsychPortAudio('LatencyBias', obj.pHandle, obj.latency);
                        obj.latency = PsychPortAudio('LatencyBias', obj.pHandle);
                    catch ME
                        disp(ME.message);
                    end
                    
                otherwise
                    obj.hwType = 1;
            end
        end
        
        function prepare(obj, data)
            switch obj.hwType 
                case 1
                    PsychPortAudio('FillBuffer', obj.pHandle, data);
                    disp(['Prepare Audio: obj.pHandle = ', int2str(obj.pHandle)]);
                case 2
                    putdata(obj.pHandle,data');
                    start(obj.pHandle);
            end
            
        end
        
        function onsetTime = present(obj,startTime)
            if nargin < 2
                startTime = 0;
            end
            switch obj.hwType
                case 1
                    onsetTime = 0;
                    PsychPortAudio('Start', obj.pHandle, 1, startTime, 0);
                    disp(['Present Audio: obj.pHandle = ', int2str(obj.pHandle)]);
%dead loop during checking the status of the PsychPortAudio handler
%                     active = 0;
%                     while active == 0
%                         s = PsychPortAudio('GetStatus',obj.pHandle);
%                         active = s.Active;
%                     end
                    s = PsychPortAudio('GetStatus',obj.pHandle);
                    if s.Active ~= 0
                        onsetTime = s.StartTime; %get true onset time
                    end 
                case 2        
                    if startTime > 0
                        waitTime = max(0,startTime - GetSecs);
                        WaitSecs(waitTime);
                    end
                    trigger(obj.pHandle);
                    onsetTime = getSecs;
            end
        end
        
        function stop(obj,afterFinish)
            if nargin < 2
                afterFinish = 0; %immediately stop
            end
            switch obj.hwType
                case 1
                    PsychPortAudio('Stop', obj.pHandle, afterFinish);
                    disp(['Stop Audio: obj.pHandle = ', int2str(obj.pHandle)]);
                case 2
                    stop(obj.pHandle);
            end
        end
        function open(obj)
            obj.pHandle = PsychPortAudio('Open', [],[],2,obj.freq,obj.channels);
            PsychPortAudio('LatencyBias', obj.pHandle, obj.latency);
            obj.latency = PsychPortAudio('LatencyBias', obj.pHandle);
        end
        function close(obj)
            switch obj.hwType
                case 1
                    PsychPortAudio('Stop', obj.pHandle);
                    % Close the audio device:
                    PsychPortAudio('Close', obj.pHandle);
            end    
        end
        
        function tone = genTone(obj, freq, duration)
            t = MakeBeep(freq, duration, obj.freq);
            tone = repmat(t(1:end-1),obj.channels(1),1);
        end
        
        function tone = genNoise(obj,duration)   
            % change to randn (normal distribution) - Oct.14. 2013
            num_samp = obj.freq*duration;
            t = randn(1,num_samp);
            tone = repmat(t,obj.channels(1),1);
            % rescale
            tone = tone./3; % 3 sigma
        end
                        
        
        function [tone, pn] = genPinkNoise(obj, f0, duration, dbc_per_hz, num_taps)
            % based on code phase_noise by Jeff Schenck 11/21/95
            % f0 reference frequency (must be in Hz.)
            % dbc_per_hz power per hertz relative to carrier at ref. freq.
            % num_taps number of filter taps in AR 1/f filter
            % (optional; default = 100)

            % pn phase-modulated 1/f process
            % theta 1/f process (before phase modulation)

            % Check input.

            if dbc_per_hz >= 0
                error('Power per Hz. must be negative.');
            elseif f0 <= 0
                error('Reference frequency must be positive.');
            end

            if nargin < 5
                num_taps = 100;
            end

            num_samp = obj.freq * duration;
            % Generate white noise. Apply gain for desired dBc/Hz. Warn user
            % if gain is too large (gain thresholds have been chosen somewhat
            % arbitrarily -- needs work).

            gain = sqrt(2*pi * f0 * 10^(dbc_per_hz/10));
            wn = gain * randn(1,num_samp);

            fprintf('Gain applied to white noise = %f.\n', gain);
            if gain >= 1
                fprintf('WARNING: Narrowband approximation no longer valid.\n');
            elseif gain >= .5
                fprintf('WARNING: Narrowband approximation on the verge of collapse.\n');
            end


            % Generate 1/f AR filter and apply to white noise to produce 1/f
            % noise.

            a = zeros(1,num_taps);
            a(1) = 1;
            for ii = 2:num_taps
                a(ii) = (ii - 2.5) * a(ii-1) / (ii-1);
            end
            theta = filter(1,a,wn);

            tone = repmat(theta,obj.channels(1),1);
            % Phase modulate.

            pn = exp(i*theta);

            % normalize
            sd3 = std(tone(1,:))*3; % 3 sigma
            tone = tone./sd3; 
        end
        
    end
end

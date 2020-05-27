classdef CInput < handle
    % A class of peripheral input device
    % Input device:
    % 1. keyboard: CInput('k',response values, {KeyName List}). Keys are
    % mapped into response values. For example, 
    %   kb = CInput('k',[1 2],{'LeftArrow','RightArrow'); % Left arrow 1,   right arrow 2
    %   by default: kb=CInput('k'); % Left and Right Arrow
    % 2. Mouse: CInput('m', [1 2]);% left/right mouse button to 1 ,2
    % 3. Parallel port key: CInput('p', [1 2]);
    %       note: parallel keypad encoding depends on pins, so be checked
    %       before use this method
    % 4. NI input: CInput('n',defaultValues, pins, portname);
    %       e.g. nk = CInput('n', [0 0], [8 9],'AOut');
    %           Card 'AOut', pin 8 and 9, default values [0 0]
    % Methods:
    % 1. wait(); %wait for input
    % 2. [kinput ktime x y] = response(bWait); %by default, wait for response. 
    %   when bWait == 0, it samples reponse values only. 
    % 3. releaseTime = keyRelease();
    %      return releaseTime
    % last modify: 15. March, 2012
    % First created: 2nd April, 2001
    % created by: Z. Shi, shi@lmu.de
    properties
        esckey = 41;
        key = [80 79];
        kval = [1 2];
        ktype = 1;
        wantStop = 0;
        dio ;
    end
    
    methods 
        function obj = CInput(keytype,keyvalues,keys,portname)
            try
                KbName('UnifyKeyNames');
                obj.esckey = KbName('ESCAPE');
                switch lower(keytype)
                    case {'keyboard','k'}
                        obj.ktype = 1;
                        if nargin < 3
                            obj.key = [KbName('LeftArrow') KbName('RightArrow')];
                        end
                        if nargin < 2
                            keyvalues = [1 2];
                        end
                        if nargin == 3
                            for k = 1:length(keys)
                                obj.key(k) = KbName(keys{k});
                            end
                        end
                            obj.kval= keyvalues;
                    case {'mouse','m'}
                        obj.ktype = 2;
                        if nargin < 2
                            keyvalues = [1 2];
                        end
                        obj.kval = keyvalues;

                    case {'parallel','p'}
                        obj.ktype = 3;
                        % todo
                        obj.key = keys; % actually parallel port pins
                        obj.kval = keyvalues;
                        obj.dio = digitalio('parallel', 'LPT1');
                        addline(obj.dio,obj.key,'in');

                    case {'n','nidaq'}
                        obj.ktype = 4;
                        obj.key = keys; % pins
                        obj.kval = keyvalues; % default key values
                        obj.dio = digitalio('nidaq',portname);
                        addline(obj.dio,obj.key,'in');
                        
                    otherwise % default as keyboard
                        obj.ktype = 1;
                        obj.key(1) = KbName('LeftArrow');
                        obj.key(2) = KbName('RightArrow');
                        obj.kval = [1 0];
                end
            catch ME
                disp(ME.message);
            end
        end
        
        function wait(obj)
            % halt program for key input 
            switch obj.ktype
                case 1 %keyboard
                    while KbCheck; end %clear buffer
                    KbWait;
                case 2 % mouse
                    GetClicks;
                case 3 % parallel port keypad
                    % need to test for different parallel port pad
                    %here assume: initial status: 1111 ...
                    k = getvalue(obj.dio);
                    while sum(k) == length(k)
                        k = getvalue(obj.dio);
                    end
                case 4 %ni daq
                    k = getvalue(obj.dio);
                    while sum(xor(obj.kval,k)) == 0 % sum(k) == length(k)
                        k = getvalue(obj.dio);
                    end                    
            end
        end
        
        %
        function [kinput ktime x y] = response(obj,bWaitPress)
            % wait for response (by default, it waits for response)
            % if bWaitPress == 0, it samples response values
            % return input keys, time, and x, y
            x = -1;
            y = -1;
            if nargin < 2
                bWaitPress = 1;
            end
            validKey = 0;
            kinput = -1;
            ktime = GetSecs;
            switch obj.ktype
                case 1 %keyboard
                    while 1
                        [ keyIsDown, seconds, keyCode ] = KbCheck;
                        if keyIsDown
                            if keyCode(obj.esckey)
                                kinput = -1;
                                ktime = GetSecs;
                                obj.wantStop = 1;
                                break;
                            end
                            
                            for k = 1: length(obj.key)
                                if keyCode(obj.key(k))
                                    kinput = obj.kval(k);
                                    ktime = GetSecs;
                                    validKey = 1;
                                end
                            end
                            if validKey
                                break;
                            end
                            while KbCheck; end %clear buffer
                        end
                        if ~bWaitPress
                            break;
                        end
                    end
                case 2 % mouse
                    while 1
                      [x,y,buttons] = GetMouse;
                      if buttons(1)
                          kinput = obj.kval(1);
                          ktime = GetSecs;
                          break;
                      end
                      if length(buttons)>=2 & buttons(end) 
                          kinput = obj.kval(2);
                          ktime = GetSecs;
                          break;
                      end
                      % stop running
                      [ keyIsDown, seconds, keyCode ] = KbCheck;
                        if keyIsDown && keyCode(obj.esckey)
                            kinput = -1;
                            ktime = GetSecs;
                            obj.wantStop = 1;
                            break;
                        end
                        if ~bWaitPress
                            break;
                        end
                    end 
                case 3 % parallel port keypad
                    k = getvalue(obj.dio);
%                    while sum(k) == length(k)
                    while sum(xor(obj.kval,k)) == 0
                        k = getvalue(obj.dio);
                        if ~bWaitPress
                            break;
                        end
                      % stop running
                      [ keyIsDown, seconds, keyCode ] = KbCheck;
                        if keyIsDown && keyCode(obj.esckey)
                            kinput = -1;
                            ktime = GetSecs;
                            obj.wantStop = 1;
                            break;
                        end
                    end
%                    kidx = find(k == 0);
                    ktime = GetSecs;
%                    kinput = obj.kval(kidx(1));
                    kinput = sum(find(xor(obj.kval,k)));
                case 4 %NI
                    k = getvalue(obj.dio);
                    while sum(xor(obj.kval,k)) == 0 %sum(k) == length(k)
                        k = getvalue(obj.dio);
                        if ~bWaitPress
                            break;
                        end
                    end
                  % stop running
                  [ keyIsDown, seconds, keyCode ] = KbCheck;
                    if keyIsDown && keyCode(obj.esckey)
                        kinput = -1;
                        ktime = GetSecs;
                        obj.wantStop = 1;
                    end
                    %kidx = find(k == 0);
                    ktime = GetSecs;
                    kinput = sum(find(xor(obj.kval,k)));
                    %kinput = obj.kval(kidx(1));
            end % end of switch

            
        end % end of function response
 
        function releaseTime = keyRelease(obj)
         %return release time  when key is released,    
            switch obj.ktype
                case 1 %keyboard
                    while 1
                        [ keyIsDown, seconds, keyCode ] = KbCheck;
                        if ~keyIsDown
                                break;
                        end
                    end
                case 2 % mouse
                    while 1
                      [x,y,buttons] = GetMouse;
                      if sum(buttons) == 0
                          break;
                      end
                    end
                case 3 % parallel port keypad
                    k = getvalue(obj.dio);
                    while sum(k) < length(k)
                        k = getvalue(obj.dio);
                    end
                case 4
                    %keys are too sensitive, we need to apply filter mechanisms
                    noise_ticks = 0;
                    while 1
                        k = getvalue(obj.dio);
                        if sum(xor(obj.kval,k))==0 % key releases
                            noise_ticks = noise_ticks +1;
                        else 
                            noise_ticks = 0;
                        end
                        if noise_ticks >10 % 10 ms key releases
                            break;
                        end
                        WaitSecs(0.001);
                    end
            end % end of switch
            releaseTime = GetSecs-0.01; %- 10 ms
            
            % stop running
            [ keyIsDown, seconds, keyCode ] = KbCheck;
            if keyIsDown && keyCode(obj.esckey)
                obj.wantStop = 1;
            end
            
        end %end of function keyResease
        
        
    end
end
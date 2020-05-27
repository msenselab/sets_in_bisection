classdef CDisplay < handle
    % a display class, manage the most common display functions, such as
    % display text, fixation, flip, degree of visual angle. It also stores
    % display informations, such as resolution, center x,y, ifi ...
    % methods:
    %   1. constructor
    %       obj=CDisplay('monitorSize',size,'viewDistance',viewDistance,'bgColor',bgColor,'Color',color,'fontSize',fz,'skipSync',bSync);
    %   2. display text
    %       dispText(txt[,flip,clearBackground]);
    %   3. close screen
    %       close();
    %   4. calculate visual angle (deg) to pixels
    %       deg2pix(degree)
    %   5. display fixation
    %       dispFixation(size [,type, flip, clearBackground]);
    %   6. flip: bring back buffer to front (display)
    %       flip(obj [, clearBackground])    
    % Screen class of PTB
    % Last modify: 7.7.2011
    % Created by: Z. Shi, shi@lmu.de
    
    % 14.07.2011 add createShape function to create simple shapes
    % 07.07.2014 add debug option - fullWindow, 
    % 2015 add gabor and gratings

    properties
        wnd = -1;   %window handle
        bSkipSync = 0;
        fullWindow = 1; % full screen
        ifi = -1;   %inter-frame interval (refresh rate)
        cx = -1;    % center x
        cy = -1;    % center y
        pdeg;       % pixels per degree
        bgColor = 0;     % background color
        color = 255;     % front color
        fontSize = 14;
        lineWidth = 60;
        lineSpace = 1.5;
        inch = 20;
        viewDistance = 57;
        nItem = 0;
        % image parameters
        % for gabor: phase , freq , sigma, contrast, aspectRatio, ...
        % for grating: [phase, freq, contrast, 0]
        iPara = []; % 
        items; 
    end
    
    methods
        function obj = CDisplay(varargin) 
        	% constructure
            % inch,viewDistance, bgColor, color, fontsize
            p = inputParser;
            p.addParamValue('monitorSize',20,@isnumeric);
            p.addParamValue('viewDistance',57,@isnumeric);
            p.addParamValue('bgColor',[0 0 0],@isnumeric);
            p.addParamValue('Color',[255 255 255],@isnumeric);
            p.addParamValue('fontSize',14,@isnumeric);
            p.addParamValue('lineWidth',60,@isnumeric);
            p.addParamValue('lineSpace',1.5,@isnumeric);
            p.addParamValue('skipSync',0,@isnumeric);
            p.addParamValue('fullWindow',1,@isnumeric);
            p.parse(varargin{:});
            
            %init screens
            obj.inch = p.Results.monitorSize;
            obj.viewDistance = p.Results.viewDistance;
            obj.bgColor = p.Results.bgColor;
            obj.color = p.Results.Color;
            obj.fontSize = p.Results.fontSize;
            obj.lineWidth = p.Results.lineWidth;
            obj.lineSpace = p.Results.lineSpace;
            obj.bSkipSync = p.Results.skipSync;
            obj.fullWindow = p.Results.fullWindow;
           try
                InitializeMatlabOpenGL;
                AssertOpenGL;
%                 priority = MaxPriority('KbCheck');
%                 Priority(priority);
               HideCursor; 
                if obj.bSkipSync
                    Screen('Preference','SkipSyncTests',1);
                else
                    Screen('Preference','SkipSyncTests',0);
                end

                screens=Screen('Screens');
                screenNumber=max(screens);
                PsychImaging('PrepareConfiguration');
                PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
                if obj.fullWindow
%                 [obj.wnd wsize] = Screen('OpenWindow',screenNumber); % Open On Screen window, mainWnd
                    [obj.wnd, wsize] = PsychImaging('OpenWindow', screenNumber, obj.bgColor);
                else
                    [obj.wnd, wsize] = PsychImaging('OpenWindow', screenNumber, obj.bgColor, [0 0 600 400]);
%                    [obj.wnd wsize] = Screen('OpenWindow',screenNumber,[],[0 0 600 400]); 
                end
                obj.ifi=Screen('GetFlipInterval', obj.wnd); 
                Screen('BlendFunction', obj.wnd, 'GL_ONE', 'GL_ZERO');
                
                obj.cx = wsize(3)/2; %center x
                obj.cy = wsize(4)/2; %center y
                pix = obj.inch*2.54/sqrt(1+9/16)/wsize(3);  % calculate one pixel in cm
                obj.pdeg = round(2*tan((1/2)*pi/180) * obj.viewDistance / pix); 
            catch ME
                Screen('CloseAll');
                Priority(0);
                disp(ME.message);
                disp('error in initial display');
            end
        end % end of constructor

        function nframes = sec2frames(obj,secs)
            %convert seconds to number of frames
            nframes = round(secs/obj.ifi);
        end
        
        function dispText(obj,txt, flip,clearBackground)
        % dispText
            try
                if nargin < 3  
                    flip = 1;
                    clearBackground = 1;
                end
                if nargin == 3
                    clearBackground = 1;
                end
                if clearBackground 
                	Screen('FillRect',obj.wnd,obj.bgColor);
                end
                Screen('TextSize',obj.wnd,obj.fontSize);
                DrawFormattedText(obj.wnd,txt,'center','center',obj.color,obj.lineWidth,[],[],obj.lineSpace);
                if flip
                    Screen('Flip', obj.wnd);
                end  
            catch ME
                Screen('CloseAll');
                Priority(0);
                disp(ME.message);
            end
        end % end of dispText method
        
        function itemIndex = createItem(obj,itemData)
            %create items (texture, image, etc.)
             itemIndex= Screen('MakeTexture', obj.wnd, itemData);
             obj.nItem = obj.nItem + 1;
             obj.items(obj.nItem) = itemIndex;
        end
        
        function updateGratingPara(obj,n, varargin)
        % To modify parameters efficiently, directly access the
        % obj.iPara;
            p = inputParser;
            p.addParameter('phase',0);
            p.addParameter('freq',1);
            p.addParameter('contrast',0.5);
            p.parse(varargin{:});
            if length(p.Results.phase(:))< n
                phase = repmat(p.Results.phase(:),n,1);
                phase = phase(1:n);
            else
                phase = p.Results.phase(1:n);
            end
            if length(p.Results.freq(:))< n
                freq = repmat(p.Results.freq(:),n,1);
                freq = freq(1:n);
            else
                freq = p.Results.freq(1:n);
            end
            
            if length(p.Results.contrast(:))< n
                contrast = repmat(p.Results.contrast(:),n,1);
                contrast = contrast(1:n);
            else
                contrast = p.Results.contrast(1:n);
            end
            
            % convert frequency from degree to pix
            obj.iPara = [phase(:), freq(:), contrast(:), zeros(n,1)]';
        end
        
        function updateGaborPara(obj,n, varargin)
        % note: inputParser is relative slow, 5 parameters requires 15 ms
        % etc. So use if carefully, only for initialization. 
        % To modify parameters efficiently, directly access the
        % obj.iPara;
        % when sigma --> inf, gabor becomes grating. 
            p = inputParser;
            p.addParameter('phase',0);
            p.addParameter('freq',1);
            p.addParameter('sigma',2);
            p.addParameter('contrast',0.5);
            p.addParameter('aspectRatio',1);
            p.parse(varargin{:});
            if length(p.Results.phase(:))< n
                phase = repmat(p.Results.phase(:),n,1);
                phase = phase(1:n);
            else
                phase = p.Results.phase(1:n);
            end
            if length(p.Results.freq(:))< n
                freq = repmat(p.Results.freq(:),n,1);
                freq = freq(1:n);
            else
                freq = p.Results.freq(1:n);
            end
            if length(p.Results.sigma(:))< n
                sigma = repmat(p.Results.sigma(:),n,1);
                sigma = sigma(1:n);
            else
                sigma = p.Results.sigma(1:n);
            end
            
            if length(p.Results.contrast(:))< n
                contrast = repmat(p.Results.contrast(:),n,1);
                contrast = contrast(1:n);
            else
                contrast = p.Results.contrast(1:n);
            end
            if length(p.Results.aspectRatio(:))< n
                aspectRatio = repmat(p.Results.aspectRatio(:),n,1);
                aspectRatio = aspectRatio(1:n);
            else
                aspectRatio = p.Results.aspectRatio(1:n);
            end
            
            obj.iPara = [phase(:), freq(:), sigma(:), contrast(:), ...
                aspectRatio(:), zeros(n,3)]';
        end
        
        function itemIndex = createGrating(obj,x,y,varargin)
            p = inputParser;
            p.addRequired('x',@(x) x>0);
            p.addRequired('y',@(x) x>0);
            p.addParameter('bgOffset',[0.5, 0.5, 0.5, 0]);
            p.parse(x,y,varargin{:});
            bgColorOffset = p.Results.bgOffset;
            xp = round(p.Results.x * obj.pdeg);
            yp = round(p.Results.y * obj.pdeg); %convert to pixels
            
            itemIndex = CreateProceduralSineGrating(obj.wnd, xp, yp,...
                bgColorOffset, [], 0.5);
        end
            
        function itemIndex = createGabor(obj, x,y, varargin)
            p = inputParser;
            p.addRequired('x',@(x) x>0);
            p.addRequired('y',@(x) x>0);
            p.addParameter('bgOffset',[0.5, 0.5, 0.5, 0]);
            p.parse(x,y,varargin{:});
            bgColorOffset = p.Results.bgOffset;
            xp = round(p.Results.x * obj.pdeg);
            yp = round(p.Results.y * obj.pdeg); %convert to pixels
            
            itemIndex = CreateProceduralGabor(obj.wnd, xp, yp,...
                [], bgColorOffset, 1, 0.5);
        end
 
        function itemIndex = createGaborTex(obj,width,freq,tilt)
            % this function is backward compatible if psychimaging is not
            % working
            if nargin<4
                tilt = 0;
            end
            if nargin<3
                freq = 2;
            end
              % visual gabor stimuli  
            r = obj.pdeg * width;
            widthArray = -round(r/2): round(r/2);
            pixelsPerPeriod = 1/freq * obj.pdeg; % How many pixels will each period/cycle occupy?
            spatialFrequency = 1 / pixelsPerPeriod; % How many periods/cycles are there in a pixel?
            radiansPerPixel = spatialFrequency * (2 * pi); % = (periods per pixel) * (2 pi radians per period)
            gaussianSpaceConstant = r/4; %1.5  * pixelsPerPeriod;
            [x, y] = meshgrid(widthArray, widthArray);
            tiltInRadians = tilt * pi / 180;
            a=cos(tiltInRadians)*radiansPerPixel;
            b=sin(tiltInRadians)*radiansPerPixel;

            gratingMatrix = sin(a*x+b*y);
            circularGaussianMaskMatrix = exp(-((x .^ 2) + (y .^ 2)) / (gaussianSpaceConstant ^ 2));
            imageMatrix = gratingMatrix .* circularGaussianMaskMatrix;
            grayscaleImageMatrix = 128 + 128 * imageMatrix;

            %create texture
            itemIndex= Screen('MakeTexture', obj.wnd, grayscaleImageMatrix);
            obj.nItem = obj.nItem + 1;
            obj.items(obj.nItem) = itemIndex;
            
        end

        
        function itemIndex = createShape(obj,name,x,y,varargin)
            %create simple shape 
            %be sure to do this before the trial starts
            p = inputParser;
            p.addRequired('name', @(x) any(strcmpi(x,{'rectangle','circle'})));
            p.addRequired('x',@(x) x>0);
            p.addRequired('y',@(x) x>0);
            p.addParameter('fill',1,@isnumeric);
            p.addParameter('border',0.1,@(x) x>0);
            p.addParameter('bgColor',obj.bgColor,@isnumeric);
            p.addParameter('color',obj.color,@isnumeric);
            p.parse(name,x,y,varargin{:});
            
            xp = round(p.Results.x * obj.pdeg)/2;
            yp = round(p.Results.y * obj.pdeg)/2; %convert to pixels
            bp = round(p.Results.border * obj.pdeg);
            bc = p.Results.bgColor;
            fc = p.Results.color;
            if length(bc) == 1
                bc = repmat(bc,3,1);
            end
            if length(fc) == 1
                fc = repmat(fc,3,1);
            end
            data = zeros(xp*2,yp*2,3); %store in rgb format
            switch p.Results.name
                case {'circle'}
                    if p.Results.fill == 1 %fill
                        for ix = 1:xp*2
                            for iy = 1:yp*2
                                if (ix-xp)*(ix-xp)/xp/xp + (iy-yp)*(iy-yp)/yp/yp < 1 
                                    data(ix,iy,:) = fc;
                                else
                                    data(ix,iy,:) = bc;
                                end
                            end
                        end
                    else % frame
                        for ix = 1:xp*2
                            for iy = 1:yp*2
                                if (ix-xp)*(ix-xp)/xp/xp + (iy-yp)*(iy-yp)/yp/yp < 1 && ...
                                        (ix-xp)*(ix-xp)/(xp-bp)/(xp-bp) + (iy-yp)*(iy-yp)/(yp-bp)/(yp-bp) >= 1 
                                    data(ix,iy,:) = fc;
                                else
                                    data(ix,iy,:) = bc;
                                end
                            end
                        end
                    end
                case {'rectangle'}
                    if p.Results.fill == 1 %fill
                        for ix = 1:xp*2
                            for iy = 1:yp*2
                                    data(ix,iy,:) = fc;
                            end
                        end
                    else %frame
                        for ix = 1:xp*2
                            for iy = 1:yp*2
                                if abs(ix-xp)>xp-bp || abs(iy-yp)>yp-bp
                                    data(ix,iy,:) = fc;
                                else
                                    data(ix,iy,:) = bc;
                                end
                            end
                        end
                    end
            end %end of switch
            %create texture
            itemIndex= Screen('MakeTexture', obj.wnd, data);
            obj.nItem = obj.nItem + 1;
            obj.items(obj.nItem) = itemIndex;
            
        end
        function dispItems(obj, xys, itemIndex, itemSizes,rotations, flip)
            %disp items at xys (in visual angle, center 0,0)
            if nargin < 6
                flip = 1;
            end
            if nargin < 5
                rotations = [];
            end
            if nargin < 4
                itemSizes = [obj.cx / obj.pdeg, obj.cy/obj.pdeg]*2;
            end
            destRects = zeros(4, length(itemIndex));

            for iObj = 1: size(itemSizes,1)
                if size(itemSizes,1) == 1
                    itemRect = [0 0 itemSizes*obj.pdeg];
                else
                    itemRect = [0 0 itemSizes(iObj,:)*obj.pdeg];
                end
                rect = CenterRectOnPoint(itemRect,obj.pdeg*xys(iObj,1)+obj.cx, obj.pdeg*xys(iObj,2)+obj.cy);
                destRects(:,iObj) = rect';
            end

            % convert freq, gama from deg to pixel
            switch size(obj.iPara,1) 
                case 8 % gabor
                    param = obj.iPara;
                    param(2, :) = param(2,:)./obj.pdeg;
                    param(3, :) = param(3,:).*obj.pdeg;
                    Screen('DrawTextures',obj.wnd,itemIndex,[],destRects,rotations, ...
                        [], [], [], [], kPsychDontDoRotation, param);
                case 4 % grating
                    param = obj.iPara;
                    param(2, :) = param(2,:)./obj.pdeg;
                    Screen('DrawTextures',obj.wnd,itemIndex,[],destRects,rotations, ...
                        [], [], [], [], [], param);
                otherwise
                    Screen('DrawTextures',obj.wnd,itemIndex,[],destRects,rotations);
            end                    
            if flip
                Screen('Flip', obj.wnd);
            end
        end
        
        function close(obj)
        %% dispClose: close display
            %close display and delete Screen
            Screen('CloseAll');
            ShowCursor;
            Priority(0);
            obj.wnd = -1;
        end % end of closeDisp
        
        function pixs=deg2pix(obj, degree) 
        %% deg2pix: calculate degree to pixels
            screenWidth = obj.inch*2.54/sqrt(1+9/16);  % calculate screen width in cm
            pix=screenWidth/obj.cx/2;  %calculates the size of a pixel in cm 
            pixs = round(2*tan((degree/2)*pi/180) * obj.viewDistance / pix); 
        end %end of deg2pix
        
        function deg = pix2deg(obj, pixels)
            deg = 1/obj.pdeg*pixels;
        end
        
        function dispFixation(obj,sz,type,flip,clearBackground)
        %% dispCrossFix: display cross fixation
            if nargin < 3
                type = 1;
                flip = 1;
                clearBackground = 1;
            end
            if nargin == 3
                flip = 1;
                clearBackground = 1;
            end
            if nargin == 4
                clearBackground = 1;
            end
            if clearBackground
                Screen('FillRect',obj.wnd, obj.bgColor);
            end
            if type == 1
                rect = CenterRectOnPoint([0 0 2 sz],obj.cx,obj.cy);
                Screen('FillRect',obj.wnd,obj.color,rect);
                rect = CenterRectOnPoint([0 0 sz 2],obj.cx,obj.cy);
                Screen('FillRect',obj.wnd,obj.color, rect);
            else
                rect = CenterRectOnPoint([0 0 sz sz],obj.cx,obj.cy);
                Screen('FillOval',obj.wnd,obj.color,rect);
            end
            if flip
                Screen('Flip',obj.wnd);
            end
        end %end of dispCrossFix
                        
        function [vbl visonset] = flip(obj,clearBackground, when)
        %% flip: put back buffer to screen
            if nargin < 3
                when = 0;
            end
            if nargin < 2
                clearBackground = 0;
            end
            if clearBackground
                Screen('FillRect',obj.wnd,obj.bgColor);
            end
            [vbl visonset] = Screen('Flip',obj.wnd, when);
        end
        
    end
end
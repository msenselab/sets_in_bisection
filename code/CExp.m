classdef CExp < handle
% class of experiment
% store experiment data, structure, procedure etc
% 10 July, 2011
% Created by: Z. Shi, shi@lmu.de
% ver. 0.2
% 21 June, 2012
%   add option: Continue. If will continue run experiment with given subject name
%   add additional Data variable
%
%constructor: two ways to initialize: 
    % constant method: 
    %   p1, repetition, 
    %   p2, factors
    % adapative method:
    %   p1, alpha 
    %   p2, beta
    %optional parameters: blockReptition, blockFactors;  
    % adaptive method:maxTrials,maxRespRev, gamma, lambda

    % add UML method in the class 11.Feb.2015
    % initUML(obj, x0,xrange, alpha, n_alpha, beta, n_beta)
    properties
        sName = 'test'; %subject information
        sAge = 20;
        sGender = 'f';
        sPara; %additional parameter
        eType = 1; % 1 - constant stimuli; 2- bayesian adaptive method ; 3 - UML       
        % type 1
        seq;    % sequences of the trials
        resp;   % corresponding responses
        aData;  % store additional data, can be anything...
        % only for type 2
        beta;
        alpha;
        gamma;
        lambda;
        maxTrls; %maximum trial number
        maxRev; % maximum reversal number
        % index
        curTrl; %current trial index
        curIntensity;
        practice;
        % other information
        startTime;
        endTime;
        
        % uml methods
        u; % uml
        upara; % uml configuration
     end
    
    methods
        function obj = CExp(p1,p2,varargin)
            % constant method: 
            %   p1, repetition, 
            %   p2, factors
            % adapative method:
            %   p1, alpha 
            %   p2, beta
            %optional parameters: blockRepetition, blockFactors;  
            % adaptive method:maxTrials,maxRespRev, gamma, lambda
            p = inputParser;
            p.addOptional('eType','c', @(x) any(strcmpi(x,{'c','a'})));
            p.addParameter('blockRepetition',1,@(x) x>0);
            p.addParameter('blockFactors', [], @isnumeric);
            p.addParameter('maxTrials', 40, @isnumeric);
            p.addParameter('maxRespRev',8, @isnumeric);
            p.addParameter('gamma',0.025,@(x) min(x)>=0);
            p.addParameter('lambda',0.025, @(x) x>=0);
            p.addParameter('practice',5, @(x) x>=0);
            p.parse(varargin{:});
            
            
            switch p.Results.eType 
                case {'c'} %constant method or uml method
                    obj.eType = 1;
                    obj.curTrl = 1;
                    obj.seq = obj.genTrials(p1,p2,p.Results.blockRepetition,p.Results.blockFactors);
                    obj.resp = [];
                    obj.maxTrls = length(obj.seq);
                

                case 'a' % adaptive method
                    obj.alpha = p1;
                    obj.beta = p2;
                    obj.gamma = p.Results.gamma;
                    obj.lambda = p.Results.lambda;
                    obj.maxTrls = p.Results.maxTrials;
                    obj.maxRev = p.Results.maxRespRev;
                    obj.practice = p.Results.practice;
                    obj.eType = 2;
                    obj.curTrl = 1;
                    obj.seq = zeros(obj.maxTrls,2);
                    obj.resp = [];

            end
        end
        
        function initUML(obj, x0,xrange, n_samples )
            % initialize UML parameters.
            % x0: initial guess
            % xrange: the limits to the signal strength [0,1]
            % n_samples: number of samples for alphas and beta(e.g., 101)
            
            par.model = 'logit';    
            par.ndown = 1;  % the parameter for the up-down sweetpoint selection rule
            par.method = 'mean';
            par.x0 = x0;    % the initial signal strength
            par.x_lim = xrange;   % the limits to the signal strength
            
            beta0 = [0.5, 100]/(xrange(2)-xrange(1)); % a simplify beta specification, including almost all flat to steep slopes

            par.alpha = struct(...
                'limits',xrange,...       %range of the parameter space for alpha
                'N',n_samples,...                %number of alpha values. If this value is set to 1, then the first element of alpha.limits would be the assumed alpha and the alpha parameter is not estimated.
                'scale','lin',...         %the linear or log spacing. Choose between 'lin' and 'log'.
                'dist','flat',...         %prior distribution of the alpha parameter. Choose between 'norm' and 'flat'.
                'mu',10,...                %mean of the prior distribution.
                'std',5 ...              %standard deviation of the prior distribution.  
                );

            par.beta = struct(...
                'limits',beta0,...      %range of the parameter space for beta
                'N',n_samples,...                %number of beta values. If this value is set to 1, then the first element of beta.limits would be the assumed beta and the beta parameter is not estimated.
                'scale','log',...         %the linear or log spacing. Choose between 'lin' and 'log'.
                'dist','flat',...         %prior distribution of the beta parameter. Choose between 'norm' and 'flat'.
                'mu',1,...                %mean of the prior distribution.
                'std',5 ...               %standard deviation of the prior distribution.
                );

            par.gamma = 0;

            par.lambda = struct(...
                'limits',[0 0.1],...      %range of the parameter space for lambda
                'N',1,...                 %number of lambda values. If this value is set to 1, then the first element of lambda.limits would be the assumed lambda and the lambda parameter is not estimated.
                'scale','lin',...         %the linear or log spacing. Choose between 'lin' and 'log'.
                'dist','flat',...         %prior distribution of the lambda parameter. Choose between 'norm' and 'flat'.
                'mu',0,...                %mean of the prior distribution.
                'std',0.1 ...             %standard deviation of the prior distribution.  
                );
           
            obj.upara = par; % need parameters from exp_config();
            %create multiple UML instances
            nfactors = 1;
            for i=1:size(obj.seq,2)
               nfactors = nfactors * length(unique(obj.seq(:,i))); 
            end
            for i=1:nfactors
                obj.u{i} = UML(obj.upara);
            end
        end
        
        function obj = guessThreshold(obj,x)
            obj.seq(1,1:2) = [x 0];
            obj.curIntensity = x;
        end
        
        function obj = setResp(obj,resp)
            % store responses data
            % when UML method is used, resp should be 2 values [idx_uml, key]
            obj.resp(obj.curTrl,1:size(resp,2)) = resp;
            obj.curTrl = obj.curTrl+1;
            
            if ~isempty(obj.u) % uml method
                %update adaptive function
                % resp should contain (uml_id, key) % key should be 1, 0
                obj.u{resp(1)}.update(resp(2));
            end
            
        end
        
        function curSeq = getCondition(obj)
            curSeq = obj.seq(obj.curTrl,:);
        end
        
        function x = xNext(obj, id)
            if ~isempty(obj.u)
            % get next x from uml object (only valid for uml method
                x = obj.u{id}.xnext;
            end
        end
        
            
        function obj = updateThreshold(obj, p_target)
            % logistic function
            % p = Logistic(x, alpha, beta, gamma, lambda)
            %       alpha - threshold
            %       beta - slope
            %       gamma - chance performance / fa
            %       lambda - lapsing rate
            % Based on MLP toolbox
            % updated intensity and fa are stored in obj.seq
            if obj.curTrl < obj.practice
                % practice
                ra = randperm(length(obj.alpha));
                obj.seq(obj.curTrl,1)= obj.alpha(ra(1));
                obj.seq(obj.curTrl,2) = obj.gamma(1);
            else
                ll=zeros(length(obj.alpha), 1);
                x = obj.seq(obj.practice:max(obj.curTrl-1,1),1);
                responses = obj.resp(obj.practice:max(obj.curTrl-1,1));

                % calculate the likelihood of each psychometric function
                for i=1:length(obj.alpha)
                    for j=1:length(obj.gamma)
                        ll(i, j)=CalculateLikelihood(obj,x, responses, obj.alpha(i), obj.gamma(j));
                    end
                end

                % find the most likely psychometric function
                [i, j]=find(ll==max(max(ll)));
                if length(i)+length(j) > 2
                    i = i(1);
                    j = j(1);
                end;
                % calculate the level of the stimulus at p_target performance
                % using inverse logistic function
                obj.curIntensity = obj.alpha(i)-(1/obj.beta)*log(((1-obj.lambda-obj.gamma(j))./(p_target-obj.gamma(j)))-1);
                obj.seq(obj.curTrl,1)= obj.curIntensity;
                obj.seq(obj.curTrl,2) = obj.gamma(j);
            end
        end

        function bFinish = canStop(obj)
            if obj.eType == 2
                bFinish = 0;
                revnum = sum(abs(diff(obj.resp(obj.practice:obj.curTrl-1,1))));
                if revnum > obj.maxRev
                    bFinish = 1;
                end
                if obj.curTrl > obj.maxTrls
                    bFinish = 1;
                end
            else
                disp('This is only available for adaptive method');
            end
        end
        
        function ll=CalculateLikelihood(obj,x, responses, alpha, gamma)
            if obj.eType == 2
                warning off
                ll = 0;
                % calculate logistic probablity
                p=gamma+((1-obj.lambda-gamma).*(1./(1+exp(obj.beta.*(alpha-x)))));

                ll = sum(log(p(responses ==1)))+ sum(log(1-p(responses == 0)));
            else
                disp('This is only available for adaptive method');
            end
        end
        
        function h = plotLogistic(obj, a, g)
            %plot a logistic function for specify parameter, only available
            %for adaptive method
            if obj.eType == 2
                x = obj.alpha;
                if nargin < 2
                    a = median(obj.alpha);
                    g = median(obj.gamma);
                end
                y = g+((1-obj.lambda-g).*(1./(1+exp(obj.beta.*(a-x)))));
                h = figure; 
                plot(x,y);
            else
                disp('This is only available for adaptive method');
                h = -1;
            end
        end
        
        function trials=genTrials(obj,withinblock_repetitions, withinblock_factors, betweenblock_repetitions, betweenblock_factors)
            % syntax: genTrials(withinblock_repetition, betweenblock_repetition, withinblock_factors, [ betweenblock_factors])
            % eg: genTrials(2, [2 3],10); generate two factors (2 levels and 3 levels resp. ) with within block repetition 2 and between block repetition 10
            % coded by: strongway
            rand('seed',sum(clock*100));
            if  nargin < 3
                    error('Incorrect number of arguments for genTrials function');
            elseif nargin == 3
                    betweenblock_repetitions = 1;
                    betweenblock_factors  = [];
            elseif nargin == 4
                % when there's only betweenblock_repetitions, which is equal to swap
                % between block repetitions and factors -strongway 14. Dec. 2006
                    betweenblock_factors = betweenblock_repetitions;
                    betweenblock_repetitions  = 1;
            end

            trials = [];
            block_design = [];
            numblock = betweenblock_repetitions;
            if ~isempty(betweenblock_factors)
                block_design = fullfact(betweenblock_factors);
                block_design = repmat(block_design, betweenblock_repetitions,1);
                idxb = randperm(length(block_design));
                block_design = block_design(idxb,:);
                numblock = length(block_design);
            end

            for iblock = 1:numblock
                %generate within block trials
                inblock_trials = fullfact(withinblock_factors);
                inblock_trials = repmat(inblock_trials,withinblock_repetitions,1);
                idx=randperm(size(inblock_trials,1));
                inblock_trials = inblock_trials(idx,:);
                if ~isempty(block_design)
                    %add between factors
                    blockwise_factors = repmat(block_design(iblock,:),length(inblock_trials),1);
                    inblock_trials = [inblock_trials blockwise_factors];
                end
                trials = [trials; inblock_trials];
            end
        end

        function obj = subInfo(obj,sp1,sp2)
            % Get Subject information, sp1 and sp2 for additional parameters caption and parameters.
            promptParameters = {'Subject Name', 'Age', 'Gender (F or M?)'};
            defaultParameters = {'test', '20','F'};
            if nargin == 2
                promptParameters = [promptParameters, sp1];
                defaultParameters = {'test', '20','F','[]'};
            end
            if nargin == 3
                promptParameters = [promptParameters, sp1];
                defaultParameters = {'test', '20','F',sp2};
            end
            sub = inputdlg(promptParameters, 'Subject Info  ', 1, defaultParameters); 
            if isempty(sub)
                obj.sName = 'cancelled';
                return;
            end
            
            obj.sName = sub{1};
            obj.sAge = eval(sub{2});
            obj.sGender = sub{3};
            if nargin >= 2
                obj.sPara = eval(sub{4});
            end
            
            obj.startTime = now;
            % if it is continue, read the data
            if nargin == 3 && strcmp(sp1,'continue') 
                if obj.sPara == 1
                    subfile = ['data' filesep obj.sName '.mat'];
                    if ~exist(subfile)
                        error('No such data file');
                    else
                        load(subfile); % trials, expInfo, seq, resp
                        %restore the data
                        obj.seq = seq;
                        obj.resp = resp;
                        obj.aData = aData;
                        obj.curTrl = length(obj.resp)+1;
                        clear trials expInfo aData seq resp; 
                    end
                else % new subject, check if there is already a file
                    subfile = ['data' filesep obj.sName '.mat'];
                    if exist(subfile)
                        error('Data file already exist');
                    end
                end
            end
        end
        
         function saveData(obj,filename)
            if nargin<2
                filename = obj.sName;
            end
            obj.endTime = now;
            
            %create a directory for storing files
            datapath = [pwd filesep 'data'];
            if exist(datapath, 'dir') ~= 7
                mkdir(datapath)
            end
            
            seq = obj.seq;
            resp = obj.resp;
            trials = [seq(1:length(obj.resp),:), resp];
            expInfo.name = obj.sName;
            expInfo.age = obj.sAge;
            expInfo.sex = obj.sGender;
            expInfo.para = obj.sPara;
            expInfo.time = [obj.startTime, obj.endTime];
            aData = obj.aData;
            u = obj.u;
            
            id = 0;
            while(exist([datapath filesep filename '.mat']) == 2)
                filename = [filename num2str(id)];
                id = id + 1;
            end
            
            save([datapath filesep filename '.mat'],'trials','expInfo','aData','seq','resp','u');
            csvwrite([datapath filesep filename '.txt'],trials);
            
        end
        
        
    end
end
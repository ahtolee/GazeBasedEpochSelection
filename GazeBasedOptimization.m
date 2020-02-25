function [Epochs_optimized] = GazeBasedOptimization(EpochN_select,EpochTime,Reversals,blockind,ValidityLeftEye,ValidityRightEye,Artefacts,Fs,FsGaze) 

% Function GazeBasedOptimization.m
% Used in steady-state VEP analysis to perform optimal epoch selection based on the eye tracking data (validity codes) collected during the steady-state
% stimulus presentation. The output of the algorithm is the epoch division that maximized overall gaze quality for the given number of epochs.
% Algorithm also takes into account predefined artefact events, as those segements are excluded from the epochs.
%
% Used in the study:
% "Use of complex visual stimuli allows controlled recruitment of cortical networks in infants"
% Authors: Eero Ahtola1,2, Susanna Stjerna1, Anton Tokariev1,3, Sampsa Vanhatalo1,3
%
% Eero Ahtola. 2020.
% eero.ahtola@hus.fi

% %INPUTS VARIABLES:
% EpochN_select     = Give the number epochs of that will be searched (target, e.g. 60)
% EpochTime         = Duration of one epoch in sec (e.g. 1)
% Reversals         = Matrix array of information of all reversal segments in the EEG data 
%                  1st column -> Start index (sample) of each reversal segment
%                  2nd column -> Duration (in samples) of each reversal segment
%                  3rd column -> The sample where the phase reversal occured in the stimulus. Should be in the middle of the reversal segment.
%                  4th column -> Block number of for which each reversal segment belongs to.
% blockind          = Separator indices (start, end) for each stimulus block (uninterrupted stimulation sequence). In gaze data samples. 
% ValidityLeftEye   = Validity codes from the eye tracker for left eye. Code 0 means that gaze is directed towards the stimulus.
% ValidityRightEye	= Validity codes from the eye tracker for right eye. Code 0 means that gaze is directed towards the stimulus.
% Artefacts         = Predefined reversal segments that are likely to contain EEG artefacts. These segments will be excluded from the final epoch selection.  
% Fs                = Sampling frequency for EEG data (e.g. for samples in Reversal matrix)
% FsGaze            = Sampling frequency for gaze data (for ValidityLeftEye and ValidityRightEye; probably 120Hz)

% %OUTPUT VARIABLES:
% Epochs_optimized = Output matrix with information about the optimal epoch division. Size: EpochN_select x 4.
%                  1st column -> Start index (sample) of each epoch (using given sampling rate of the EEG data)
%                  2nd column -> End index (sample) of each epoch (using given sampling rate of the EEG data)
%                  3rd column -> Block number of for which each epoch belongs to.
%                  4th column -> Gaze Quality index (e.g mean validity between 0-1) for each epoch. Overall Gaze Quality = mean(orind(:,4))                    

% *** Pre-processing. Phase 0.

revGzInd = [];
NBlocks = length(unique(Reversals(:,4)));

%Check equivalency between samples of reversal segments and samples of eye tracker data
for n = 1:NBlocks
    revsInBlocks(n) = sum(Reversals(:,4)==n); %Number of reversal segments in each stimulation block
    StimCyclesInBlocks(n) = revsInBlocks(n)/2;    
    revGzBlock_start = [blockind(n,1):StimCyclesInBlocks(n)*FsGaze/revsInBlocks(n):blockind(n,2)]';
    revGzBlock_end = revGzBlock_start+StimCyclesInBlocks(n)*FsGaze/revsInBlocks(n)-1;    
    revGzInd = vertcat(revGzInd,[revGzBlock_start,revGzBlock_end]); %Start and end of each reversal segment in gaze data samples
end

%Check eye tracker validity codes during each reversal
revN=length(revGzInd);
revValid = zeros(revN,1);
for i=1:revN
    ValidL=ValidityLeftEye(revGzInd(i,1):revGzInd(i,2));
    ValidR=ValidityRightEye(revGzInd(i,1):revGzInd(i,2));   
    revValid(i)=sum(ValidL<0.5 & ValidR<0.5)/length(ValidL); %Eye tracker validity codes must be 0 for both eyes
end

NrevForEpochs = round(EpochTime/(mean(Reversals(:,2))/Fs)); %Number of adjacent reversal segements that will form one each epoch (integer)


% ***Optimal epoch search: Phase 1. Apply gaze threshold. 

% Set validity of the reversals with EEG artefacts to zero (when Artefacts=1) 
revValid_modified = revValid;
revValid_modified(Artefacts==1) = 0;

GazeThresholdStep = 0.1;            %Gaze tolerance step size for the search algorithm
GazeThreshold = 1-GazeThresholdStep;          
NfoundEpochs = 0;

% Check number of valid epochs while lowering the Gaze tolerance step-by-step
% Stop when appropriate number of epochs (or more) have been found
while NfoundEpochs<EpochN_select  
    
    tolValid = revValid_modified>=GazeThreshold;     
    tolValid = tolValid>0 & Artefacts<1;   
    
    A=[]; revValidBlocks=[]; ArtefactsBlock=[];
    
    % Matrix A size is NBlocks x reversal segments in each block
    %In each column 1=gaze validity in that segement exceeds the optimization threshold (GazeThreshold) 
    for m = 1:NBlocks               
        Acol = tolValid(Reversals(:,4)==m);        
        A(:,m)=vertcat(Acol,NaN*ones(max(revsInBlocks)-length(Acol),1)); %pad with NaNs if various block sizes are used in the same recording
        
        revValidBlocks_col=revValid_modified(Reversals(:,4)==m);
        revValidBlocks(:,m)=vertcat(revValidBlocks_col,NaN*ones(max(revsInBlocks)-length(revValidBlocks_col),1)); %pad with NaNs if vaious Blocksizes
        
        ArtefactsBlockcol=Artefacts(Reversals(:,4)==m);
        ArtefactsBlock(:,m)=vertcat(ArtefactsBlockcol,NaN*ones(max(revsInBlocks)-length(ArtefactsBlockcol),1));  %pad with NaNs if vaious Blocksizes  
    end    
    
    %Use the help function to locate adjacent valid reversal segments
    B = check_gaze_valid_neighbours(A);  
    
    % In matrix B
    % B(:,1) -> Starting indices for a pattern of adjacent valid reversal segments. Counted from start of a block. 
    % B(:,2) -> Corresponding number of adjacent valid reversal segments  
    % B(:,3) -> Block number of the detection
        
    if ~isempty(B)        
        B(:,4) = floor(B(:,2)/NrevForEpochs); % B(:,4) -> Number of full adjacent epochs (each consisting of multiple reversal segments)
        B(:,5) = B(:,2)-B(:,4)*NrevForEpochs; % B(:,5) -> Remainder number of reversal segments after the selected full epochs

        NfoundEpochs=sum(B(:,4));             %Number of acceptable valid epochs
    else
        NfoundEpochs=0;                       %No epochs found
    end
    
    % Compare the number of detected acceptable epochs to goal criterion
    % Lower optimization threshold if neccessary
    if NfoundEpochs<EpochN_select
        GazeThreshold = GazeThreshold-GazeThresholdStep; 
    end 
    
    % Optimization threshold negative, not enough valid data, stop the algorithm.
    if GazeThreshold<-0.2
        disp(['NOT ENOUGH VALID DATA. REJECT RECORDING.']); 
        break; 
    end    
end

%Pick the accepted reversal segments and form epochs
ep_start_optim=[]; ep_end_optim=[]; ep_validity_optim=[]; ep_block_optim=[];

%hits -> rev.segments with valid full epochs
hits=find(B(:,4)>0);

% *** Optimal epoch search: Phase 2. Optimize epoch division.
%
% If valid consecutve segment is longer than full epoch (i.e. remainder>0
% in B matrix), we try different combinations for the epoch division within
% that segement, and select the one that yields highest overall validity.

for n=1:length(hits)
    
    % Check one detection candidate at a time and select corresponding rev.segment gaze validity values
    candidates = B(hits(n),:);    
    candidateValid = revValidBlocks(candidates(1):candidates(1)+candidates(2)-1,candidates(3));  
    
    % Check if there were remainder segments after the epochs division?
    % if so, look for the best combination how to make the division (taking into an account the remainder)
    if candidates(5)>0      
        
        % Create all combinations of rev.segments that could be left out of the selected epoch
        combos = combnk([1:length(candidateValid)],candidates(5)); 
        % Sort them by their gaze validity (ascending order)
        [~,I]=sort(mean(candidateValid(combos),2));    
        
        % Start search from the epoch division candidates that yield highest overall validity
        combo_allowed=0; m=1;        
        while combo_allowed<1

            % Try leaving out rev.segments with lowest gaze validity
            tried_combo=setxor([1:candidates(2)],combos(I(m),:));       
            
            % Check if this selection is valid, e.g. it forms full epochs
            % from adjacent rev.segments.
            check_combo_adj = [];
            for x=1:candidates(4)             
                check_combo_adj(x) = sum(diff(tried_combo((x-1)*NrevForEpochs+1:x*NrevForEpochs))==1); %Number of adjacent rev.segments in each epoch
            end                         
            if sum(check_combo_adj==NrevForEpochs-1)==candidates(4) %Check if all full epochs were possible to create (no gaps)                             
                
                combo_allowed=1;    %This handle terminates the search when optimal solution is found  
                                
                %Store search resuls -> epochs within the best allowed candidate pattern       
                
                revsInBlock = Reversals(Reversals(:,4)==candidates(3),:);  % select rev.segments from one block
                for k = 1:candidates(4)                        % handle epochs in the pattern one by one
                    
                    % select rev.segments from one epoch                    
                    oneEpoch = revsInBlock(candidates(1)+tried_combo((k-1)*NrevForEpochs+1)-1:candidates(1)+tried_combo((k-1)*NrevForEpochs+1)-1+NrevForEpochs-1,:);                     
                    ep_start_optim(end+1) = oneEpoch(1,1);                                                         %Start ind of an epoch
                    ep_end_optim(end+1) = oneEpoch(NrevForEpochs,1)+oneEpoch(NrevForEpochs,2)-1;                   %End ind of an epoch
                    candidateValidSelected = candidateValid(tried_combo);
                    ep_validity_optim(end+1) = mean(candidateValidSelected((k-1)*NrevForEpochs+1:k*NrevForEpochs));   %Avg validity for an epoch
                    ep_block_optim(end+1) = candidates(3);                                                         %Block number for an epoch
                end               
            else
                m=m+1;
            end
        end    
    else
        
        %Store epochs of the candidate pattern (cases when remainder segments do not exist)          
        revsInBlock = Reversals(Reversals(:,4)==candidates(3),:);  % select rev.segments from one block
        for k = 1:candidates(4)        
            oneEpoch = revsInBlock(candidates(1)+(k-1)*NrevForEpochs:candidates(1)+k*NrevForEpochs-1,:);  % select rev.segments from one epoch
            ep_start_optim(end+1) = oneEpoch(1,1);
            ep_end_optim(end+1) = oneEpoch(NrevForEpochs,1)+oneEpoch(NrevForEpochs,2)-1;
            ep_validity_optim(end+1) = mean(candidateValid((k-1)*NrevForEpochs+1:k*NrevForEpochs));
            ep_block_optim(end+1) = candidates(3);
        end
    end
end
ep_start_optim=ep_start_optim'; ep_end_optim=ep_end_optim'; ep_validity_optim=ep_validity_optim'; ep_block_optim=ep_block_optim';  


% If the first phase of the algorithm produced more epochs than was needed...
if NfoundEpochs>EpochN_select
    [~,I] = sort(ep_validity_optim,'descend');                    %Sort foundEpochs by average gaze validity     
    ep_start_optim = ep_start_optim(sort(I(1:EpochN_select)));
    ep_end_optim = ep_end_optim(sort(I(1:EpochN_select)));
    ep_validity_optim = ep_validity_optim(sort(I(1:EpochN_select)));
    ep_block_optim = ep_block_optim(sort(I(1:EpochN_select)));    
end
   

%Print the results info:
disp('*** Results from the epoch selection optimization ***');
disp(['Epochs found in first phase: ' num2str(NfoundEpochs) ', with gaze threshold of: ' num2str(100*GazeThreshold) '%']);
disp(['Epochs selected in second phase: ' num2str(length(ep_validity_optim))]);
disp(['Average Gaze quality in the selected epochs: ' num2str(100*mean(ep_validity_optim)) '%']);


% Function output variable
Epochs_optimized(:,1) = ep_start_optim;
Epochs_optimized(:,2) = ep_end_optim;
Epochs_optimized(:,3) = ep_block_optim;
Epochs_optimized(:,4) = ep_validity_optim;


% %Force all epochs to even length (defined by EpochLength)
% ep_extra=round((ep_end_optim-ep_start_optim-EpochLength)/2);
% Epochs_optimized=[];
% if ~isempty(ep_start_optim)
%     Epochs_optimized(:,1)=ep_start_optim+ep_extra;
%     Epochs_optimized(:,2)=Epochs_optimized(:,1)+EpochLength-1;
%     Epochs_optimized(:,3)=ep_block_optim;
%     Epochs_optimized(:,4)=ep_validity_optim;
% else
%     Epochs_optimized=NaN*ones(1,4);
% end


% *** Optional graph output:
figure; set(gcf,'position',[100 100 1450 630]); hold on; set(gcf,'Color','w'); box off;

% A) All rev. segments and their avg gaze validity
graph_x = [Reversals(1,1): Reversals(end,1)+Reversals(end,2)];
graph_y1 = zeros(size(graph_x));
for n=1:length(revValid)
    graph_y1(Reversals(n,1)-Reversals(1,1)+1:Reversals(n,1)+Reversals(n,2)-Reversals(1,1)) = revValid(n);   
end
graph_x1A = Reversals(Artefacts,3); %mark artefacts
graph_y1A = ones(size(graph_x1A))*1.05;
graph_x1B = (Reversals(find(diff(Reversals(:,4))>0)+1,1) + Reversals(find(diff(Reversals(:,4))>0)+1,1))/2; %mark stimulation blocks
graph_y1B = [0,105];

subplot(311); hold on;
%line([0 graph_x(end)-graph_x(1)], 100*GazeThreshold*[1 1],'Color',[0.2 0.2 0.2],'LineStyle','--','LineWidth',1.5);
line([0 graph_x(end)-graph_x(1)], 100*mean(revValid)*[1 1],'Color',[0.46 0.67 0.19],'LineWidth',1.5);
area(graph_x-graph_x(1),100*graph_y1,'FaceColor',[0 0.55 0.85],'EdgeColor',[0 0.40 0.70],'LineWidth',0.1);
for n=1:length(graph_x1B), line(graph_x1B(n)*[1 1]-graph_x(1), [0,110],'Color',[0.7 0.7 0.7],'LineStyle',':','LineWidth',2); end
plot(graph_x1A-graph_x(1),100*graph_y1A,'r.','MarkerSize',5);
axis([0 graph_x(end)-graph_x(1) 0 105]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.005 .005] , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 1         );
set(gca,'FontName','Helvetica','FontSize',14);
set(gca,'XTickLabels',{});
set(gca,'layer','top');   

% B) Rev. segments above the optimization threshold and their avg gaze validity
graph_y2 = graph_y1;
graph_y2(graph_y2<GazeThreshold)=0;
subplot(312); hold on;
line([0 graph_x(end)-graph_x(1)], 100*GazeThreshold*[1 1],'Color',[0.2 0.2 0.2],'LineStyle','--','LineWidth',1.5);
line([0 graph_x(end)-graph_x(1)], 100*mean(revValid(revValid>GazeThreshold))*[1 1],'Color',[0.46 0.67 0.19],'LineWidth',1.5);
area(graph_x-graph_x(1),100*graph_y2,'FaceColor',[0 0.55 0.85],'EdgeColor',[0 0.40 0.70],'LineWidth',0.1);
for n=1:length(graph_x1B), line(graph_x1B(n)*[1 1]-graph_x(1), [0,110],'Color',[0.7 0.7 0.7],'LineStyle',':','LineWidth',2); end
axis([0 graph_x(end)-graph_x(1) 0 105]);
ylabel('Gaze quality (%)');
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.005 .005] , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 1         );
set(gca,'FontName','Helvetica','FontSize',14);
set(gca,'XTickLabels',{});
set(gca,'layer','top');   

% C) Epoch selected by the optimization search and their avg gaze validity
graph_y3 = zeros(size(graph_x));
for n=1:length(ep_validity_optim)
    graph_y3(Epochs_optimized(n,1)-Reversals(1,1)+1:Epochs_optimized(n,2)-Reversals(1,1)) = ep_validity_optim(n);   
end

subplot(313); hold on;
line([0 graph_x(end)-graph_x(1)], 100*GazeThreshold*[1 1],'Color',[0.2 0.2 0.2],'LineStyle','--','LineWidth',1.5);
line([0 graph_x(end)-graph_x(1)], 100*mean(ep_validity_optim)*[1 1],'Color',[0.46 0.67 0.19],'LineWidth',1.5);
area(graph_x-graph_x(1),100*graph_y3,'FaceColor',[0 0.55 0.85],'EdgeColor',[0 0.40 0.70],'LineWidth',0.1);
for n=1:length(graph_x1B), line(graph_x1B(n)*[1 1]-graph_x(1), [0,110],'Color',[0.7 0.7 0.7],'LineStyle',':','LineWidth',2); end
axis([0 graph_x(end)-graph_x(1) 0 105]);
xlabel('Samples');
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.005 .005] , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 1         );
set(gca,'FontName','Helvetica','FontSize',14);
set(gca,'layer','top');   

end

function B=check_gaze_valid_neighbours(A)

% Help function for gaze optimization in eye tracker assisted SSVEP epoch selection algorithm. 
% Eero Ahtola. 2020.
% eero.ahtola@hus.fi

%-Input-
% A=binary matrix

%-Output-
% Adjacency matrix:
% B(x,y,z)= detection of length y at index x of z:th column

[~,N] = size(A);
B = [];
for n = 1:N
    a = A(:,n);
    ind = find(a>0);
    m = 1;
    while m<length(ind)+1
       range = 0;
       foundsize = 0;
       hasneighbour = 1;
       while hasneighbour==1
           if ismember(ind(m)+range,ind)
              range = range+1;
           else               
              foundsize = foundsize+range;
              hasneighbour = 0;
           end           
       end       
       B(end+1,:) = [ind(m),foundsize,n];
       m = m+foundsize;
    end    
end

end




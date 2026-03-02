close all; clearvars; clc

% load files
addpath(genpath('\\XXX\HierarchicalPriors\code\behavior'));
indir = '\\XXX\HierarchicalPriors\data\MEG\PTBoutput';

% subject loop
subfilz = dir(fullfile(indir,'*block1*'));

data = struct;
RTpersub = zeros(length(subfilz),1);
ACCpersub = zeros(length(subfilz),1);
RTpercon = zeros(length(subfilz),2,3);
ACCpercon = zeros(length(subfilz),2,3);
for isub = 1:length(subfilz)
    
    runfilz = dir(fullfile(indir,['HierarchicalPriors_subject' subfilz(isub).name(end-12:end-11) '_block*']));
    
    for irun = 1:length(runfilz)
    
        load(fullfile(indir,runfilz(irun).name));
    
        if irun == 1
        	data = output;
        else
            data = [data output];
        end
        
    end
    
    catchtrials = logical(extractfield(data,'catch'));
    data = data(catchtrials);
    
    % sanity check
    if length(data) ~= length(runfilz)*24
        error('Not right catch trial amount...');
    end
    
    RT = extractfield(data,'RT');
    response = extractfield(data,'correct_response');
    
    % Identify no-response trials

    % For fixation cross condition a 'no response' is the only way of
    % getting an error, because invalidly fast responses to a fixation
    % cross are very unlikely as this is a random and unexpected event.
    % Therefore, no-responses are directly reflected in the error rate for
    % this condition). We only need to remove these trials from the further
    % RT analysis (we need them for accuracy). no-response is defined as
    % a negative value, or a value larger than 7.6 sec max (5 sec video + 
    % 3 sec after end of the video to respond minus first possible catch 
    % trial (a fixation) at 400 msec into trial.
    noresponseFIX = ~extractfield(data,'catchFixOrOcclude') & (RT < 0 | RT > 7.6);
    RT(noresponseFIX) = NaN;

    % For occlusion condition this is any response faster than 600 msec
    % (i.e., 400 msec occlusion + 200 msec minimum time to make response).
    % We want to remove those from further RT AND accuracy analysis, plus
    % store the percentage per subject
    noresponseOCC = extractfield(data,'catchFixOrOcclude') & (RT < 0.6 | RT > 7.6);
    RT(noresponseOCC) = NaN;
    response(noresponseOCC) = NaN;
    occRT2fast(isub) = sum(noresponseOCC)/sum(extractfield(data,'catchFixOrOcclude'));

    
    RTpersub(isub) = nanmean(RT);
    ACCpersub(isub) = nanmean(response);
    
    fixVSocc = extractfield(data,'catchFixOrOcclude')' == [0 1];

    % video condition: 1 = normal, 2 = upside down, 3 = scrambled
    conID = extractfield(data,'condition')' == [1 2 3];
    
    % Now determine performance per catch trial type
    for icon = 1:size(conID,2)
        for itype = 1:2
            ID = logical(conID(:,icon) + fixVSocc(:,itype) == 2);% select video condition and catch trial type
        
            ACCpercon(isub,itype,icon) = nanmean(response(ID));

            % for RT only use correct trials
            ID = logical(ID+response' == 2);

            RTpercon(isub,itype,icon) = nanmean(RT(ID));
            
    
        end
    end

end

% subjects 2 exclude
% based on no-responses to occlusion task:
subject2exclude_occ2fast = occRT2fast > 0.7;
% based on average accuracy in occlusion task
temp = squeeze(mean(ACCpercon(:,2,:),3));
subject2exclude_lowacc = temp < mean(temp)-2.5*std(temp);
subject2exclude = logical(subject2exclude_lowacc + subject2exclude_occ2fast');

% remove subjects 2 exclude
ACCpercon(subject2exclude,:,:) = [];
RTpercon(subject2exclude,:,:) = [];

% RT on the occlusion task should be counted from occlusion offset instead of onset, so subtract 400 msec for the occlusion
% Similarly, RT on fixation task is now relative to onset of color change, but that's correct because from that moment onwards it makes sense to respond
RTpercon(:,2,:) = RTpercon(:,2,:) - .4;

% compute balanced integration score (BIS)
Zpc = (ACCpercon - mean(ACCpercon(:))) ./ std(ACCpercon(:));
Zrt = (RTpercon - mean(RTpercon(:))) ./ std(RTpercon(:));
BIS = Zpc - Zrt;

% reshape data for stats in JASP
data4JASP = [reshape(ACCpercon,40,[]) reshape(RTpercon,40,[]) reshape(BIS,40,[])];

% store data for correlation plots with dRSA peak results
save(fullfile(indir,'behavioural_results'),'ACCpercon','RTpercon','BIS');

%% Figures behavioural performance
% per condition
fs = 8;

close(figure(1));
figure
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [5 5 11 6]);

%accuracy
subplot(2,2,1)
hold on
plotSpread(squeeze(ACCpercon(:,1,:)),'spreadWidth',1,'distributionColors',[0 0 0 ; 0 0 0 ; 0 0 0],'distributionMarkers','.','showMM',2);
hold off
set(gca,'ylim',[0.25 1]);
set(gca,'xlim',[0.5 3.5]);
set(gca,'XTick',[1 2 3])
set(gca,'XTickLabel','')
set(gca,'YTick',[0 0.25 0.5 0.75 1])
set(gca,'YTickLabel',{'0', '25', '50','75', '100'})
set(gca, 'FontSize', fs,'Fontname','Arial');
ylabel('Performance [%]','FontSize',fs);
title('fixation cross task')

subplot(2,2,2)
hold on
plotSpread(squeeze(ACCpercon(:,2,:)),'spreadWidth',1,'distributionColors',[0 0 0 ; 0 0 0 ; 0 0 0],'distributionMarkers','.','showMM',2);
hold off
set(gca,'ylim',[0.25 1]);
set(gca,'xlim',[0.5 3.5]);
set(gca,'XTick',[1 2 3])
set(gca,'XTickLabel','');
set(gca,'YTick',[0 0.25 0.5 0.75 1])
set(gca,'YTickLabel',{'0', '25', '50','75', '100'})
set(gca, 'FontSize', fs,'Fontname','Arial');
title('occlusion task')

%RT
subplot(2,2,3)
hold on
plotSpread(squeeze(RTpercon(:,1,:)),'spreadWidth',1,'distributionColors',[0 0 0 ; 0 0 0 ; 0 0 0],'distributionMarkers','.','showMM',2);
hold off
set(gca,'ylim',[0.42 1.51]);
set(gca,'xlim',[0.5 3.5]);
set(gca,'XTick',[1 2 3])
set(gca,'XTickLabel',{'normal','invert','scram'});
set(gca,'YTick',[.5 1 1.5 2])
set(gca,'YTickLabel',{'0.5', '1.0','1.5', '2.0'})
set(gca, 'FontSize', fs,'Fontname','Arial');
ylabel('RT [sec]','FontSize',fs);

subplot(2,2,4)
hold on
plotSpread(squeeze(RTpercon(:,2,:)),'spreadWidth',1,'distributionColors',[0 0 0 ; 0 0 0 ; 0 0 0],'distributionMarkers','.','showMM',2);
hold off
set(gca,'ylim',[0.42 1.51]);
set(gca,'xlim',[0.5 3.5]);
set(gca,'XTick',[1 2 3])
set(gca,'XTickLabel',{'normal','invert','scram'});
set(gca,'YTick',[.5 1 1.5 2])
set(gca,'YTickLabel',{'0.5', '1.0','1.5', '2.0'})
set(gca, 'FontSize', fs,'Fontname','Arial');
% xlabel('catch trial type','FontSize',8);

cd('XXX\HierarchicalPriors\figures');
% % set(gcf,'renderer','Painters')
% print -depsc -tiff -r600 -painters HierarchicalPriors_FigureX_behavior.eps


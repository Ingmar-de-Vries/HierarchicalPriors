function HierarchicalPriors_DynamicModelRDMs_eyeTracker(cfg,eyeSub)

if isfield(cfg,'cluster')
    indirEye = '//XXX/HierarchicalPriors/data/MEG/ELpreprocessed';
else
    indirEye = '\\XXX\HierarchicalPriors\data\MEG\ELpreprocessed';
end

% Euclidean distance using matlab's norm function
fn= sprintf('%s%cSUB%02d',indirEye, filesep, eyeSub);
load(fn,'data');

nstim = 14;
eyePos = data;
clear data

% only select eye positions
eyePos.trial = eyePos.trial(:,[1 2 4 5],:);

% select video time
toi = dsearchn(eyePos.time',[0 5]');
toi = toi(1):toi(end);
eyePos.trial = eyePos.trial(:,:,toi);
eyePos.time = eyePos.time(toi);

% interpolate missing data due to blinks or reject trial if gap was too large (e.g., if eye tracker really lost eye for a while)
temp = eyePos.trial;
ntrials = size(temp,1);
ntime = size(temp,3);
temp = reshape(permute(temp,[2 3 1]),size(temp,2),ntrials*ntime);
missingdata = find(logical(sum(temp==0 + isnan(temp))));
gapends = [find(diff(missingdata)~=1) length(missingdata)];
gapstarts = [1 gapends(1:end-1)+1];
padding = 5;% in samples
paddedgaps = [];
for igap = 1:length(gapstarts)
    paddedgaps = [paddedgaps missingdata(gapstarts(igap))-padding:missingdata(gapends(igap))+padding];
end% gap loop
paddedgaps(paddedgaps<1)=[];
idx = 1:size(temp,2);
idx2interp = idx;
idx2interp(paddedgaps) = [];
temp(:,paddedgaps) = [];
temp = interp1(idx2interp,temp',idx,'pchip','extrap')';
eyePos.trial = permute(reshape(temp,size(temp,1),ntime,ntrials),[3 1 2]);

% average over repetitions of same video, per condition
temp = zeros(nstim,3,4,length(toi));
% eyeMov = temp;
for icon = 1:3
    for istim = 1:nstim
        currentstim = istim+(icon-1)*20;
        temp(istim,icon,:,:) = squeeze(mean(eyePos.trial(eyePos.trialinfo == currentstim,:,:)));
        %         eyeMov(istim,icon,:,:) = interp1(t4eyeMov,diff(squeeze(temp(istim,icon,:,:)),[],2)',eyePos.time,'pchip','extrap')';
    end
end
eyePos = temp;

RDMeyeNORMAL = zeros(nstim,nstim,size(eyePos,4),size(eyePos,4));
RDMeyeINVERT = zeros(nstim,nstim,size(eyePos,4),size(eyePos,4));
RDMeyeSCRAM = zeros(nstim,nstim,size(eyePos,4),size(eyePos,4));
for istim1 = 1:nstim
    for istim2 = 1:nstim
        for itime1 = 1:length(toi)
            for itime2 = 1:length(toi)
                RDMeyeNORMAL(istim1,istim2,itime1,itime2) = norm(squeeze(eyePos(istim1,1,:,itime1))-squeeze(eyePos(istim2,1,:,itime2)));
                RDMeyeINVERT(istim1,istim2,itime1,itime2) = norm(squeeze(eyePos(istim1,2,:,itime1))-squeeze(eyePos(istim2,2,:,itime2)));
                RDMeyeSCRAM(istim1,istim2,itime1,itime2) = norm(squeeze(eyePos(istim1,3,:,itime1))-squeeze(eyePos(istim2,3,:,itime2)));
            end
        end% time 1 loop
    end% stim 2 loop
end% stim 1 loop

% save correlation matrices
save([fn '_dynRDM_eyeTracker'],'RDMeyeNORMAL','RDMeyeINVERT','RDMeyeSCRAM');

end

function HierarchicalPriors_checkAtlases(cfg)

addpath('\\XXX\HierarchicalPriors\code\dRSA');
cfg.path = '\\XXX\HierarchicalPriors'; 

for iSub = cfg.SubVec%[1:11 13:21]%
    
    % source data
    indirSource = fullfile(cfg.path,'data','MEG','brainstorm_database','HierarchicalPriors','data',cfg.subjectKEY{iSub,2},cfg.subjectKEY{iSub,2});
    indirAtlas = fullfile(cfg.path,'data','MEG','brainstorm_database','HierarchicalPriors','anat',cfg.subjectKEY{iSub,2});
    
    % atlases
    sourcefile = dir(fullfile(indirSource,'*KERNEL*'));
    fn2load = sprintf('%s%c%s',indirSource,filesep,sourcefile.name);
    kernel = load(fn2load);
    cortex = kernel.SurfaceFile(6:end-4);
    
    atlasfile = dir(fullfile(indirAtlas,'*cortex*'));
    atlasID = contains(extractfield(atlasfile,'name'),cortex);
    fn2load = sprintf('%s%c%s',indirAtlas,filesep,atlasfile(atlasID).name);
    
    % check if the here loaded source file used the same atlas as the here
    % loaded atlas file:
    if ~contains(fn2load,cortex)
        error(['cortical surface used for source reconstruction does not match the one in atlas for subject ' num2str(iSub)]);
    end
    
    atlas = load(fn2load);
    atlasIdx = logical(contains(extractfield(atlas.Atlas,'Name'),cfg.atlas));
    atlas = atlas.Atlas(atlasIdx).Scouts;
    
    idx2remove = contains(extractfield(atlas,'Label'),'Background');
    atlas(idx2remove) = [];
    
    % set atlas of subject 1 as reference atlas and compare all others
    % against it to make sure we're using the same atlas for all subjects
    if iSub == cfg.SubVec(1)
        referenceAtlas = atlas;
    else
        if ~all(strcmp(extractfield(referenceAtlas,'Label'),extractfield(atlas,'Label')))
            error(['The atlas of subject ' num2str(iSub) ' does not match the reference atlas'])
        end
    end
    
end

function Unpredict_cluster_shell(parms)
% determine whether to send job to cluster, or run script locally

matlab_version = sscanf(version('-release'),'%d%s');
if matlab_version(1) ~= 2020
    error('The package for sending jobs to the cluster only works on matlab 2020, please use that');
end

cluster = parcluster;
job = createJob(cluster);

% determine what to parallelise over
subs = 1;
rois = 1;
cons = 1;
pthresh = 1;
sim = 1;

if any(contains(parms.parallel,'sub'))
    subs = parms.subjects;
end
if any(contains(parms.parallel,'roi'))
    rois = parms.ROI;
end
if any(contains(parms.parallel,'con'))
    cons = parms.conditions;
end
if any(contains(parms.parallel,'pthresh'))
    pthresh = parms.pthresh;
end
if any(contains(parms.parallel,'sim'))
    sim = parms.models2sim;
end
    
% create tasks
for isub = subs
    for iroi = rois      
        for icon = cons
            for ithresh = pthresh
                for isim = sim
                    createTask(job,str2func(parms.script2run), 0,{parms, isub, iroi, icon, ithresh, isim});%
                end
            end
        end
    end
end

% submit job
submit(job);

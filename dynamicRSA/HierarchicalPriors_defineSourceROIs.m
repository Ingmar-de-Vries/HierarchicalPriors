function HierarchicalPriors_defineSourceROIs()

% visual areas
ROIdefinition.parcels{1} = {'V1'}; ROIdefinition.names{1} = {'V1'}; 
ROIdefinition.parcels{2} = {'V2'}; ROIdefinition.names{2} = {'V2'};
ROIdefinition.parcels{3} = {'V3','V4'}; ROIdefinition.names{3} = {'V3V4'};% combined to keep #vertices similar across ROIs and because MEG can't clearly distinguish between V3 and V4 anyway 

% action observation network (AON)
% ROIdefinition.parcels{5} = {'STV', 'TPOJ1','TPOJ2', 'PHT', 'TE1p'}; ROIdefinition.names{5} = {'aLOTC'};% anterior Lateral Occipital Temporal Cortex
ROIdefinition.parcels{4} = {'V4t', 'FST', 'MT', 'MST', 'LO1', 'LO2', 'LO3', 'PH', 'PHT', 'TPOJ2', 'TPOJ3'}; ROIdefinition.names{4} = {'LOTC'};%posterior Lateral Occipital Temporal Cortex
ROIdefinition.parcels{5} = {'PF', 'PFt', 'AIP', 'IP2'}; ROIdefinition.names{5} = {'aIPL'};% anterior Inferior Parietal Lobule
ROIdefinition.parcels{6} = {'IFJa', 'IFJp', '6r', '6v', 'PEF', 'IFSp', '44', '45'}; ROIdefinition.names{6} = {'VentPremot'};%

% dorsal stream / dorsal attention network
ROIdefinition.parcels{7} = {'V3A', 'V3B', 'V6', 'V6A', 'V7', 'IPS1', 'IP0', 'DVT'}; ROIdefinition.names{7} = {'PostDorsStream'};%
ROIdefinition.parcels{8} = {'MIP', 'VIP', 'LIPv', 'LIPd', 'IP1', 'AIP', '7PC'}; ROIdefinition.names{8} = {'AntDorsStream'};% 
ROIdefinition.parcels{9} = {'6d', '6a', 'FEF', '55b', '6ma', '6mp'}; ROIdefinition.names{9} = {'DorsPremot'};%
ROIdefinition.parcels{10} = {'s6-8', 'i6-8', 'SFL', '8Ad', '8Av', '8BL', '8C'}; ROIdefinition.names{10} = {'pDLPFC'};% posterior Dorso-Lateral Pre-Frontal Cortex

save('\\XXX\HierarchicalPriors\code\dRSA\ROIdefinitions','ROIdefinition');

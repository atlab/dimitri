%{
aodi.Joseph (computed) # dumping data for Joseph
-> vis2p.MaskGroup
-> vis2p.MaskGroupRaw
-----
%}

classdef Joseph < dj.Relvar & dj.AutoPopulate
    
    properties(Constant)
        pmtOffset = -495424
        savePath = '~/Google Drive/Neurophotonics/aod-data'
    end
    
    properties
        popRel  = vis2p.MaskGroup * pro(vis2p.MaskGroupRaw) & vis2p.VisStims & aodi.TrialTraces & vis2p.MaskTracesRaw
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            disp 'fetching raw traces...'
            [tuple.cellnums, tuple.celltypes, x, y, z, traces] = ...
                fetchn(pro(vis2p.MaskTraces,'mask_type')*pro(vis2p.MaskCells,'img_x','img_y','img_z')*pro(vis2p.MaskTracesRaw,'calcium_trace') & key, ...
                'masknum', 'mask_type', 'img_x', 'img_y', 'img_z', 'calcium_trace');
            times = fetch1(vis2p.MaskGroupRaw & key, 'frame_timestamps');
            traces = double([traces{:}]) - aodi.Joseph.pmtOffset;
            folder = fullfile(aodi.Joseph.savePath,sprintf('%s-%02d',key.exp_date,key.scan_idx));
            disp saving..
            mkdir(folder)
            csvwrite(fullfile(folder,'traces.csv'),traces)
            csvwrite(fullfile(folder,'times.csv'),uint64(times))
            csvwrite(fullfile(folder,'xyz.csv'),[x y z])
            self.insert(key)
        end
    end
    
end
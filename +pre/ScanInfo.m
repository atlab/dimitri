%{
pre.ScanInfo (imported) # header information
-> rf.Scan
---
nframes_requested           : int                           # number of valumes (from header)
px_width                    : smallint                      # pixels per line
px_height                   : smallint                      # lines per frame
um_width                    : float                         # width in microns
um_height                   : float                         # height in microns
bidirectional               : tinyint                       # 1=bidirectional scanning
fps                         : float                         # (Hz) frames per second
zoom                        : decimal(4,1)                  # zoom factor
dwell_time                  : float                         # (us) microseconds per pixel per frame
nchannels                   : tinyint                       # number of recorded channels
nslices                     : tinyint                       # number of slices
slice_pitch                 : float                         # (um) distance between slices
nframe_avg                  : smallint                      # number of averaged frames
fill_fraction=0             : float                         # raster scan fill fraction (see scanimage)
%}

classdef ScanInfo < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = rf.Scan
    end
    
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            reader = pre.getReader(key);
            
            key.nframe_avg = reader.header.acqNumAveragedFrames;
            assert(strcmp(reader.header.fastZImageType, 'XY-Z'), ...
                'we assume XY-Z scanning')
            if reader.header.fastZActive
                key.nframes_requested = reader.header.fastZNumVolumes;
            else
                key.nframes_requested = reader.header.acqNumFrames;
            end
            sz = size(reader);
            key.px_height = sz(2);
            key.px_width  = sz(1);
                       
            %%%% Find appropriate field of view
            fov=560;
            %%%%
            
            zoom = reader.header.scanZoomFactor;
            key.um_height = fov/zoom*reader.header.scanAngleMultiplierSlow;
            key.um_width  = fov/zoom*reader.header.scanAngleMultiplierFast;
            if reader.header.fastZActive
                key.fps =  1/reader.header.fastZPeriod;
                key.slice_pitch = reader.header.stackZStepSize;
            else
                key.fps = reader.header.scanFrameRate;
                key.slice_pitch = 0;
            end
            
            key.bidirectional = ~strncmpi(reader.header.scanMode, 'uni', 3);
            key.zoom = zoom;
            key.dwell_time = reader.header.scanPixelTimeMean*1e6;
            key.nchannels = length(reader.header.channelsSave);
            key.nslices = reader.header.stackNumSlices;
            if key.nslices == 1
                key.fps = reader.header.scanFrameRate;
            end
            key.fill_fraction = reader.header.scanFillFraction;
            
            self.insert(key)
        end
    end
end

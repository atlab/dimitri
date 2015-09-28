%{
carfs.Subfield (computed) # subfield location in x,y,z coordinates
-> carfs.SubfieldSet
subfield  : tinyint   # subfield number.  0=fake subfield
-----
sub_x=0: float   # (degrees) relative to center
sub_y=0: float   # (degrees) relative to center
latency: float   # (ms) latency of subfield's peak
amplitude=0: float  # expressed in sigmas
%}

classdef Subfield < dj.Relvar
        
    methods
        
        function makeTuples(self, key)
            map = fetch1(carfs.NoiseMap & key, 'noise_map');
            sz = size(map);
            
            % compute pixel pitch
            stim = fetchn(vis2p.VisStims & key, 'stim_file');
            monitorSize = stim{1}.params.constants.monitorSize;
            monitorDistance = stim{1}.params.constants.monitorDistance;
            angle = 2*atan(norm(monitorSize)/2/monitorDistance);  % angle across diagonal
            pixelPitch = 180/pi*angle/norm(sz(1:2)/2);  % degrees/pixel
            dt = fetch1(carfs.NoiseMapSet & key, 'bin_size');
            
            sigma = std(map(:));
            r=round(2/pixelPitch);
            f = fspecial('gaussian', 2*r+1, r/2);
            
            thresh = 4.5;  %sigma
            
            subfield = 0;
            while true
                smooth = max(imfilter(abs(map)/sigma,f),[],3);
                amp = max(smooth(:));
                if max(smooth(:)) < thresh
                    break
                end
                
                tuple = key;
                subfield = subfield + 1;
                tuple.subfield = subfield;
                
                % find the peak
                [~,i] = max(abs(map(:)));
                [y, x, z] = ind2sub(sz,i);
                
                % fit gaussian blob
                tuple.sub_x = (x-sz(2)/2) * pixelPitch;
                tuple.sub_y = (y-sz(1)/2) * pixelPitch;
                tuple.latency = z*dt;
                tuple.amplitude = amp;
                self.insert(tuple)
                break
            end
            
        end
    end
    
end
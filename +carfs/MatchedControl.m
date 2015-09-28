%{
carfs.MatchedControl (computed) # nearest tuned non-red cell for each tuned red cell
-> carfs.TraceSet
pair_id  : smallint     # pair id
mask_type : enum('red','neuron')  # role in the pair
-----
-> carfs.Trace
distance   : float   # (um) distsance between the cells
%}

classdef MatchedControl < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = pro(carfs.TraceSet) ...
            & (carfs.Trace*carfs.VonMises & 'mask_type="red"' & 'von_p_value<0.05')
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            cells = pro(carfs.Trace,'mask_type')*carfs.CellXYZ & (carfs.VonMises & 'von_p_value<0.05');
            
            sisters = cells & 'mask_type="red"';
            controls = pro(cells & 'mask_type="neuron"', 'masknum->masknum2', 'x->x2', 'y->y2', 'z->z2');
            
            pairs = pro(sisters*controls, 'pow(x-x2,2)+pow(y-y2,2)+pow(z-z2,2)->d');
            pairs = pairs & key & 'd>400';  % not closer than 10 um to avoid matching duplicate labels
            
            [d, red, neuron] = dj.struct.tabulate(pairs.fetch('d'), 'd', 'masknum', 'masknum2');
            
            iPair = 0;
            % save closest neighbor
            tuple = key;
            [distance, j] = nanmin(d,[],2);
            
            nearestNeuron = arrayfun(@(i) neuron(j(i)), 1:length(j));
            
            for iSister = find(~isnan(distance))'
                iPair = iPair + 1;
                tuple.pair_id = iPair;
                tuple.mask_type = 'red';
                tuple.masknum = red(iSister);
                tuple.distance = sqrt(distance(iSister));
                self.insert(tuple)
                
                tuple.mask_type = 'neuron';
                tuple.masknum = nearestNeuron(iSister);
                self.insert(tuple)
            end
        end
        
    end
end

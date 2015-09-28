%{
carfs.CellXYZ (computed) # position in microns
-> carfs.Trace
-----
x  : float
y  : float
z  : float

%}

classdef CellXYZ < dj.Relvar & dj.AutoPopulate

	properties
		popRel = carfs.TraceSet
	end

	methods(Access=protected)

		function makeTuples(self, key)
            tuples = fetch(carfs.Trace*vis2p.MaskCells & key, ...
                'img_x->x','img_y->y','img_z->z');
			self.insert(tuples)
		end
	end

end
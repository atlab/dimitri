%{
carfs.ImageMapResponses (imported) # responses to images
# add primary key here
-----
# add additional attributes
%}

classdef ImageMapResponses < dj.Relvar & dj.AutoPopulate

	properties
		popRel  % !!! update the populate relation
	end

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			self.insert(key)
		end
	end

end
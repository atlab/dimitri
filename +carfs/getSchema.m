function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'carfs', 'dimitri_carfs');
end
obj = schemaObject;
end

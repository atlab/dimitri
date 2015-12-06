function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    vis2p.getSchema;
    schemaObject = dj.Schema(dj.conn, 'carfs', 'dimitri_carfs');
end
obj = schemaObject;
end

function obj = getSchema
persistent schemaObject

if isempty(schemaObject)
    vis2p.getSchema;
    schemaObject = dj.Schema(dj.conn, 'aodi', 'dimitri_aodi');
end

obj = schemaObject;
end
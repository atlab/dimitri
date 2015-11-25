function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'monet', 'dimitri_monet');
end
obj = schemaObject;
end

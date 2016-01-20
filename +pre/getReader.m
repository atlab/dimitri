function reader = getReader(key, cache_dir)


[path, basename, scanIdx] = fetch1(rf.Session*rf.Scan & key, ...
    'scan_path', 'file_base', 'scan_idx');

if nargin > 1
    slpos =  strfind(path,'/');
    stump = path(slpos(end):end);
    cache_path = [cache_dir stump];
    try
        fprintf('Loading data %s\n', cache_path);
        tpath = getLocalPath(fullfile(cache_path, sprintf('%s_%03u.tif', basename, scanIdx)));
        reader = ne7.scanimage.Reader4(tpath);
    catch
        fprintf('Copying data from %s to %s\n', path, cache_path);
        copyfile(getLocalPath(path), cache_path);
        path = getLocalPath(fullfile(cache_path, sprintf('%s_%03u.tif', basename, scanIdx)));
        reader = ne7.scanimage.Reader4(path);
    end
else
    fprintf('Loading data %s\n', path);    
    path = getLocalPath(fullfile(path, sprintf('%s_%03u.tif', basename, scanIdx)));
    reader = ne7.scanimage.Reader4(path);
end
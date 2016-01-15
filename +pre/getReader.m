function reader = getReader(key)
[path, basename, scanIdx] = fetch1(rf.Session*rf.Scan & key, ...
    'scan_path', 'file_base', 'scan_idx');
path = getLocalPath(fullfile(path, sprintf('%s_%03u.tif', basename, scanIdx)));
reader = ne7.scanimage.Reader4(path);
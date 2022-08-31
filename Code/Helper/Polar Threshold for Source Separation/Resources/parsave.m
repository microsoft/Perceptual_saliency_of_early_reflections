function parsave(fname, data)
    % Written by Tobias Jüterbock
    % Helper function for writing variables to files inside parfor loops
    save(fname, 'data')
end
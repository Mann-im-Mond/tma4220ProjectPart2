function mesh = getCakeMesh(folder, neumannIdentifier)
    addpath('../utils/')
    
    tetrFile = join([folder '/cake_tetr.m'])
    fid=fopen(tetrFile);
    s=textscan(fid,'%f %f %f %f %f %f');
    fclose(fid);
    tetr = [s{1},s{2},s{3},s{4}];
    
    ptsFile = join([folder '/cake_nodes.m'])
    fid=fopen(ptsFile);
    s=textscan(fid,'%f %f %f %f');
    fclose(fid);
    pts = [s{2},s{3},s{4}];
    
    mesh = FullMesh(tetr,pts,neumannIdentifier);
end


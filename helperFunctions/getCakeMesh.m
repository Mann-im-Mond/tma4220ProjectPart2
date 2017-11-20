function mesh = getCakeMesh(folder, neumannIdentifier, radius)
    addpath('../utils/')
    
    tetrFile = [folder '/cake_tetr.m'];
    fid=fopen(tetrFile);
    s=textscan(fid,'%f %f %f %f %f %f');
    fclose(fid);
    tetr = [s{1},s{2},s{3},s{4}];
    rod = [s{5}];

    
    ptsFile = [folder '/cake_nodes.m'];
    fid=fopen(ptsFile);
    s=textscan(fid,'%f %f %f %f');
    fclose(fid);
    pts = [s{2},s{3},s{4}];
    pts = pts*(radius);
    
    mesh = FullMesh(tetr,pts,rod,neumannIdentifier);
end


function Widefield_WriteTiff(outFile,Frames)
% Code to write 3D matrix that is a stack of images into a .tif file.
% Usage: Widefield_WriteTiff(outFile,Frames). 
% outFile is the filename with or without path, excluding .tif at the end.
% Frames is the input matrix that should be written down.

% Write first frame then do the rest of the stack
imwrite(uint16(Frames(:,:,1)), [outFile '.tif'], 'TIF', ...
    'Resolution', [size(Frames, 2) size(Frames, 1)], 'Compression', 'none');

if size(Frames,3) > 1
    for f = 2:size(Frames,3)
        % Report on progress
        if mod(f, 100) == 0
            fprintf('%d ', f);
        end
        if mod(f, 100) == 0
            fprintf('\n');
        end
        
        % Write next frame
        imwrite(uint16(Frames(:,:,f)), [outFile '.tif'], 'TIF', ...
            'Resolution', [size(Frames, 2) size(Frames, 1)], 'Compression', 'none', ...
            'WriteMode', 'append');
    end
end
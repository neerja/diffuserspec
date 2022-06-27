function [data] = readSPF2(filename)

    fid = fopen(filename);
    fseek(fid, 524, 'bof');
    Length = fread(fid, 1, 'uint32');                %file stores little-endian
    fseek(fid, -Length*4*2,'eof');      %4 bytes per entry, 2 columns
    data = reshape( fread(fid, '*single'), [], 2);   %file stores little-endian
    fclose(fid);
    
end
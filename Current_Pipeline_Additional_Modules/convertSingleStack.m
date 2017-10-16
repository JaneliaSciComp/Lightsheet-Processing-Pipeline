inputFileName  = 'Stack.tif';
outputFileName = 'Stack.inr';
resolution     = [6.5/16 6.5/16 6.5/16*5];
machineFormat  = 'l'; % use 'l' for correct byte order in block-matching software (default)
                      % use 'b' for correct byte order in Fiji

stack = readImage(inputFileName);
[xSize, ySize, zSize] = size(stack);

fid = fopen(outputFileName, 'w');
headerLength = 0;

headerA = '#INRIMAGE-4#{';
fwrite(fid, headerA);
fprintf(fid, '\n');
headerLength = headerLength + length(headerA) + 1;

headerB = ['XDIM=' num2str(xSize)];
fwrite(fid, headerB);
fprintf(fid, '\n');
headerLength = headerLength + length(headerB) + 1;

headerC = ['YDIM=' num2str(ySize)];
fwrite(fid, headerC);
fprintf(fid, '\n');
headerLength = headerLength + length(headerC) + 1;

headerD = ['ZDIM=' num2str(zSize)];
fwrite(fid, headerD);
fprintf(fid, '\n');
headerLength = headerLength + length(headerD) + 1;

headerE = 'TYPE=unsigned fixed'; % can also be float
fwrite(fid, headerE);
fprintf(fid, '\n');
headerLength = headerLength + length(headerE) + 1;

headerF = 'PIXSIZE=16 bits'; % can be 8, 16, 32, 64 bit
fwrite(fid, headerF);
fprintf(fid, '\n');
headerLength = headerLength + length(headerF) + 1;

headerG = 'CPU=decm';
fwrite(fid, headerG);
fprintf(fid, '\n');
headerLength = headerLength + length(headerG) + 1;

headerH = ['VX=' num2str(resolution(1))];
fwrite(fid, headerH);
fprintf(fid, '\n');
headerLength = headerLength + length(headerH) + 1;

headerI = ['VY=' num2str(resolution(2))];
fwrite(fid, headerI);
fprintf(fid, '\n');
headerLength = headerLength + length(headerI) + 1;

headerJ = ['VZ=' num2str(resolution(3))];
fwrite(fid, headerJ);
fprintf(fid, '\n');
headerLength = headerLength + length(headerJ) + 1;

headerK = '#GEOMETRY=CARTESIAN';
fwrite(fid, headerK);
fprintf(fid, '\n');
headerLength = headerLength + length(headerK) + 1;

if mod(headerLength + 4, 256) ~= 0 % total header length has to be multiple of 256
    fillerElements = 256 - mod(headerLength + 4, 256);
    for i = 1:fillerElements
        fprintf(fid, '\n');
    end;
end;

fwrite(fid, '##}');
fprintf(fid, '\n');

fwrite(fid, stack(:), 'uint16', machineFormat);

fclose(fid);
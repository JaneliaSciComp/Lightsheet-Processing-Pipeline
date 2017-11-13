%% parameters

timepoints   = 0:500;

rootFolder   = ['U:' filesep 'SiMView2.Processing' filesep '12-06-12' filesep 'Dme_E1_His2ARFP_01_20120612_175853.corrected' filesep 'Results' filesep 'TimeFused.Interpolated'];
inputHeader1 = 'Dme_E1_His2ARFP.TM';
inputHeader2 = '_timeFused_blending';
inputHeader3 = 'SPC0_TM';
inputFooter  = '_CM0_CM1_CHN00_CHN01.fusedStack';
timeStamp    = '%.4d';

embryoVector = [1209, 295, 204; 84, 250, 551]; % start-x,y,z; stop-x,y,z (Matlab coordinates)
percentRange = [17-0.5, 17+0.5];               % start/stop percentage along embryoVector defining the axial extent of slicing volume

interMode    = 0; % 0 for linear, 1 for spline, 2 for cubic

inputType    = 0; % 0: input data in KLB format
                  % 1: input data in JP2 format
                  % 2: input data in TIF format

saveStacks   = 0;
                       
%% main loop

switch inputType
    case 0
        inputExtension = '.klb';
    case 1
        inputExtension = '.jp2';
    case 2
        inputExtension = '.tif';
end;

if interMode == 0
    interpolationKeyword = 'linear';
elseif interMode == 1
    interpolationKeyword = 'spline';
else
    interpolationKeyword = 'cubic';
end;

outputFolder = [rootFolder '.Sliced'];
if exist(outputFolder, 'dir') ~= 7
    mkdir(outputFolder);
end;

disp(' ');
disp('Computing slicing geometry');
disp(' ');

inputName = [rootFolder filesep '' inputHeader1 num2str(timepoints(1), timeStamp) inputHeader2 filesep '' inputHeader3 num2str(timepoints(1), timeStamp) inputFooter inputExtension];

stack = readImage(inputName);
xSize = size(stack, 1);
ySize = size(stack, 2);
zSize = size(stack, 3);

clear stack;

% determine new coordinate system (u, v, w)

u = [embryoVector(2, 1) - embryoVector(1, 1), embryoVector(2, 2) - embryoVector(1, 2), embryoVector(2, 3) - embryoVector(1, 3)];
u = u ./ sqrt(u * u');

v = [0, 1, - u(2)/u(3)];
v = v ./ sqrt(v * v');

w = cross(u, v);

a = embryoVector(1, :) + (mean(percentRange) / 100) * (embryoVector(2, :) - embryoVector(1, :));

disp('Slicing coordinate system:');
disp(['u = (' num2str(u(1), '%.4f') ', ' num2str(u(2), '%.4f') ', ' num2str(u(3), '%.4f') ')']);
disp(['v = (' num2str(v(1), '%.4f') ', ' num2str(v(2), '%.4f') ', ' num2str(v(3), '%.4f') ')']);
disp(['w = (' num2str(w(1), '%.4f') ', ' num2str(w(2), '%.4f') ', ' num2str(w(3), '%.4f') ')']);
disp(' ');
disp('Anchor vector:');
disp(['a = (' num2str(a(1), '%.4f') ', ' num2str(a(2), '%.4f') ', ' num2str(a(3), '%.4f') ')']);

embryoVectorRelative = embryoVector(2, :) - embryoVector(1, :);
uRadius = sqrt(embryoVectorRelative * embryoVectorRelative') * (percentRange(2) - percentRange(1)) / 200;

cu = ceil(max(abs((ySize - a(2)) / v(2)), abs((a(2) - 1) / v(2))));
cv = ceil(uRadius);
cw = ceil(max(abs((zSize - a(3)) / w(3)), abs((a(3) - 1) / w(3))));

[CU, CV, CW] = meshgrid(-cv:cv, -cu:cu, -cw:cw);
XYZ = bsxfun(@plus, bsxfun(@times, CU(:), u) + bsxfun(@times, CV(:), v) + bsxfun(@times, CW(:), w), a);

% Visualize XYZ points:
% plot3(XYZ(:, 1), XYZ(:, 2), XYZ(:, 3), '.'); xlabel('x'); ylabel('y'); zlabel('z'); axis equal;

[inputX, inputY, inputZ] = meshgrid(1:ySize, 1:xSize, 1:zSize);
sz = [2*cu + 1, 2*cv + 1, 2*cw + 1];

disp(' ');

tic;

for currentTP = 1:numel(timepoints)
    disp(['Slicing time point ' num2str(timepoints(currentTP))]);
    
    fileName = [...
        rootFolder filesep '' ...
        inputHeader1 num2str(timepoints(currentTP), timeStamp) inputHeader2 filesep '' ...
        inputHeader3 num2str(timepoints(currentTP), timeStamp) inputFooter inputExtension];
    
    stack = readImage(fileName);
    
    slice = uint16(interp3(inputX, inputY, inputZ, single(stack), reshape(XYZ(:, 2), sz), reshape(XYZ(:, 1), sz), reshape(XYZ(:, 3), sz), interpolationKeyword));
    
    rotatedSlice = zeros(size(slice, 3), size(slice, 1), size(slice, 2), 'uint16');
    for z = 1:size(slice, 2)
        rotatedSlice(:, :, z) = rot90(squeeze(slice(:, z, :)));
    end;
    
    outputSliceName = [outputFolder filesep '' inputHeader3 num2str(timepoints(currentTP), timeStamp) inputFooter '_sliced' inputExtension];
    outputProjectionName = [outputFolder filesep '' inputHeader3 num2str(timepoints(currentTP), timeStamp) inputFooter '_sliced.projection' inputExtension];
    
    if saveStacks
        writeImage(rotatedSlice, outputSliceName);
    end;
    writeImage(squeeze(max(rotatedSlice, [], 3)), outputProjectionName);
end;

elapsedTime = toc;

disp(' ');

disp(['Total processing time for stack slicing: ' num2str(elapsedTime / 60, '%.2f') ' min']);

disp(' ');
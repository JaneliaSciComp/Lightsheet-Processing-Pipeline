function Iout = imTranslate_Z(I, x, y, z)
sizX = size(I, 1);
sizY = size(I, 2);
sizZ = size(I, 3);
x = round(x);
y = round(y);

if x<0
    Iout = [zeros(-x, sizY, sizZ); I(1:end+x, :, :)];
else
    Iout = [I(x+1:end, :, :); zeros(x, sizY, sizZ)];
end

if y<0
    Iout = [zeros(sizX, -y, sizZ) Iout(:, 1:end+y, :)];
else
    Iout = [Iout(:, y+1:end, :) zeros(sizX, y, sizZ)];
end

ROI = (1:sizZ) + z;
Iout = permute(Iout, [3 1 2]);
Iout = interp1(1:sizZ, double(Iout), ROI);
Iout(isnan(Iout)) = 0;
Iout = uint16(permute(Iout, [2 3 1]));
end
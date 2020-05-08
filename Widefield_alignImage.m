function [error, shift, cPic] = Widefield_alignImage(alignPic, cPic, rotAngle, cropRange)

cPic = imrotate(cPic,rotAngle,'crop', 'bilinear');

if cropRange > 0
imCenter = round(size(alignPic) / 2);
alignPic = alignPic(imCenter(1)-cropRange+1 : imCenter(1)+cropRange, imCenter(2)-cropRange+1 : imCenter(2)+cropRange);
cPic = cPic(imCenter(1)-cropRange+1 : imCenter(1)+cropRange, imCenter(2)-cropRange+1 : imCenter(2)+cropRange);
end

[stats, cPic] = Widefield_dftregistration(fft2(alignPic), fft2(cPic), 10);
stats = gather(stats);
error = stats(1);
shift = stats(3:4);
cPic = abs(ifft2(cPic));


end
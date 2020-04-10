function result = paddingimage(image,multiple)%multiple为扩充对应的倍数，这里为8
%图像扩充，将图像宽与高扩充到8的倍数，以便于分割
[height,width] = size(image);
if mod(height,multiple) ~= 0
    addheight = multiple - mod(height,multiple);%扩充高数
else
    addheight = 0;
end
if mod(width,multiple) ~= 0
    addwidth = multiple - mod(width,multiple);%扩充宽数
else
    addwidth = 0;
end
result = padarray(image,[addheight,addwidth],'replicate','post');%复制外边界，右下填充保证数目正确
end

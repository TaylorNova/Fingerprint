function enhance_fingerprint3()
%%
%读取图片并扩展
I = preprocess();%图像预处理
multiple = 4;%扩展对应的倍数
windowsize = 8;%DFT窗口大小
padsize = (windowsize-multiple)/2;%扩展大小
centerpoint = windowsize/2 + 1;%DFT中心点
I = paddingimage(I,multiple);%扩展到倍数
Ipad = padarray(I,[padsize,padsize],'replicate','both');%扩展以便于更大的DFT
Ipad2 = padarray(I,[18,18],'replicate','both');%扩展以便于更大范围的灰度均值和方差的计算
%%
%分割
[height,width] = size(I);%图像尺寸
heinum = height/multiple;%分块数目，宽高两个方向
widnum = width/multiple;
divideimage = zeros(heinum,widnum,windowsize,windowsize);%分割的图像
dftmatrix = zeros(heinum,widnum,windowsize,windowsize);%幅度谱
sign = zeros(heinum,widnum);%分割标志
threshold_lowfre = 0.8;%分割低频部分幅度占比的阈值
threshold_meangray = 0.45;%灰度均值的阈值
threshold_stdgray = 0.136;%灰度方差的阈值
%%
%通过幅度谱低频部分占比进行分割
for i = 1:heinum
    for j = 1:widnum
        divideimage(i,j,:,:) = Ipad((i-1)*multiple+1:i*multiple+windowsize-multiple,(j-1)*multiple+1:j*multiple+windowsize-multiple);%图像分块存储
        dftmatrix(i,j,:,:) = abs(fftshift(fft2(Ipad((i-1)*multiple+1:i*multiple+windowsize-multiple,(j-1)*multiple+1:j*multiple+windowsize-multiple))));%图像DFT
        persent_lowfre = sum(sum(dftmatrix(i,j,centerpoint-1:centerpoint+1,centerpoint-1:centerpoint+1)))/sum(sum(dftmatrix(i,j,:,:)));%低频部分占比
        meangray = mean2(Ipad2((i-1)*multiple+1:i*multiple+36,(j-1)*multiple+1:j*multiple+36));%40*40灰度均值
        stdgray = std2(Ipad2((i-1)*multiple+1:i*multiple+36,(j-1)*multiple+1:j*multiple+36));%40*40灰度方差
        if persent_lowfre > threshold_lowfre
            sign(i,j) = 0;
        else
            if meangray < threshold_meangray && stdgray < threshold_stdgray
                sign(i,j) = 1;
            else
                sign(i,j) = 0;
            end
        end
    end
end
%形态学处理
%求最大连通域
CC = bwconncomp(sign);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,num] = size(numPixels);
sign2 = false(size(sign));
for i = 1:num
    if numPixels(i)>5000%如果连通域像素个数大于5000
        sign2(CC.PixelIdxList{i}) = 1;
    end
end
%膨胀以扩张边缘
circle = strel('disk',3);%圆形mask
sign3 = imdilate(sign2,circle);
%查看分割情况
divide = I;
for i = 1:heinum
    for j = 1:widnum
        if sign3(i,j) == 0
            divide((i-1)*multiple+1:i*multiple,(j-1)*multiple+1:j*multiple) = 0;%背景部分置零
        end
    end
end
figure(1);
I0 = imread('3.bmp');
subplot(2,3,1);imshow(I0);title('原图');
subplot(2,3,2);imshow(divide);title('分割结果');
%%
%方向图的确定
cita = zeros(heinum,widnum);%方向图
for i = 1:heinum
    for j = 1:widnum
        if sign3(i,j) == 1
        dftcopy = zeros(windowsize,windowsize);%拷贝
        for m = 1:windowsize
            for n = 1:windowsize
               dftcopy(m,n) =  dftmatrix(i,j,m,n);
            end
        end
        dftcopy(centerpoint,centerpoint) = 0;%直流分量置零
        maxnum = max(max(dftcopy));
        maxpoint = zeros(2,2);%幅度最大值点
        tempnum = 1;
        for m = 1:windowsize
            for n = 1:windowsize
                if dftcopy(m,n) == maxnum %寻找最大值点
                    maxpoint(tempnum,1) = m;
                    maxpoint(tempnum,2) = n;
                    tempnum = tempnum + 1;
                end
            end
        end
        cita(i,j) = 0.5*pi + atan((maxpoint(1,2)-maxpoint(2,2))/(maxpoint(1,1)-maxpoint(2,1)));%方向计算
        end
    end
end
subplot(2,3,4);quiver(flip(cos(cita - 0.5*pi)), flip(sin(cita - 0.5*pi)));title('方向图');
%%
%频率图的确定
frequency = zeros(heinum,widnum);%频率
nmbda = zeros(heinum,widnum);%波长
for i = 1:heinum
    for j = 1:widnum
        if sign3(i,j) == 1%判断是否为指纹区域
        dftcopy = zeros(windowsize,windowsize);%拷贝
        for m = 1:windowsize
            for n = 1:windowsize
               dftcopy(m,n) =  dftmatrix(i,j,m,n);
            end
        end
        dftcopy(centerpoint,centerpoint) = 0;%删除直流分量
        maxnum = max(max(dftcopy));
        maxpoint = zeros(2,2);%幅度最大值点，对应指纹方向与频率
        tempnum = 1;
        for m = 1:windowsize
            for n = 1:windowsize
                if dftcopy(m,n) == maxnum %寻找最大值点
                    maxpoint(tempnum,1) = m;
                    maxpoint(tempnum,2) = n;
                    tempnum = tempnum + 1;
                end
            end
        end
        %计算频率与波长
        distanceofmax = ((maxpoint(1,2)-centerpoint)^2+(maxpoint(1,1)-centerpoint)^2)^0.5;%幅度最大值点到（9,9）的距离
        frequency(i,j) = 2*pi*distanceofmax/windowsize;%频率
        nmbda(i,j) = windowsize/distanceofmax;%波长
        end
    end
end
%%
%方向图平滑
coscita = zeros(heinum,widnum);
sincita = zeros(heinum,widnum);
for i = 1:heinum
    for j = 1:widnum
        coscita(i,j) = cos(2*cita(i,j));
        sincita(i,j) = sin(2*cita(i,j));
    end
end
w0 = fspecial('gaussian',3,0.5);
%分别采用高斯均值滤波
coscitasmooth = imfilter(coscita,w0);
sincitasmooth = imfilter(sincita,w0);
citasmooth = 0.5*atan2(sincitasmooth,coscitasmooth);
subplot(2,3,5);quiver(flip(cos(citasmooth - 0.5*pi)), flip(sin(citasmooth - 0.5*pi)));title('平滑后的方向图');
%频率图平滑，实际是波长图的平滑
nmbdasmooth = imfilter(nmbda,w0);
subplot(2,3,3);imshow(nmbdasmooth,[]);title('平滑后的波长图');
%%
%滤波
padlen = 1;
Ipad2 = padarray(I,[padlen,padlen],'replicate','both');%扩展
[height2,width2] = size(Ipad2);
filterresult = zeros(height2,width2);
imageresult = zeros(height,width);
citause = citasmooth.*180./pi;%角度变换为角度制
for i = 1:heinum
    for j = 1:widnum
        if sign3(i,j) == 1
        [mag,phase] = imgaborfilt(Ipad2((i-1)*multiple+1:i*multiple+padlen*2,(j-1)*multiple+1:j*multiple+padlen*2), nmbdasmooth(i,j),citause(i,j));%Gabor滤波
        temp = mag .*cos(phase);
        filterresult((i-1)*multiple+1:i*multiple+padlen*2,(j-1)*multiple+1:j*multiple+padlen*2) = temp;
        end
    end
end
imageresult = filterresult(padlen+1:height2-padlen,padlen+1:width2-padlen);
subplot(2,3,6);imshow(imageresult);title('处理结果');
end
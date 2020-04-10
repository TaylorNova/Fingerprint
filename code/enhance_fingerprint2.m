function enhance_fingerprint2()
%%
%读取图片并扩展
I = imread('23_2.bmp');
multiple = 8;%扩展对应的倍数
I = paddingimage(I,multiple);%扩展到8的倍数
Ipad = padarray(I,[12,12],'replicate','both');%扩展以便于32*32的DFT
%%
%分割
[height,width] = size(I);%图像尺寸
heinum = height/multiple;%分块数目，宽高两个方向
widnum = width/multiple;
divideimage = zeros(heinum,widnum,32,32);%分割的图像
dftmatrix = zeros(heinum,widnum,32,32);%幅度谱
sign = zeros(heinum,widnum);%分割标志
threshold_lowfre = 0.36;%分割低频部分幅度占比的阈值
%%
%通过幅度谱低频部分占比进行分割
for i = 1:heinum
    for j = 1:widnum
        divideimage(i,j,:,:) = Ipad((i-1)*multiple+1:i*multiple+24,(j-1)*multiple+1:j*multiple+24);%图像分块存储
        dftmatrix(i,j,:,:) = abs(fftshift(fft2(Ipad((i-1)*multiple+1:i*multiple+24,(j-1)*multiple+1:j*multiple+24))));%图像DFT
        persent_lowfre = sum(sum(dftmatrix(i,j,15:20,15:20)))/sum(sum(dftmatrix(i,j,:,:)));%低频部分占比
        if persent_lowfre > threshold_lowfre
            sign(i,j) = 0;
        else
            sign(i,j) = 1;
        end
    end
end
%查看分割情况
divide = I;
for i = 1:heinum
    for j = 1:widnum
        if sign(i,j) == 0
            divide((i-1)*multiple+1:i*multiple,(j-1)*multiple+1:j*multiple) = 0;%背景部分置零
        end
    end
end
figure(1);
subplot(2,3,1);imshow(I);title('原图');
subplot(2,3,2);imshow(divide);title('分割结果');
%%
%方向图的确定
cita = zeros(heinum,widnum);%方向图
for i = 1:heinum
    for j = 1:widnum
        if sign(i,j) == 1
        dftcopy = zeros(32,32);%拷贝
        for m = 1:32
            for n = 1:32
               dftcopy(m,n) =  dftmatrix(i,j,m,n);
            end
        end
        dftcopy(17,17) = 0;%直流分量置零
        maxnum = max(max(dftcopy));
        maxpoint = zeros(2,2);%幅度最大值点
        tempnum = 1;
        for m = 1:32
            for n = 1:32
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
        if sign(i,j) == 1%判断是否为指纹区域
        dftcopy = zeros(32,32);%拷贝
        for m = 1:32
            for n = 1:32
               dftcopy(m,n) =  dftmatrix(i,j,m,n);
            end
        end
        dftcopy(17,17) = 0;%删除直流分量
        maxnum = max(max(dftcopy));
        maxpoint = zeros(2,2);%幅度最大值点，对应指纹方向与频率
        tempnum = 1;
        for m = 1:32
            for n = 1:32
                if dftcopy(m,n) == maxnum %寻找最大值点
                    maxpoint(tempnum,1) = m;
                    maxpoint(tempnum,2) = n;
                    tempnum = tempnum + 1;
                end
            end
        end
        %计算频率与波长
        distanceofmax = ((maxpoint(1,2)-17)^2+(maxpoint(1,1)-17)^2)^0.5;%幅度最大值点到（17,17）的距离
        frequency(i,j) = 2*pi*distanceofmax/32;%频率
        nmbda(i,j) = 32/distanceofmax;%波长
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
Ipad2 = padarray(I,[4,4],'replicate','both');%扩展
[height2,width2] = size(Ipad2);
filterresult = zeros(height2,width2);
imageresult = zeros(height,width);
citause = citasmooth.*180./pi;%角度变换为角度制
for i = 1:heinum
    for j = 1:widnum
        if sign(i,j) == 1
        [mag,phase] = imgaborfilt(Ipad2((i-1)*multiple+1:i*multiple+8,(j-1)*multiple+1:j*multiple+8), nmbdasmooth(i,j),citause(i,j));%Gabor滤波
        temp = mag .*cos(phase);
        filterresult((i-1)*multiple+1:i*multiple+8,(j-1)*multiple+1:j*multiple+8) = temp;
        end
    end
end
imageresult = filterresult(5:height2-4,5:width2-4);
w1 = fspecial('gaussian', [8,8],1);%高斯滤波，使得效果更佳
imageresult = imfilter(imageresult,w1);
subplot(2,3,6);imshow(imageresult);title('处理结果');
end
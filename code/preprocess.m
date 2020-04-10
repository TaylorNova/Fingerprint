function result = preprocess()
%对第三幅图进行预处理
I = imread('3.bmp');
multiple = 8;%扩展对应的倍数
I = paddingimage(I,multiple);%按倍数扩展
dftI = fftshift(fft2(I));%傅里叶变换
margin = log(1+abs(dftI));%幅度谱
center = [385,401];%幅度谱中心点
%手动找出的噪声点坐标
x = [112,295,202,295,21,203,111,21];
y = [306,215,494,679,120,30,772,586];
len = 7;%置零区域大小
% figure(1);
% subplot(1,2,1);imshow(margin,[]); %显示原幅度谱
for i = 1:8
    dftI(x(i)-len:x(i)+len,y(i)-len:y(i)+len) = 0;
    dftI(2*center(1)-x(i)-len:2*center(1)-x(i)+len,2*center(2)-y(i)-len:2*center(2)-y(i)+len) = 0;
end
% margin = log(1+abs(dftI));
% subplot(1,2,2);imshow(margin,[]);%显示处理后的幅度谱
Itemp = ifft2(fftshift(dftI));%傅里叶反变换
% figure(2);
% imshow(Itemp,[]);%显示频域处理结果
Idivide = Itemp(220:539,306:505);%切割出指纹图像
MinValue = min(min(Idivide));
MaxValue = max(max(Idivide));
Idivide =(Idivide-MinValue)/(MaxValue-MinValue);
% figure(3);
% imshow(Idivide,[]);%显示指纹图像
result = Idivide;
end
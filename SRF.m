clear;clc;close all;warning off
%% 平均最后n个光谱数据文件
file_list1 = dir('F:\课题\数据\武\SO2标准谱\11.19-50cm-ppb\bd\*.txt');
file_list2 = dir('F:\课题\数据\武\SO2标准谱\11.19-50cm-ppb\12ppm\*.txt');
file_list3 = dir('F:\课题\数据\武\SO2标准谱\11.19-50cm-ppb\1000\*.txt');
u0=196;v0=225;%测量波段
u1=196;v1=225;%拟合波段
%标注浓度，测量浓度
S_C=12000;M_C=1000;xishu=S_C/M_C;
len1 = size(file_list1, 1);
len2 = size(file_list2, 1);
len3 = size(file_list3, 1);
Q1=[];Q2=[];Q3=[];
%背底谱采用最后10个文件
for i1=(len1-9):1:len1
    file_list1(i1).name;
    %     fileID=fopen(file_list(i).name);(路径不全,使用strcat连接文本)
    fileID1=fopen(strcat('F:\课题\数据\武\SO2标准谱\11.19-50cm-ppb\bd\',file_list1(i1).name));
    A1=textscan(fileID1,'%f%f');
    %     figure
    %     plot(A1{1,1},A1{1,2});
    %     hold on;循环打开n个光谱图
    Q1=[Q1;A1];
    fclose all;
end
%平均了最后10个文件的背底谱数据
beidibochang=Q1{1,1};
a1=zeros(length(beidibochang),1);
for i2=1:length(Q1)
    a1=a1+Q1{i2,2};
end
beidiguangqiang=a1/length(Q1);
beidishuju=[beidibochang beidiguangqiang];
%差分谱1采用最后10个
for i3=(len2-9):1:len2
    file_list2(i3).name;
    fileID2=fopen(strcat('F:\课题\数据\武\SO2标准谱\11.19-50cm-ppb\12ppm\',file_list2(i3).name));
    A2=textscan(fileID2,'%f%f');
    %     figure
    %     plot(A2{1,1},A2{1,2});
    %     hold on;
    Q2=[Q2;A2];
    fclose all;
end
celiangbochang1=Q2{1,1};
a2=zeros(length(celiangbochang1),1);
for i4=1:length(Q2)
    a2=a2+Q2{i4,2};
end
%平均了最后10个文件的测量谱1数据
celiangguangqiang1=a2/length(Q2);
celiangshuju1=[celiangbochang1 celiangguangqiang1];
%差分谱2采用最后10个
for i5=(len3-29):1:len3
    file_list3(i5).name;
    fileID3=fopen(strcat('F:\课题\数据\武\SO2标准谱\11.19-50cm-ppb\1000\',file_list3(i5).name));
    A3=textscan(fileID3,'%f%f');
    %     figure
    %     plot(A3{1,1},A3{1,2});
    %     hold on;
    Q3=[Q3;A3];
    fclose all;
end
celiangbochang2=Q3{1,1};
a3=zeros(length(celiangbochang2),1);
for i6=1:length(Q3)
    a3=a3+Q3{i6,2};
end
%平均了最后10个文件的测量谱2数据
celiangguangqiang2=a3/length(Q3);
celiangshuju2=[celiangbochang2 celiangguangqiang2];
fclose all;
%% 测量波长与拟合波长选取
array=beidishuju(:,1);
value1=u0;value2=v0;value3=u1;value4=v1;
u00 = findClosestNum(array, value1);
v00 = findClosestNum(array, value2);
u11 = findClosestNum(array, value3);
v11 = findClosestNum(array, value4);
wave1=[u00,v00];
wave2=[u11,v11];
% fprintf('吸收波长: ');
% fprintf('%f ', wave1);
% fprintf('\n');
% fprintf('拟合波长: ');
% fprintf('%f ', wave2);
% fprintf('\n');  % 可选，添加换行符
%% 差分谱1
xishoupu2=celiangguangqiang1./beidiguangqiang;
u01 = find(celiangbochang1==u00);
v01 = find(celiangbochang1==v00);
yongdexishoupu1=xishoupu2(u01:v01,1);
yongdebochang1=celiangbochang1(u01:v01,1);%用的吸收谱波长选取

u12 = find(celiangbochang1==u11);
v12 = find(celiangbochang1==v11);
nihebochang1=celiangbochang1(u12:v12,1);
nihexishoupu1=xishoupu2(u12:v12,1);
manbianxishou=polyfit(nihebochang1,nihexishoupu1,4);
manbianxishou1=polyval(manbianxishou,nihebochang1);
manbianpu1=[nihebochang1,manbianxishou1];
u22 = find(nihebochang1==u00);
v22 = find(nihebochang1==v00);
yongdemanbian1=manbianpu1(u22:v22,:);%用的慢变拟合波段选取

kuaibianxishoupu2=yongdexishoupu1./yongdemanbian1(:,2);
chafenpu1=log(kuaibianxishoupu2);%取对数
%% 差分谱2
xishoupu2=celiangguangqiang2./beidiguangqiang;
yongdexishoupu2=xishoupu2(u01:v01,1);
yongdebochang2=celiangbochang1(u01:v01,1);%用的吸收谱波长选取

nihebochang2=celiangbochang1(u12:v12,1);
nihexishoupu2=xishoupu2(u12:v12,1);
manbianxishou=polyfit(nihebochang1,nihexishoupu2,4);
manbianxishou2=polyval(manbianxishou,nihebochang1);
manbianpu2=[nihebochang2,manbianxishou2];%用的慢变拟合波段选取

yongdemanbian2=manbianpu2(u22:v22,:);
kuaibianxishoupu2=yongdexishoupu2./yongdemanbian2(:,2);
chafenpu2=log(kuaibianxishoupu2);%取对数
%% 图1两个差分谱
y1=chafenpu1;
y2=chafenpu2;
plot(yongdebochang1,chafenpu1,'r')
hold on
plot(yongdebochang2,chafenpu2,'b')
title('两组差分谱');
xlabel('波长λ(nm)');
ylabel('差分吸收光谱(a.u)');
legend('标准浓度差分谱','待测浓度差分谱')
wen_01=sum(y1.*y2);
wen_02=sum(y1.*y1);
wen_k2=wen_01/wen_02;
wen_c2=wen_k2*S_C;
% fprintf('%0.20f\n', wen_k2)
%% 图2 插值
y_0=y1;
y_1=y2;
%% 差分谱选取排序
bochang=yongdebochang1;
m=length(bochang);
a=1:m;
Y=[y_0 y_1 a'];
Y0=[y_0 y_1];
Y1=sortrows(Y,1);%标准谱从小到大排序，测量谱跟着变动
Y01=sortrows(Y0,1);
%% 图3 标准谱线性，测量谱变动
X=[];
X(1)=0;
for p1=2:length(Y1)
    X(p1)=1*(Y1(p1,1)-Y1(p1-1,1))+X(p1-1);
    p1=p1+1;
end
XX=X';
lunwen02=[X' Y01];
% figure
% plot(X,Y01,'-o');
%% 图4 线性标准谱，拟合测量谱
%拟合变动后的测量谱
cc=Y1(:,2)';
yn=polyfit(X,cc,1);
ynn=polyval(yn,X);
%Y1排序
%Y2标准谱，测量谱，
%Y3标准谱，拟合后的测量谱，
Y2=[Y1(:,1) cc'];
Y3=[Y1(:,1) ynn'];
ynn0=ynn';
figure
plot(X,Y3);
%% 图5 标准谱横坐标，测量谱做纵坐标
figure
plot(Y1(:,1),cc','r');
hold on
plot(Y1(:,1),ynn','b');%拟合完
%% 图6 逆重构
Y31=[Y1(:,1) ynn' Y1(:,3)];
Y32=sortrows(Y31,3);
Y33=[Y32(:,1) Y32(:,2)];
bochang=bochang';
figure
plot(bochang,y_1,'r');hold on
% plot(bochang,Y33(:,1),'-ko');hold on
plot(bochang,Y33(:,2),'b');
title('三组差分谱');
xlabel('波长λ(nm)');
ylabel('差分吸收光谱(a.u)');
legend('待测浓度差分谱','逆重构差分谱')
fclose all;

wen_S=Y33(:,1);
wen_M=Y33(:,2);

yn2=polyfit(wen_S,wen_M,1);


% wen_11=sum(Y33(:,1).*Y33(:,2))
% wen_12=sum(Y33(:,1).^2)
% wen_k1=wen_11/wen_12
wen_k=sum(wen_S.*wen_M)/sum(wen_S.^2);
% wen_c=S_C*wen_k
% fprintf('%0.20f\n', wen_k)
% error=(abs(wen_c-M_C)/M_C)*100
% OP1=sum(abs(y1));
% OP_01=sum(abs(y2));
% OP_02=sum(abs(Y33(:,2)));
% A_OP=[OP1 OP_01 OP_02];
format compact
A_C1=yn(1)*S_C
% fprintf('%0.20f\n', A_C1)
% A_C2=yn2(1)*S_C
% fprintf('%0.20f\n', A_C2)
error=(abs(A_C1-M_C)/M_C)*100
A_A=[A_C1 error];

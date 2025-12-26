clc
clear
ff=fopen('C:\Users\Administrator\Desktop\文章各个图\supplementary figures\supplefig5\新增加的动物数据\DMS文件\DMS-annotated\all_transcripts_UTR5_human.fa.out');
data=cell(1,1);
t=0;
while ~feof(ff)
    f=fgetl(ff);
    t=t+1;data{t,1}=f;
    f=fgetl(ff);
    data{t,2}=f;
    f=fgetl(ff);
    data{t,3}=f;
    data{t,4}=(f((length(data{t,2})+3):(length(f)-1)));
    data{t,5}=data{t,3}(1:length(data{t,2}));
    f=fgetl(ff);
    f=fgetl(ff);
    f=fgetl(ff);
    data{t,5}=data{t,3}(1:length(data{t,2}));
end
fclose(ff)

for i=1:size(data,1)
w=strfind(data{i,1},' ');
s=data{i,1}(1:w(1)-1);
data{i,1}=s;

end



ff=fopen('C:\Users\Administrator\Desktop\文章各个图\supplementary figures\supplefig5\新增加的动物数据\DMS文件\DMS-annotated\humaninvivo5utr.out');
datainvivo=cell(1,1);
t=0;
while ~feof(ff)
    f=fgetl(ff);
    t=t+1;datainvivo{t,1}=f;
    f=fgetl(ff);
    datainvivo{t,2}=f;
    f=fgetl(ff);
    datainvivo{t,3}=f;
    datainvivo{t,4}=(f((length(datainvivo{t,2})+3):(length(f)-1)));
    datainvivo{t,5}=datainvivo{t,3}(1:length(datainvivo{t,2}));
    f=fgetl(ff);
    f=fgetl(ff);
    f=fgetl(ff);
    datainvivo{t,5}=datainvivo{t,3}(1:length(datainvivo{t,2}));
end
fclose(ff)



% 一共有5003个
[x,w1,w2]=intersect(datainvivo(:,1),data(:,1));

part1=datainvivo(w1,5);
part2=data(w2,5);
stemratio=zeros(length(part1),2);
for i=1:length(part1)
stemratio(i,1)=1-length(strfind(part1{i},'.'))/length(part1{i});
stemratio(i,2)=1-length(strfind(part2{i},'.'))/length(part2{i});

end

xlswrite('stemratio_human_invivo_insilico_5utr.xlsx',stemratio)

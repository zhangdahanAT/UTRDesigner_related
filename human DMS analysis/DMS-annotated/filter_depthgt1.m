clc
clear
load("humanfiveprimeutrmutations2.mat")



a=cell2table(fiveprimeutrmutations2(:,1));
transcriptids=a.Variables;
transcriptidsunique = unique(transcriptids);
output=zeros(size(transcriptidsunique,1),3);
for i=1:size(transcriptidsunique,1)
    i
    w=find(ismember(transcriptids,transcriptidsunique(i)));
    mt=sum(cell2mat(fiveprimeutrmutations2(w,9)));
    depth=fiveprimeutrmutations2{w(1),4};
    output(i,1) = mt;
    output(i,2) = depth;
    output(i,3) = mt/depth; 

end

length(find(output(:,3)>1))
% 5198个 mutation site sum / utr length >1 
transcriptids_qualified=transcriptidsunique(find(output(:,3)>1));
[fiveutrhead,fiveutrseq]=fastaread('C:\Users\Administrator\Desktop\文章各个图\supplementary figures\supplefig5\新增加的动物数据\DMS文件\DMS-annotated\all_transcripts_5UTR_annotated.fa');
% 筛选出DMS-seq支持的5'UTR,五千多个。
[~,w1,w2]=intersect(fiveutrhead,transcriptids_qualified);






a1=cell2table(fiveprimeutrmutations2(:,4));
utrlen=a1.Variables;

a2=cell2table(fiveprimeutrmutations2(:,6));
pos=a2.Variables;

a3=cell2table(fiveprimeutrmutations2(:,10));
mutrate=a3.Variables;
threshold_mt=0.03;

DMSfile=cell(size(transcriptids_qualified,1),1);
for i=1:size(transcriptids_qualified,1)
i
   w=find(ismember(transcriptids,transcriptids_qualified(i)));
   dmsdata=[utrlen(w),pos(w),mutrate(w)];
   dmsoutput=zeros(utrlen(w(1)),2);
   dmsoutput(:,1)=1:utrlen(w,1);
   wpass=find(dmsdata(:,3)>=threshold_mt);
   dmsoutput(dmsdata(wpass,2),2)=1;
   DMSfile{i,1}=dmsoutput;
   DMSfile{i,2}=fiveutrhead{w1(i)};
   DMSfile{i,3}=fiveutrseq{w1(i)};

end



for i=1:size(DMSfile,1)
    writematrix(DMSfile{i,1},strcat(DMSfile{i,2},'.txt'),'Delimiter',' ')
end

for i=1:size(DMSfile,1)
    fastawrite(strcat(DMSfile{i,2},'.fa'),DMSfile{i,2},DMSfile{i,3})
end

kong=' ';
mingling=cell(size(DMSfile,1),1);
for i=1:size(DMSfile,1)
    s=append('/home/zdh1996/software/ViennaRNA-2.4.13/src/bin/RNAfold -p -d2 --noLP --shape=',DMSfile{i,2},'.txt -i',kong,DMSfile{i,2},'.fa >',kong,DMSfile{i,2},'.constraint.out');
    mingling{i,1}=s;

end

writetable(mingling,'runhumanconstrantRNAfold')







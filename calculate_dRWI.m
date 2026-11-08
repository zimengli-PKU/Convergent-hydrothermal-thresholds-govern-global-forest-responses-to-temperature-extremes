%read the data
cd 'C:'\Users\18000\Desktop\
[num,txt] = xlsread('\RWI_info.xlsx');
source=length(num);
set(0,'defaultfigurecolor','w')
%% lat&lon-map
tm0=nan(source,12*122);
ncfile =('G:\data and documents\cru_ts4.07.1901.2022.tmp.dat.nc');
tmp=ncread(ncfile,'tmp');
col=ceil((num(:,2)+90)/0.5);
row=ceil((num(:,3)+180)/0.5);
for i = 1:source
  tm0(i,:)=tmp(row(i),col(i),1:12*122);
  if numel(~isnan(tm0(i,:)))==0
        A=tmp(row(i)-2:row(i)+2,col(i)-2:col(i)+2,1:12*122);
        tm0(i,:) = mean(mean(A, 1, 'omitnan'), 2, 'omitnan');
  end
end
clear ncfile;

 load('squeeze_GS_mon.mat');
 som=nan(source,1);eom=nan(source,1);
 for i=1:source
 som(i)=mon(360-col(i),row(i),1);
 eom(i)=mon(360-col(i),row(i),2);
 if isnan(som(i))
 som(i)=round(tsnanmean(mon(360-col(i)-2:360-col(i)+2,row(i)-2:row(i)+2,1),'all'));
 end
 if isnan(eom(i))
 eom(i)=round(tsnanmean(mon(360-col(i)-2:360-col(i)+2,row(i)-2:row(i)+2,2),'all'));
 end

 end
 clear mon;
 gs_len=nan(source,1);
 for i=1:source
    if som(i)<=eom(i)%north?
        gs_len(i)=eom(i)-som(i)+1;
    else
        gs_len(i)=eom(i)+12-som(i)+1;
    end
 end
tm_gs=nan(source,122);

for i=1:source
    if ~isnan(som(i))
    for j=1:122
        if 12*(j-1)+som(i)+gs_len(i)-1<=12*122
        tm_gs(i,j)=mean(tm0(i,(12*(j-1)+som(i)):(12*(j-1)+som(i)+gs_len(i)-1)));
    
        end
    end
    end
end
disp('GST done');

%% extreme identification
smooth_tm=nan(source,122);
for i=1:source
    if numel(find(~isnan(tm_gs(i,:))))>0
smooth_tm(i,:)=smooth(tm_gs(i,:),10)';
    end
end

tm_gs_detrended=tm_gs-smooth_tm;
tm_gs_detrended=tm_gs_detrended(:,60:122);
tm_gs_de_jiduan=nan(source,63);
tm_gs_std=nan(source,1);
tm_gs_std_value=nan(source,63);
for i=1:source
    seq=tm_gs_detrended(i,:);
    seq=seq(~isnan(seq));
        tm_gs_std(i)=std(seq);
        findp=find(tm_gs_detrended(i,:)>1.5*tm_gs_std(i));
        findn=find(tm_gs_detrended(i,:)<-1.5*tm_gs_std(i));
        tm_gs_de_jiduan(i,findp)=1;
        tm_gs_de_jiduan(i,findn)=-1;
        std_tmp=tm_gs_detrended(i,:)/tm_gs_std(i);
        tm_gs_std_value(i,findp)=std_tmp(findp);
        tm_gs_std_value(i,findn)=std_tmp(findn);      
end
clear bp finp findn std_tmp;
% percentage of events
perc=nan(source,2);
for i=1:source
for k=1:2
    perc(i,k)=numel(find(tm_gs_de_jiduan(i,:)==(3-2*k)))/63;
end
end
perc(perc==0)=nan;
tsnanmean(perc(:,1))
tsnanmean(perc(:,2))
% independent
extre_tmp=nan(source,63);
for i=1:source
x=tm_gs_de_jiduan(i,:);
findjiduan=find(~isnan(x));
finddep=find(diff(findjiduan)==1);
x(findjiduan(finddep+1))=nan;
extre_tmp(i,:)=x;
end
clear x;
% select severity for independent events
tm_gs_std_value(isnan(extre_tmp))=nan;

% RWI
rwi=nan(source,63);

cd('C:\Users\18000\Desktop\RWI\RWI0\');
for i=1:source
dataname=char(txt(i+1,1));
A=readtable([lower(dataname),'.csv']);
A=A(table2array(A(:,4)>0),:);
yr=table2array(A(:,1));
b=table2array(A(:,2));
if ~isnumeric(b)
    b=cellfun(@str2num,b);
end

for j=1:63
    if ~isempty(find(yr==j+1959,1))
rwi(i,j)=b(find(yr==j+1959,1));
    end
end

end
clear table dataname title;


r_mean=nan(source,1);
for i=1:source
r_mean(i)=tsnanmean(rwi(i,isnan(tm_gs_de_jiduan(i,:))));
end

r_hc=nan(source,63,2);
for i=1:source
r_hc(i,extre_tmp(i,:)==1,1)=rwi(i,extre_tmp(i,:)==1);
r_hc(i,extre_tmp(i,:)==-1,2)=rwi(i,extre_tmp(i,:)==-1);
end

% influence//resistance
resistance=r_hc;
for i=1:source
resistance(i,:,:)=resistance(i,:,:)/r_mean(i)-1;
end

Rt=squeeze(tsnanmean(resistance,2));
disp('resistance&k_for_sd done');

sd_hc=nan(source,63,2);
for i=1:source
sd_hc(i,extre_tmp(i,:)==1,1)=tm_gs_std_value(i,extre_tmp(i,:)==1);
sd_hc(i,extre_tmp(i,:)==-1,2)=tm_gs_std_value(i,extre_tmp(i,:)==-1);
end
sd0=squeeze(tsnanmean(sd_hc,2));

col={'lat','lon','hotr','coldr','hotsd','coldsd'}; 
result_table=table(num(:,2),num(:,3),Rt(:,1),Rt(:,2),sd0(:,1),sd0(:,2),'VariableNames',col);
writetable(result_table,'C:\Users\18000\Desktop\response.csv');
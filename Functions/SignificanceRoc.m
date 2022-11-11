%ROC for binaryfeatures (sex) the testimony of the credibility of RIP_S
%input RIP_S True_sex

function[auc]=SignificanceRoc_predictNose(RIP_S,True)
data=[RIP_S True];
data=sortrows(data,1);
RIP_test=zeros(size(RIP_S,1),1);
Roc_value=zeros(size(data,1),2);
auc=0;
for i=1:size(data,1)
    index_m=find(RIP_S>data(i,1));
    index_f=find(RIP_S<=data(i,1));
    RIP_test(index_m)=1;
    RIP_test(index_f)=0;
    %confound matrix row(observe male female) col(predict male female)
    confound_matrix=zeros(2,2);
    confound_matrix(1,1)=sum(RIP_test==True & RIP_test==1);
    confound_matrix(1,2)=sum(RIP_test~=True & RIP_test==1);
    confound_matrix(2,1)=sum(RIP_test~=True & RIP_test==0);
    confound_matrix(2,2)=sum(RIP_test==True & RIP_test==0);
    FPR=confound_matrix(1,2)/sum(confound_matrix(:,2));
    TPR=confound_matrix(1,1)/sum(confound_matrix(:,1));
    Roc_value(i,1)=FPR;
    Roc_value(i,2)=TPR;
end
dif_FRP=abs(diff(Roc_value(:,1)));
auc=sum(dif_FRP.*Roc_value([1:end-1],2));
%plot
plot(Roc_value(:,1),Roc_value(:,2),'LineWidth',2,'Color',[0 0.4470 0.7410]);
h=area(Roc_value(:,1),Roc_value(:,2));
h.FaceColor = [1 0.98 0.94];
title('Receiver Operating Characteristic Curve');
set(gca, 'Fontname', 'Times New Roman','FontSize',14); %Font style
xlabel('False Positive Rate','fontsize',14);
ylabel('True Positive Rate','fontsize',14);
text(0.5,0.5,['Area Under Curve=',num2str(auc)],'horiz','center','color',[0.7,0.13,0.13],'FontSize',14); 

end
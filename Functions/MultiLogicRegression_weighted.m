function [y_test]=MultiLogicRegression_weighted(SNP_tag_noNA,NumType,NumTestingSample)
%% regression
% Multivariate logistic regression-Refer to Andrew Wu 
% Reference-https://blog.csdn.net/jizhidexiaoming/article/details/80386360
%Here use noNA data
%% Data Format
%SNP_tag_noNA[x1,x2,x3```,y] 
%limitation of y:only one column 
%% 1-load data
y=zeros(size(SNP_tag_noNA,1),NumType); 
%1,2,3列分别是标签为1,2,3类
for i=1:NumType
    y(find(SNP_tag_noNA(:,size(SNP_tag_noNA,2))==i),i)=1;
end
x=SNP_tag_noNA(:,1:size(SNP_tag_noNA,2)-1);
%% 2-Parameters set
% % % iteration = 10000;
% % % sample_num = length(x); % 样本个数
% % % x = [ones(sample_num, 1), x];
% % % theta = zeros(size(x, 2), 1); % 参数
% % % alpha = 0.1;
% % % %% 3-特征归一化
% % % for i=2:42   %exclude the first 1 column
% % %     x(:,i) = (x(:,i)- mean(x(:,i)))./ std(x(:,i));
% % % end
% % % %% 4-iteration
% % % %暂时将y变为n*1
% % % y=y(:,3);
% % % 
% % % for i = 1:iteration
% % %     h = 1 ./ (1 + exp(-x * theta)); % 通过假设函数得到预测值
% % %     J(i,1) = -1/sample_num * (y' * log(h+eps) + (1-y)'*log(1-h+eps)); % 当前参数下的损失值  %%%%%%% +eps FOR WHAT
% % %     for j=1:length(theta)
% % %         theta(j,1) = theta(j,1) - alpha * sum((h - y) .* x(:,j));  % 更新参数
% % %     %theta = theta - alpha * x'*(h-y); % 同时更新所有参数
% % %     end
% % % end

%% 5 综合上述步骤＋validation(training/testing) 
%300 testing samples
test_index=randperm(size(x,1),NumTestingSample);
train_index=setdiff([1:size(x,1)],test_index);
y_test=y(test_index,:);
x_test=x(test_index,:);
y_train=y(train_index,:);
x_train=x(train_index,:);

iteration = 10000;
sample_num = length(x_train); % 样本个数
x_train = [ones(sample_num, 1), x_train];
x_test = [ones(length(x_test), 1), x_test];
alpha = 0.1;

for i=2:size(x_train,2)   %exclude the first 1 column
    x_train(:,i) = (x_train(:,i)- mean(x_train(:,i)))./ std(x_train(:,i));
    x_test(:,i) = (x_test(:,i)- mean(x_test(:,i)))./ std(x_test(:,i));
end

for k=1:NumType
    y_train_1=y_train(:,k); 
    theta = zeros(size(x_train, 2), 1); % 参数
    for i = 1:iteration
        h = 1 ./ (1 + exp(-x_train * theta)); % 通过假设函数得到预测值
        %J(i,1) = -1/sample_num * (y_train_1' * log(h+eps) + (1-y_train_1)'*log(1-h+eps)); % 当前参数下的损失值  %%%%%%% +eps FOR WHAT
        %加权
        tmp=h - y_train_1;
        tmp(find(y_train_1==1))=tmp(find(y_train_1==1))*sum(y_train_1==0)/sample_num;
        tmp(find(y_train_1==0))=tmp(find(y_train_1==0))*sum(y_train_1==1)/sample_num;
        for j=1:length(theta)
            theta(j,1) = theta(j,1) - alpha * sum(tmp .* x_train(:,j));  % 更新参数
        %theta = theta - alpha * x'*(h-y); % 同时更新所有参数
        end
    end
    y_pre=1 ./ (1 + exp(-x_test * theta));
    y_test=[y_test,y_pre];
end
end
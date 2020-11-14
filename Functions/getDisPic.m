function getDisPic()
% % lm=zeros(27,3);
% % lm(:,1)=1:27;
% % lm(1,[2:3])=[173,300];
% % lm(2,[2:3])=[173,304];
% % lm(3,[2:3])=[143,185];
% % lm(4,[2:3])=[203,185];
% % lm(5,[2:3])=[132,288];
% % lm(6,[2:3])=[215,288];
% % lm(7,[2:3])=[99,187];
% % lm(8,[2:3])=[173,347];
% % lm(9,[2:3])=[247,187];
% % lm(10,[2:3])=[173,253];
% % lm(11,[2:3])=[173,326];
% % lm(12,[2:3])=[173,288];
% % lm(13,[2:3])=[101,221];
% % lm(14,[2:3])=[247,221];
% % lm(15,[2:3])=[102,149];
% % lm(16,[2:3])=[248,149];
% % lm(17,[2:3])=[173,163];
% % lm(18,[2:3])=[173,236];
% % lm(19,[2:3])=[149,249];
% % lm(20,[2:3])=[197,249];
% % lm(21,[2:3])=[173,277];
% % lm(22,[2:3])=[164,273];
% % lm(23,[2:3])=[182,273];
% % lm(24,[2:3])=[167,252];
% % lm(25,[2:3])=[179,252];
% % lm(26,[2:3])=[118,175];
% % lm(27,[2:3])=[227,175];
% % fid=fopen(['D:\guolu\features\CartoonpicLm.txt'], 'wt');
% % for i = 1 : size(lm,1)
% % 	fprintf(fid, '%g\t%g\t%g\t', lm(i,:));
% % 	fprintf(fid, '\n');
% % end
PointsCombinations=combntns([1:27],2); % combination of all the distance conditions
hf=imshow('D:\guolu\pic_ori.jpg');
hf=gca;
hold on;
lm=importdata('D:\guolu\features\CartoonpicLm.txt');
colorset= hsv(351);
for i=1:size(PointsCombinations,1)
    x1=lm(PointsCombinations(i,1),2);
    x2=lm(PointsCombinations(i,2),2);
    y1=lm(PointsCombinations(i,1),3);
    y2=lm(PointsCombinations(i,2),3);
    %line([x1,x2],[y1,y2],'color',[0.8,0.36,0.36]);
    line([x1,x2],[y1,y2],'color',colorset(i,:));

    hold on;
end
scatter(lm(:,2),lm(:,3),'MarkerEdgeColor',[0.7,0.13,0.13],...
              'MarkerFaceColor',[0.7,0.13,0.13],...
              'LineWidth',1);
hold on;
saveas(gca,['D:\guolu\features\pic_dist.jpeg']);





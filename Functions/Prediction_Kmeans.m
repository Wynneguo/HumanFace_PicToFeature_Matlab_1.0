function Prediction_Kmeans(data,N_class)
pattern=data;
[m,n]=size(pattern);
center=zeros(N_class,n);
for x = 1 : N_class
    center(x,:) = pattern(randi(m,1),:); % 第一次随机产生聚类中心 randi返回N_class个从m中随机选择的数
end
while true
distance = zeros(1,N_class);   % 产生1行N列的零矩阵
num = zeros(1,N_class);        % 产生1行N列的零矩阵-保存各类数量
new_center = zeros(N_class,n); % 产生N行n列的零矩阵
%% 将所有的点打上标签1 2 3...N
for x = 1 : m
    for y = 1 : N_class
        distance(y) = norm(pattern(x,1:n) - center(y,:)); % norm函数计算到每个类的距离
    end
    [~,temp] = min(distance); %求最小的距离 ~是距离值，temp是第几个
    pattern(x,n + 1) = temp;   % 增加一列-打标签-与N_class中的某分类tmp距离最近，则其标签为tmp
end
k = 0;
%% 将所有在同一类里的点坐标全部相加，计算新的中心坐标
for y = 1 : N_class
    for x = 1 : m
        if pattern(x,n + 1) == y
           new_center(y,:) = new_center(y,:) + pattern(x,1:n);
           num(y) = num(y) + 1;
        end
    end
    new_center(y,:) = new_center(y,:) / num(y);
    if norm(new_center(y,:) - center(y,:)) < 0.1 %新旧差距不大-k记为确认是中心
        k = k + 1;
    end
end
if k == N_class %即class数未增加
     break;
else
     center = new_center;
end
end
[m, n] = size(pattern); 
%% visualization
figure;
hold on;
for i = 1 : m
    if pattern(i,n) == 1 
         plot(pattern(i,1),pattern(i,2),'r*');
         plot(center(1,1),center(1,2),'ko');
    elseif pattern(i,n) == 2
         plot(pattern(i,1),pattern(i,2),'g*');
         plot(center(2,1),center(2,2),'ko');
    elseif pattern(i,n) == 3
         plot(pattern(i,1),pattern(i,2),'b*');
         plot(center(3,1),center(3,2),'ko');
    elseif pattern(i,n) == 4
         plot(pattern(i,1),pattern(i,2),'y*');
         plot(center(4,1),center(4,2),'ko');
    else
         plot(pattern(i,1),pattern(i,2),'m*');
         plot(center(5,1),center(5,2),'ko');
    end
end
function softmax_multiclass_grad_cw()


% load in data
[X,y] = load_data();


% initializations
N = size(X,2);
C = length(unique(y));
X = [ones(size(X,1),1), X]';
W0 = randn(N+1,C);
alpha = 0.1;
P=40;
x_new=[0.5,0.5];
X_new=[1,x_new]';
% find separators via multiclass softmax classifier
W = softmax_multiclass_grad(X,y,W0,alpha);
accuracy=1-accu(X,y,P,W,C)/P;
accuracy
predict_new=Class(X_new,W,C);
predict_new
% plot the separators as well as final classification
plot_separators(W, X, y)


%%%%%%%%%%%%%%%%% subfunctions %%%%%%%%%%%%%%%%%%%
function W = softmax_multiclass_grad(X,y,W0,alpha)
  
    % initialize
    max_its = 10000; 
    [N,P] = size(X);
    C = length(unique(y));
    W = W0;
    k = 1;
    %grad=zeros(3,4);
    
    %%% main %%%
    while k <= max_its
    %while norm(grad) >= 1*10^(-7)
%%% fill in your codes below
        grad=zeros(3,4);
        for c=1:1:C
            for p=1:1:P
                sum_sub=subsum(p,c,W);
                b=belongs(p,c);
                grad(1:end,c)=grad(1:end,c)+(1/sum_sub-b)*X(1:end,p);
            end
        end
%%% fill in your codes above

        W = W - alpha*grad;
        
        % update counter
        k = k + 1;
        
    end
end

function [X,y] = load_data()
    data = csvread('4class_data.csv');
    X = data(:,1:end - 1);
    y = data(:,end);       
end

function plot_separators(W,X,y,deg)
    
    red =  [ 1 0 .4];
    blue =  [ 0 .4 1];
    green = [0 1 0.5];
    cyan = [1 0.7 0.5];
    grey = [.7 .6 .5];
    colors = [red;blue;green;cyan;grey]; 
    % plot data
    subplot(1,3,1)
    plot_data(X,y)
    title('Distribution of points')
    subplot(1,3,2)
    title('Linear separators')
    plot_data(X,y)
    subplot(1,3,3)
    plot_data(X,y)
    title('Classification')

    %%% plot all linear separators %%%
    subplot(1,3,2)
    num_classes = length(unique(y));
    x = [0:0.01:1];
    for j = 1:num_classes
        hold on
        w = W(:,j);    
        plot (x,(-w(1)-w(2)*x)/w(3),'Color',colors(j,:),'linewidth',2);
    end

    %%% generate max-separator surface %%%
    s = [0:0.005:1];
    [s1,s2] = meshgrid(s,s);
    s1 = reshape(s1,numel(s1),1);
    s2 = reshape(s2,numel(s2),1);
    
    % compute criteria for each point in the range [0,1] x [0,1]
    square = [s1(:), s2(:)];
    p = [ones(size(s1(:),1),1),s1(:), s2(:)];
    f = W'*p';         
    [f,z] = max(f,[],1);
    
    % fill in appropriate regions with class colors
    subplot(1,3,3)
    ind = find(z == 1);
    k = boundary(square(ind,:));
    v = [square(ind(k),:), ones(length(k),1)];
    f = 1:length(k);
    patch('Faces',f,'Vertices',v,'FaceColor',red,'FaceAlpha',0.3)

    ind = find(z == 2);
    k = boundary(square(ind,:));
    v = [square(ind(k),:), 2*ones(length(k),1)];
    f = 1:length(k);
    patch('Faces',f,'Vertices',v,'FaceColor',blue,'FaceAlpha',0.3)
    
    ind = find(z == 3);
    k = boundary(square(ind,:));
    v = [square(ind(k),:), 3*ones(length(k),1)];
    f = 1:length(k);
    patch('Faces',f,'Vertices',v,'FaceColor',green,'FaceAlpha',0.3)
    
    ind = find(z == 4);
    k = boundary(square(ind,:));
    v = [square(ind(k),:), 4*ones(length(k),1)];
    f = 1:length(k);
    patch('Faces',f,'Vertices',v,'FaceColor',cyan,'FaceAlpha',0.3)
    
    ind = find(z == 5);
    k = boundary(square(ind,:));
    v = [square(ind(k),:), 4*ones(length(k),1)];
    f = 1:length(k);
    patch('Faces',f,'Vertices',v,'FaceColor',grey,'FaceAlpha',0.3)
    
    % produce decision boundary
    s1 = reshape(s1,[length(s),length(s)]);
    s2 = reshape(s2,[length(s),length(s)]);   
    z = reshape(z,[length(s),length(s)]);   
    num_classes = length(unique(z));
    subplot(1,3,3)
    for i = 1:num_classes - 1
       hold on
       contour(s1,s2,z,[i + 0.5,i + 0.5],'Color','k','LineWidth',2)
    end
    
    % make plot real nice lookin'
    for i = 1:3
        subplot(1,3,i)
        axis([0 1 0 1])
        axis square
        xlabel('x_1','FontName','cmmi9','Fontsize',18)
        ylabel('x_2','FontName','cmmi9','Fontsize',18)
        set(get(gca,'YLabel'),'Rotation',0)
        zlabel('y','FontName','cmmi9','Fontsize',18)
        set(get(gca,'ZLabel'),'Rotation',0)

        set(gca,'XTick',[0,1])
        set(gca,'YTick',[0,1])
        set(gca,'ZTick',[0:1:num_classes])
        set(gcf,'color','w');
    end
end

function plot_data(X,y)
    red =  [ 1 0 .4];
    blue =  [ 0 .4 1];
    green = [0 1 0.5];
    cyan = [1 0.7 0.5];
    grey = [.7 .6 .5];
    colors = [red;blue;green;cyan;grey]; 
    % how many classes in the data? maximum 4 here.
    class_labels = unique(y);           % class labels
    num_classes = length(class_labels);

    % plot data
    for i = 1:num_classes
        class = class_labels(i);
        ind = find(y == class);
        hold on
        scatter3(X(2,ind),X(3,ind),class*ones(length(ind),1),'Linewidth',2,'Markeredgecolor',colors(i,:),'markerFacecolor','none');
        hold on
        scatter3(X(2,ind),X(3,ind),class*ones(length(ind),1),'Linewidth',2,'Markeredgecolor',colors(i,:),'markerFacecolor','none');
    end
    axis([0 1 0 1])
    axis square
    box on
end

function sum_sub=subsum(i,j,W)
    sum_sub = exp(X(1:end,i)'*(W(1:end,1)-W(1:end,j)))+exp(X(1:end,i)'*(W(1:end,2)-W(1:end,j)))...
        +exp(X(1:end,i)'*(W(1:end,3)-W(1:end,j)))+exp(X(1:end,i)'*(W(1:end,4)-W(1:end,j)));
end



function b=belongs(p,c)
    if y(p)==c
        b=1;
    else
        b=0;
    end
end

function idx=Class(x,w,C)
    for c=1:1:C
        yp(1,c)=x'*w(1:end,c);
    end
    [unused_value,idx] = max(yp);
end

function a=accu(x,y,P,w,C)
    a=0;
    for p=1:1:P
        original_class=y(p);
        xt=x(1:end,p);
        predict_class=Class(xt,w,C);
        if predict_class ~= original_class
            a=a+1;
        else
            a=a;
        end
    end
end


end
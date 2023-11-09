function ohm_wrapper_cw()

% load in data
data = csvread('ohms_data.csv');
x = data(:,1);
y = data(:,2);

% fit model in transformed feature space
x_tilde = [ones(length(x),1) x];  % compact notation

%%% fill in your codes below
y_test = [y.^(-1),-sqrt(y),-exp(y),-y.^2,y]; %possible solutions
% calculate correct answer by Correlation Coefficient
coef_matrix = corrcoef([x,y_test]);
coef_matrix
[corr,idx1] = max(coef_matrix(1,2:end));
y_transformed = y_test(1:end,idx1);
w_tilde = linsolve(x_tilde,y_transformed);  % compact notation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_tilde
b = w_tilde(1);
w = w_tilde(2);
b
w
% plot model in the original space
figure(1)
hold on
model = [0:0.01:100];
%%% fill in your codes below
out = (b+ w*model).^-1; 
input_length = 60;
predict = 1/(b+w*input_length);
predict_transformed = 1/predict;
predict
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(model,out,'m','LineWidth',2)
plot(x,y,'o','MarkerEdgeColor','k','MarkerFaceColor','k')
plot(input_length,predict,'o','MarkerSize', 10, 'MarkerEdgeColor','g','MarkerFaceColor','g')
box on
axis tight

% make labels correct
title('Original Data Space','Fontsize',18)
xlabel('length of wire (cm)','Fontsize',14,'FontName','cmr10')
ylabel('current-y (A)','Fontsize',14,'FontName','cmr10')
set(get(gca,'YLabel'),'Rotation',90)
set(gca,'FontSize',12); 
hold off
legend('model','original','predict')
%%% fill in your codes below
% plot model in the transformed space
figure(2)
hold on
model = [0:0.01:100];

out = b+ w*model; 

plot(model,out,'m','LineWidth',2)
plot(x,y.^-1,'o','MarkerEdgeColor','k','MarkerFaceColor','k')
plot(input_length,predict_transformed,'o','MarkerSize', 10, 'MarkerEdgeColor','g','MarkerFaceColor','g')
box on
axis tight

% make labels correct
title('Transformed Data Space','Fontsize',18)
xlabel('length of wire (cm)','Fontsize',14,'FontName','cmr10')
ylabel('f(y)=1/y','Fontsize',14,'FontName','cmr10')
set(get(gca,'YLabel'),'Rotation',90)
set(gca,'FontSize',12); 
legend('model','original','predict','Location','northwest')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end






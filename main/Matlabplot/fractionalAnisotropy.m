I = eye(3);
ef = zeros(3);
ef(1,1) = cos(0);
ef(1,2) = sin(0);
theta = 0 ; %linspace(0, pi/2, 6);
kappa = linspace(0,0.33, 20);
FA = zeros(1, length(kappa));
for j = 1:length(theta)
    ef(1,1) = cos(theta(j));
    ef(1,2) = sin(theta(j));
    for i =1:length(kappa)
        H = kappa(i)*I + (1 - 3*kappa(i))*ef;
        lambda = sort(eig(H), 'descend');
        numerator = (lambda(1) -lambda(2))^2 + (lambda(2) -lambda(3))^2 + (lambda(3) -lambda(1))^2;
        denominator = (lambda(1))^2 + (lambda(2))^2 + (lambda(3))^2;
        FA(i) = sqrt(numerator/ (denominator*2));
    end
    figure(1)
    plot(kappa,FA, 'k', 'linewidth', 3 )
    hold on
end

% figure(1)
% plot(kappa,FA, 'k', 'linewidth', 1.5)
grid on
grid minor 
% legend('J_0','J_1','J_2','J_3','J_4','Location','Best')
%title('Relation between fractional anisotropy and disarray','interpreter','latex', 'FontSize',14)
%xlabel('kappa','interpreter','latex','fontweight','bold', 'FontSize',16)
%ylabel('Fractional Anisotropy','interpreter','latex','fontweight','bold', 'FontSize',16)
%title('Relationship','fontweight','bold', 'FontSize',16)
xticks([0 0.11 0.22 0.34])
yticks([0 0.25 0.5 0.75 1.0 1.25])
%xticklabels({}, 'fontweight','bold')
xlabel('Myofiber Disarray','fontweight','bold', 'FontSize',16)
ylabel('Fractional Anisotropy','fontweight','bold', 'FontSize',16)
saveas(gcf, 'Fractional Anisotropy & Kappa', 'png')
%figure
subplot(3,1,1)
semilogy(gas.General.cev_ds.*1000, table.N2.et.*100, 'b+')
line([0 15], [2 2], 'color', 'black');
ylim([1 100])
xlabel('CEV [l]')
ylabel('Cet [%]')

subplot(3,1,2)
semilogy(gas.General.cev_ds.*1000, table.N2.et.*100, 'b+')
line([gas.General.cev_ds(12) gas.General.cev_ds(13)].*1000, [table.N2.et(12) table.N2.et(13)].*100)
line([0 15], [2 2], 'color', 'black');
ylim([1 100])
xlabel('CEV [l]')
ylabel('Cet [%]')

subplot(3,1,3)
hold off
semilogy(gas.General.cev_ds.*1000, table.N2.et.*100, 'b+')
handle = @(x)80.*exp(-0.5.*x)+1.0.*exp(-0.001.*x);
a=0:0.01:15;
hold on
semilogy(a, handle(a), 'r-')
line([0 15], [2 2], 'color', 'black');
ylim([1 100])
xlabel('CEV [l]')
ylabel('Cet [%]')

% %Stop at the end of the calculatefastslow function to let this run
% 
% x=linspace(xNum(1), xNum(end));
% xelin = linspace(10, 32);
% xs = 1:10;
% xe = 10:32;
% subplot(2,1,1);
% hold off
% semilogy(xNum,yCet, 'ro')
% hold on
% semilogy(x,expFuncitonsC(7).function(expFuncitonsC(7).params,x), 'b-')
% semilogy(x,expFuncitonsC(7).fFast(expFuncitonsC(7).pFast,x), 'r-')
% semilogy(x,expFuncitonsC(7).fSlow(expFuncitonsC(7).pSlow,x), 'g-')
% line([0 35], [yCet(1)./40 yCet(1)./40], 'lineStyle', ':', 'color', 'black', 'LineWidth', 1)
% legend('N_2 end tidal concentration', 'Modelled concentration decay', 'Contribution of the fast compartment', 'Contribution of the slow compartment', 'Target concentration')
% ylabel('N_2 end tidal concentration [%]')
% xlabel('Breath number')
% set(gca, 'YTick', [0.01 0.02 0.03 0.05 0.10 0.20 0.28 0.28 0.50 0.80])
% set(gca, 'YTickLabel', [1 2 3 5 10 20 30 40 50 80])
% ylim([0.008 0.85])
% grid on
% title('Modelled completed measurement')
% 
% subplot(2,1,2);
% hold off
% semilogy(xNum(xs),yCet(xs), 'ro')
% hold on
% area(xNum(xs), yCet(xs), 'FaceColor', [0.7 0.7 0.99])
% area(xelin, expFuncitonsC(7).function(expFuncitonsC(7).params,xelin), 'FaceColor', [0.99 0.7 0.7])
% semilogy(xNum(xs),yCet(xs), 'ro')
% semilogy(x,expFuncitonsC(7).function(expFuncitonsC(7).params,x), 'b-')
% semilogy(x,expFuncitonsC(7).fFast(expFuncitonsC(7).pFast,x), 'r-')
% semilogy(x,expFuncitonsC(7).fSlow(expFuncitonsC(7).pSlow,x), 'g-')
% line([0 35], [yCet(1)./40 yCet(1)./40], 'lineStyle', ':', 'color', 'black', 'LineWidth', 1)
% legend('N_2 end tidal concentration', 'Measured values', 'Extrapolated values')
% ylabel('N_2 end tidal concentration [%]')
% xlabel('Breath number')
% set(gca, 'YTick', [0.01 0.02 0.03 0.05 0.10 0.20 0.28 0.28 0.50 0.80])
% set(gca, 'YTickLabel', [1 2 3 5 10 20 30 40 50 80])
% ylim([0.008 0.85])
% grid on
% 
% title('Predicted measurement')
% 
% % x=1:7;
% % area(xNum(x), yCet(x))
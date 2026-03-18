function Result_plot( save_sketch, save_result, DEL, BEST_D_G, BEST_NUM)
global tar
line= tar.line;

L1 = save_sketch.line1;
L2 = save_sketch.line2;
L3 = save_sketch.line3;
L4 = save_sketch.line4;
LC = save_sketch.center;
tarline = save_result.tardata;
resline = save_result.resdata;


fprintf('Best DEL=:%g    ',BEST_D_G);
fprintf('最佳个体信息 ： JOB-%d-%d \n', BEST_NUM(1), BEST_NUM(2));
%绘制适应度曲线
figure('color','w')
plot   (DEL,'r')
xlabel ('迭代次数')
ylabel ('误差值')
title  ('适应度进化曲线')
set(gca,'yticklabel',get(gca,'ytick'));

%绘制二维草图
figure('color','w')
axis equal
for k = 1:line
    plot(L1{k}(:,1),L1{k}(:,2),'b');
    hold on
    
    plot(L2{k}(:,1),L2{k}(:,2),'b');
    plot(L3{k}(:,1),L3{k}(:,2),'b');
    plot(L4{k}(:,1),L4{k}(:,2),'b');
    plot(LC{k}(:,1),LC{k}(:,2),'--r');
    plot([L3{k}(1,1),L3{k}(end,1)],[L3{k}(1,2),L3{k}(end,2)],'r')
    plot([L4{k}(1,1),L4{k}(end,1)],[L4{k}(1,2),L4{k}(end,2)],'r')
end
axis equal
hold off

%绘制三维中心线
figure('color','w')
axis equal

for k = 1:line
    plot3( tarline{k,1}, tarline{k,2}, tarline{k,3},'r');
    hold on
    plot3( resline{k,1}, resline{k,2}, resline{k,3},'b');
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
end
axis equal
legend('tar','res')
hold off

end


function plotScenario(obstacles,d_horiz,d_vert,a_i)
    hold on
    plot([0 d_horiz],[d_vert d_vert],'Color','b','LineWidth', 1.5);
    plot([d_horiz d_horiz],[0 d_vert],'Color','b','LineWidth', 1.5);
    plot([0 0],[0 d_vert],'Color','b','LineWidth', 1.5);
    plot([0 d_horiz],[0 0],'Color','b','LineWidth', 1.5);
    plot(a_i(1,:), a_i(2,:), 'ks', 'LineWidth', 2.5, 'MarkerSize', 20)
    for i = 1 : size(obstacles,3)
        plot(obstacles(:,1,i), obstacles(:,2,i),'Color','r','LineWidth', 1.5);
    end
end

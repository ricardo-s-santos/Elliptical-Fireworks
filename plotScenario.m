function plotScenario(B,a_i)
    hold on
    plot([0 B],[B B],'Color','b','LineWidth', 1.5);
    plot([B B],[0 B],'Color','b','LineWidth', 1.5);
    plot([0 0],[0 B],'Color','b','LineWidth', 1.5);
    plot([0 B],[0 0],'Color','b','LineWidth', 1.5);
    plot(a_i(1,:), a_i(2,:), 'ks', 'LineWidth', 2.5, 'MarkerSize', 20)
end
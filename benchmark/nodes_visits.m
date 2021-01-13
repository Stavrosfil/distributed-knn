num_of_processors = [1 2 4 8 16 32 48 60];

figure;
plot(num_of_processors, nv_corel_hist, '-o', 'MarkerFaceColor', 'black', 'linewidth', 1.5);
hold on; 
plot(num_of_processors, nv_corel_mom, '-o', 'MarkerFaceColor', 'black', 'linewidth', 1.5);
hold on;
plot(num_of_processors, nv_corel_tex, '-o', 'MarkerFaceColor', 'black', 'linewidth', 1.5);
xlabel('number of processors');
ylabel('nodes visits');
xticks(num_of_processors);
title('average nodes visits per query point (vpt search)');
legend('corel colors hist', 'corel colors mom', 'corel cooc tex');
grid on;
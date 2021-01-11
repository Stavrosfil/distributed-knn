figure;
speedup(corelcolorhistogram_v1);
hold on;
speedup(corelcolormoments_v1);
hold on;
speedup(corelcooctexture_v1);
hold on;
speedup(corelcolorhistogram_v2);
hold on;
speedup(corelcolormoments_v2);
hold on;
speedup(corelcooctexture_v2);
legend('corel colors hist v1', 'corel colors mom v1', 'corel cooc tex v1', 'corel colors hist v2', 'corel colors mom v2', 'corel cooc tex v2');


figure;
v1_versus_v2(corelcolorhistogram_v1, corelcolorhistogram_v2);
hold on;
v1_versus_v2(corelcolormoments_v1, corelcolormoments_v2);
hold on;
v1_versus_v2(corelcooctexture_v1, corelcooctexture_v2);
legend('corel colors hist', 'corel colors mom', 'corel cooc tex');
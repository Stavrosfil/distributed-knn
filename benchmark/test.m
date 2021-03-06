% corel

num_of_processors = [1 2 4 8 16 32 48 60];

figure;
speedup(corelcolorhistogram_v1, num_of_processors);
hold on;
speedup(corelcolormoments_v1, num_of_processors);
hold on;
speedup(corelcooctexture_v1, num_of_processors);
legend('corel color hist v1', 'corel color mom v1', 'corel cooc tex v1');


figure
speedup(corelcolorhistogram_v2, num_of_processors);
hold on;
speedup(corelcolormoments_v2, num_of_processors);
hold on;
speedup(corelcooctexture_v2, num_of_processors);
legend('corel color hist v2', 'corel color mom v2', 'corel cooc tex v2');


figure;
v1_versus_v2(corelcolorhistogram_v1, corelcolorhistogram_v2, num_of_processors);
hold on;
v1_versus_v2(corelcolormoments_v1, corelcolormoments_v2, num_of_processors);
hold on;
v1_versus_v2(corelcooctexture_v1, corelcooctexture_v2, num_of_processors);
legend('corel color hist', 'corel color mom', 'corel cooc tex');

%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% tv_news_comm

num_of_processors = [1 2 4 8 16];

figure;
speedup(bbc_v1, num_of_processors);
hold on;
speedup(cnn_v1, num_of_processors);
hold on;
speedup(cnnibn_v1, num_of_processors);
hold on;
speedup(ndtv_v1, num_of_processors);
hold on;
speedup(timesnow_v1, num_of_processors);
legend('bbc v1', 'cnn v1', 'cnnibn v1', 'ndtv v1', 'timesnow v1')


figure;
speedup(bbc_v2, num_of_processors);
hold on;
speedup(cnn_v2, num_of_processors);
hold on;
speedup(cnnibn_v2, num_of_processors);
speedup(ndtv_v2, num_of_processors);
hold on;
speedup(timesnow_v2, num_of_processors);
hold on;
legend('bbc v2', 'cnn v2', 'cnnibn v2', 'ndtv v2', 'timesnow v2');


figure;
v1_versus_v2(bbc_v1, bbc_v2, num_of_processors);
hold on;
v1_versus_v2(cnn_v1, cnn_v2, num_of_processors);
hold on;
v1_versus_v2(cnnibn_v1, cnnibn_v2, num_of_processors);
hold on;
v1_versus_v2(ndtv_v1, ndtv_v2, num_of_processors);
hold on;
v1_versus_v2(timesnow_v1, timesnow_v2, num_of_processors);
legend('bbc', 'cnn', 'cnnibn', 'ndtv', 'timesnow');


%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% fma features

num_of_processors = [1 2 4 8 16 32 48 60];

figure;
speedup(features_v1, num_of_processors);
hold on;
speedup(miniboonepid_v1, num_of_processors);
legend('fma features v1', 'miniboonepid v1');

figure;
speedup(features_v2, num_of_processors);
hold on;
speedup(miniboonepid_v2, num_of_processors);
legend('fma features v2', 'miniboonepid v2');

figure;
v1_versus_v2(features_v1, features_v2, num_of_processors);
hold on;
v1_versus_v2(miniboonepid_v1, miniboonepid_v2, num_of_processors);
legend('fma features', 'miniboonepid');






figure;
X = categorical({'corelcolorhistogram v1', 'corelcolorhistogram v2'; 'corelcolormoments v1', 'corelcolormoments v2'; 'corelcooctexture v1', 'corelcooctexture v2'});
Y = [corelcolorhistogram_v1(1), corelcolorhistogram_v2(1); corelcolormoments_v1(1), corelcolormoments_v2(1); corelcooctexture_v1(1), corelcooctexture_v2(1)];
bar(X, Y, 2);
ylabel('time (ms)');
title('corel sequential times');

figure;
X = categorical({'miniboonepid v1', 'miniboone v2'; 'fma features v1', 'fma features v2'});
Y = [miniboonepid_v1(1), miniboonepid_v2(1); features_v1(1), features_v2(1)];
bar(X, Y, 2);
ylabel('time (ms)');
title('miniboone pid & fma features sequential times');


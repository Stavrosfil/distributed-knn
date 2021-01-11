   function res = v1_versus_v2(v1, v2)

    num_of_processors = [1 2 4 8 16 32 48 60];
    speedup = zeros(1, length(v1));
    
    for i = 1:length(v1)

        speedup(i) = v1(i) / v2(i);
    end

    plot(num_of_processors, speedup, '-o', 'MarkerFaceColor', 'black', 'linewidth', 1.5);
    xlabel('number of processors');
    ylabel('speedup');
    xticks(num_of_processors);
    
    title('t_v_1 / t_v_2 speedup');
    grid on;

    res = speedup;
end
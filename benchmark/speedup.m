function res = speedup(data, num_of_processors)

    
    speedup = zeros(1, length(data));
    
    for i = 1:length(data)

        speedup(i) = data(1) / data(i);
    end

    plot(num_of_processors, speedup, '-o', 'MarkerFaceColor', 'black', 'linewidth', 1.5);
    xlabel('number of processors');
    ylabel('speedup');
    xticks(num_of_processors);
    title('t_1 / t_p speedup');
    grid on;

    res = speedup;
end
function res = speedup(data)

    num_of_processors = [1 2 4 8 16 32 48 60];
    speedup = zeros(1, length(data));
    
    for i = 1:length(data)

        speedup(i) = data(1) / data(i);
    end

    plot(num_of_processors, speedup, '-o', 'MarkerFaceColor', 'black', 'linewidth', 1.5);
    xlabel('number of processors');
    ylabel('speedup');
    xticks(num_of_processors);
    yticks(0:2:30);
    title('v1 & v2 speedup');
    grid on;

    res = speedup;
end
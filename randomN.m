function[w] = randomN(E,W)
    w = zeros(length(W),1);
    for i = 1:length(W)
        w(i) = randn*sqrt(W(i,i)) + E;
    end
end
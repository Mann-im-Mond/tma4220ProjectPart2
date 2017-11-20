function val = alphaFkt(rodValue)
    alphaAl = 9.4 * 10^(-5);
    alphaCake = 1.5 * 10^(-7);
    if rodValue==1001
        val = alphaAl;
    else
        val = alphaCake;
    end
end
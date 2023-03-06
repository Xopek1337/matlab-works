function [signal] = createSignal(snr, data, modOrder, h, delay, sps)
    dataSym = bi2de(data);
    modData = qammod(dataSym, modOrder, 'UnitAveragePower' , true);

    output = upsample(modData, sps);
    
    signal = conv(h, output);
    varSignal = var(signal);

    hd = fdesign.fracdelay(delay, 'n', 3);
    fir = design(hd, 'lagrange', 'FilterStructure', 'farrowfd');
    
    signal = filter(fir, signal);

    varNoise = varSignal*10^(-snr/10);
    signal = signal+sqrt(varNoise/2)*(randn(size(signal))+1i*randn(size(signal)) );
end

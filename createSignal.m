function [signal] = createSignal(snr, data, modOrder, h, delay, sps)
    dataSym = bi2de(data);
    modData = qammod(dataSym, modOrder, 'UnitAveragePower' , true);

    output = upsample(modData, sps);
    
    signal = conv(h, output);
    signal = awgn(signal, snr, 'measured');
    
    hd = fdesign.fracdelay(delay, 'n', 3);
    fir = design(hd, 'lagrange', 'FilterStructure', 'farrowfd');
    
    if(delay ~= 0)
        signal = filter(fir, signal);
    end
end
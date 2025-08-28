function EqualiseHMatrix(read_dir, write_dir, num_rows, num_cols)
    b = [0.9833149058539373 -1.8342605438232056 0.8613046229420136 % -2.1 dB peak @ 810 Hz Q 0.8
        0.9847343376264435 -1.6080644395882402 0.6800564732437505 % -1.4 dB low shelf @ 2 kHz
        % 0.967072031610277 -1.8417678204094883 0.8824546873377966
        % 0.9956116838874545 -1.9783220035903373 0.982793372740487
        % 0.9606108454641027 -1.2837192925854206 0.42887690419522617
        % 0.9946848248759728 -1.9002527445030875 0.9075639847631554
];
    a = [1 -1.8342605438232056 0.8446195287959508
        1 -1.6031037877693375 0.6697514626890969
        % 1 -1.8417678204094883 0.8495267189480735
        % 1 -1.9782464924773528 0.9784805677409256
        % 1 -1.2837192925854206 0.3894877496593288
        % 1 -1.9002527445030875 0.9022488096391281
        ];
    sos = dsp.SOSFilter(b, a);

    PlotBodeMag(sos);

    for row = 1:num_rows
        for col = 1:num_cols
            [ir, fs] = audioread(read_dir + "H_R"+row+"_S"+col+".wav");
            filtered_output = sos(ir);
            audiowrite(write_dir + "H_R"+row+"_S"+col+".wav", filtered_output, fs, "BitsPerSample", 32);
        end
    end
end

function PlotBodeMag(sos)
    nexttile
    [h, w] = freqz(sos);
    semilogy(20*log(abs(h)),w);
    grid on
    xlim([-9 1]);
    fs = 48000;
    ylim([pi * 20/(fs / 2) pi * 20000/(fs / 2)]);
end
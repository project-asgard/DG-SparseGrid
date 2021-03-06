% generate test data for wavelet transform operator

data_dir = strcat("generated-inputs", "/", "basis", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

out_format = strcat(data_dir, 'transform_blocks_l%i_d%i_%d.dat');
levels = [2 5];
for i = levels
    for d=2:3:5
        [FWMT, blocks] = OperatorTwoScale_wavelet2(d, i);
        for b = 1:length(blocks)
            write_octave_like_output(sprintf(out_format,i,d,b), blocks{b});
        end
    end
end

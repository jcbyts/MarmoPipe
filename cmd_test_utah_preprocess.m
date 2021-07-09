

%%

fname = 'C:\Raw\utah_sample\richardExampleVoltages.dat';
fid = fopen(fname, 'r');

DATA = fread(fid, '*int16');

DATA = (reshape(DATA, 96, []));
%%

figure(1); clf
imagesc(DATA(:,1:1e3))
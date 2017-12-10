file = fopen('TAPE27');
data = textscan(file, '%f %f');
fclose(file);

plot(data{1}, data{2});
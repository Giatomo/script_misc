mean_radius = 15;
f = 2/3;

img = getTitle();
selectWindow(img);
run("32-bit");

run("Clear Results");
run("Set Measurements...", "area mean standard min redirect=None decimal=9");
makeRectangle(0, 0, getWidth(), getHeight());
run("Measure");
sd = getResult("StdDev", 0);
k = -sd/(255-(f*sd));

	
name = "NICK thresholding k="+k;
k = "value=" + k;
mean_radius =  "radius=" + mean_radius;

selectWindow(img);
run("Duplicate...", "title=tmp_mean");
run("Duplicate...", "title=tmp_mean(squared)");

selectWindow("tmp_mean");
run("Mean...", mean_radius);
run("Duplicate...", "title=tmp_squared(mean)");

selectWindow("tmp_squared(mean)");
run("Square");

selectWindow("tmp_mean(squared)");
run("Square");
run("Mean...", mean_radius);


imageCalculator("Subtract create 32-bit", "tmp_mean(squared)","tmp_squared(mean)");
selectWindow("Result of tmp_mean(squared)");
rename("tmp_variance");
selectWindow("tmp_variance");

imageCalculator("Add create 32-bit", "tmp_variance", "tmp_squared(mean)");
selectWindow("Result of tmp_variance");
rename("tmp_temp");
selectWindow("tmp_temp");
run("Square Root");
run("Multiply...", k);
imageCalculator("Add create 32-bit", "tmp_mean", "tmp_temp");
rename("tmp_threshold");
selectWindow("tmp_threshold");
expression = "expression=A>B a=["+img+"] b=[tmp_threshold]";
print(expression);
run("Image Expression Parser (Macro)", "expression=A>B a=["+img+"] b=[tmp_threshold]");
rename(name);

selectWindow(img);
run("Set...", 0);
imageCalculator("Add 32-bit", img, name);
close("tmp_*");
close("Result of tmp_*");
close(name);
// print("\\Clear");


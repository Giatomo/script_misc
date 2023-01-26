#@ File (style="directory") path
print(path);
list = getList("image.titles");
if (list.length==0)
	print("No image windows are open");
else {
	print("Image windows:");
	for (i=0; i<list.length; i++){
    	print("   "+list[i]);
		selectWindow(list[i]);
		saveAs("Tiff", path+"/"+list[i]+".tif");
		close();
    }
}



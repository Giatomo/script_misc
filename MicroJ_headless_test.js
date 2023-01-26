importClass(Packages.com.ducret.microbeJ.BatchJ);
importClass(Packages.ij.IJ);
//importClass(Packages.microbeJ.com.ducret.microbeJ);
out = new BatchJ();



out.loadFile("/home/thomas/Bureau/default.xml");
out.selectFolder("/run/media/thomas/DATA/MICROSCOPY/Timelapses/ZapTmsfGFP_Chromosome/2022-05-25/Phase/", loadSetting = true, silent = true);
print(out);

print("Done");
params = out.getParameters();
print(params.get("DETECTION_BACTERIA"));



IJ.wait(10000);
out.startExperiment(4);
print("______________");
//out.run(0);
// /run/media/thomas/DATA/MICROSCOPY/Timelapses/ZapTmsfGFP_Chromosome/2022-05-25/Phase/
//BatchJ.loadFile("/home/thomas/Bureau/default.xml")

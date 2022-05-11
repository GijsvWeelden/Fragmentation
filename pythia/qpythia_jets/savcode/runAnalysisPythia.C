runAnalysisPythia(Char_t const *indir = "output/events_40_80_50", Char_t const *foutname = "qPythiaJetShape") {
 
  gSystem->AddIncludePath("-I$FASTJET/include"); 
  gSystem->AddIncludePath("-I/project/alice/users/marcovl/fjcontrib-1.017-install/include");

  gSystem->AddDynamicPath("/usr/lib64");
  gSystem->AddDynamicPath("/project/alice/users/marcovl/fjcontrib-1.017-install/lib");

  /*
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCDB");
  //gSystem->Load("libTENDER");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGEMCAL");
  
  gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libsiscone");
  gSystem->Load("libsiscone_spherical");
  gSystem->Load("libfastjetplugins");
  //gSystem->Load("libSISConePlugin");
  //gSystem->Load("libCDFConesPlugin");
  gSystem->Load("libfastjetcontribfragile-1.017");
  gSystem->Load("libPWGJEEMCALJetTasks");
*/
  gSystem->Load("libqpythia");
  gSystem->Load("libEGPythia6");
  gSystem->Load("liblhapdf");
  gSystem->Load("libAliPythia6");


  gROOT->LoadMacro("analysePythiaShapes.cxx+");
  //gROOT->LoadMacro("analysePythiaSoftDrop.cxx+");

  analysePythia(indir, foutname);

}

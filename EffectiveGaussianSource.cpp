#include <iostream>
#include <stdio.h>
#include <string.h>
//#include <omp.h>
#include <complex>

#include "CATS.h"
#include "CATSconstants.h"
#include "CATStools.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_Fitters.h"
#include "DLM_CppTools.h"
#include "DLM_CkDecomposition.h"
#include "CommonAnaFunctions.h"
#include "DLM_Random.h"
#include "DLM_Bessel.h"
#include "DLM_Integration.h"
#include "DLM_WfModel.h"
#include "DLM_Histo.h"
#include "DLM_CkModels.h"
#include "DLM_HistoAnalysis.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;


void EffectiveGaussianXiCock(){

  const double CoreSize = 1.2;

  //DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
  DLM_CleverMcLevyResoTM MagicSource;

  //DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
  MagicSource.InitStability(1,2-1e-6,2+1e-6);
  MagicSource.InitScale(38,0.15,2.0);
  MagicSource.InitRad(257*2,0,64);
  MagicSource.InitType(2);
  ///////////////////

  //for p-Xi, set up the amount of secondaries
  //first for the protons (64.22%)
  MagicSource.SetUpReso(0,0.6422);
  //than for the Xis, here its 0% (we have ONLY primordials)
  MagicSource.SetUpReso(1,0.0);

  //the cut off scale in k*, for which the angular distributions from EPOS
  //are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
  const double k_CutOff = 200;

  //to be used for the NTuple later on
  Float_t k_D;
  Float_t fP1;
  Float_t fP2;
  Float_t fM1;
  Float_t fM2;
  Float_t Tau1;
  Float_t Tau2;
  Float_t AngleRcP1;
  Float_t AngleRcP2;
  Float_t AngleP1P2;
  //random generator dimi style. The input is incompatible with the ROOT random generator,
  //do not mix and match, do not ask me how I know this. Ask Bernie.
  //11 is the seed, you can change that to you favorite number
  DLM_Random RanGen(11);
  //dummies to save random shit
  double RanVal1;
  double RanVal2;
  double RanVal3;

  //open the magic file from dimi with the angular distributions.
  TFile* F_EposDisto_pReso_Xim = new TFile("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles/Source/EposAngularDist/EposDisto_pReso_Xim.root");
  //set up the ntuple, do not change anything unless told so by dimi
  TNtuple* T_EposDisto_pReso_Xim = (TNtuple*)F_EposDisto_pReso_Xim->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_pReso_Xim = T_EposDisto_pReso_Xim->GetEntries();
  T_EposDisto_pReso_Xim->SetBranchAddress("k_D",&k_D);
  T_EposDisto_pReso_Xim->SetBranchAddress("P1",&fP1);
  T_EposDisto_pReso_Xim->SetBranchAddress("P2",&fP2);
  T_EposDisto_pReso_Xim->SetBranchAddress("M1",&fM1);
  T_EposDisto_pReso_Xim->SetBranchAddress("M2",&fM2);
  T_EposDisto_pReso_Xim->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_pReso_Xim->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_pReso_Xim->SetBranchAddress("AngleP1P2",&AngleP1P2);
  //iterate over the ntuple
  for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_Xim; uEntry++){
      //get each entry
      T_EposDisto_pReso_Xim->GetEntry(uEntry);
      //disregard the entry of you are outside the desired k*
      if(k_D>k_CutOff) continue;
      //overwrite the value for the lifetime. This is computed from the
      //stat. hadronization model (Vale) or thermal fist (Max)
      //this is the value for the secondary protons
      Tau1 = 1.65;
      //for primoridials (the Xis) we put 0
      Tau2 = 0;
      //put in the average mass of the resonances (again from SHM or TF)
      //this is the value for protons
      fM1 = 1362;
      //generate a random path length for the propagation of the resonances
      //nothing to change!
      RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
      //adds a single entry into the PDF for the angular distribution to be used
      MagicSource.AddBGT_RP(RanVal1,cos(AngleRcP1));


  }
  delete F_EposDisto_pReso_Xim;
  //if you have resonances contributing to both particles, we need to repeat the above procedure
  //for the prim-reso (AddBGT_PR) and reso-reso (AddBGT_RR) cases

  const unsigned NumSourceBins = 128;
  const double rMin = 0;
  const double rMax = 16;
  TFile* fOutput = new TFile("fOutput.root","recreate");
  TH1F* hSource = new TH1F("hSource","hSource",NumSourceBins,rMin,rMax);

  //fill the histo fro the source
  for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
    //get the x-axis (r value) of the current bin
    double xaxis = hSource->GetBinCenter(uBin+1);
    //an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
    double parameters[2];
    parameters[0] = CoreSize;
    parameters[1] = 2.0;
    double SourceValue = MagicSource.RootEval(&xaxis, parameters);
    hSource->SetBinContent(uBin+1,SourceValue);
    //infinite errors for now
    hSource->SetBinError(uBin+1,1000.);
  }

  //idea: fit the source distribution only in a range around its peak
  //to do this: silly idea: put very large uncertainties in the bins outside of this range
  //we can get this range automatically, by evaluating the central (median) integral of the source distribution
  //with this set up, we fit the 68% most central yield of the source distribution
  double lowerlimit;
  double upperlimit;
  GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
  unsigned lowerbin = hSource->FindBin(lowerlimit);
  unsigned upperbin = hSource->FindBin(upperlimit);
  for(unsigned uBin=lowerbin; uBin<=upperbin; uBin++){
    hSource->SetBinError(uBin+1,0.01);
  }

  printf("Core size of %.3f fm\n",CoreSize);
  printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
  //fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
  TF1* fSource = new TF1("fSource",GaussSourceTF1,rMin,rMax,1);
  fSource->SetParameter(0,CoreSize);
  fSource->SetParLimits(0,CoreSize*0.5,CoreSize*2.0);
  hSource->Fit(fSource,"S, N, R, M");
  printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));

  //get rid of weird plotting
  for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
    hSource->SetBinError(uBin+1,0.01);
  }
  hSource->Write();
  fSource->Write();

  delete hSource;
  delete fSource;
  delete fOutput;
}


void EffectiveGaussianpAL(){

//value taken from pp scaling (mT_scaling forlder)
  const double CoreSize = 0.993257;

  //DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
  DLM_CleverMcLevyResoTM MagicSource;

  //DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
  MagicSource.InitStability(1,2-1e-6,2+1e-6);
  MagicSource.InitScale(38,0.15,2.0);
  MagicSource.InitRad(257*2,0,64);
  MagicSource.InitType(2);
  ///////////////////

  //for p-antiL, set up the amount of secondaries
  //first for the protons (64.22%)
  MagicSource.SetUpReso(0,0.6422);
  //than for the Lambdas, here its 64.38%
  MagicSource.SetUpReso(1,0.6438);

  //the cut off scale in k*, for which the angular distributions from EPOS
  //are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
  const double k_CutOff = 200;

  //to be used for the NTuple later on
  Float_t k_D;
  Float_t fP1;
  Float_t fP2;
  Float_t fM1;
  Float_t fM2;
  Float_t Tau1;
  Float_t Tau2;
  Float_t AngleRcP1;
  Float_t AngleRcP2;
  Float_t AngleP1P2;
  //random generator dimi style. The input is incompatible with the ROOT random generator,
  //do not mix and match, do not ask me how I know this. Ask Bernie.
  //11 is the seed, you can change that to you favorite number
  DLM_Random RanGen(11);
  //dummies to save random shit
  double RanVal1;
  double RanVal2;
  double RanVal3;

// 1. Reso p and primordials Lambdas
  //open the magic file from dimi with the angular distributions.
  TFile* F_EposDisto_pReso_Lam = new TFile("/Users/sartozza/cernbox/SourceStudies/EPOS_AngDistrib/EposDisto_pReso_Lam.root");
   //set up the ntuple, do not change anything unless told so by dimi
  TNtuple* T_EposDisto_pReso_Lam = (TNtuple*)F_EposDisto_pReso_Lam->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_pReso_Lam = T_EposDisto_pReso_Lam->GetEntries();
  T_EposDisto_pReso_Lam->SetBranchAddress("k_D",&k_D);
  T_EposDisto_pReso_Lam->SetBranchAddress("P1",&fP1);
  T_EposDisto_pReso_Lam->SetBranchAddress("P2",&fP2);
  T_EposDisto_pReso_Lam->SetBranchAddress("M1",&fM1);
  T_EposDisto_pReso_Lam->SetBranchAddress("M2",&fM2);
  T_EposDisto_pReso_Lam->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_pReso_Lam->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_pReso_Lam->SetBranchAddress("AngleP1P2",&AngleP1P2);
  //iterate over the ntuple
  for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_Lam; uEntry++){
      //get each entry
      T_EposDisto_pReso_Lam->GetEntry(uEntry);
      //disregard the entry of you are outside the desired k*
      if(k_D>k_CutOff) continue;
      //overwrite the value for the lifetime. This is computed from the
      //stat. hadronization model (Vale) or thermal fist (Max)
      //this is the value for the secondary protons
      Tau1 = 1.65;
      //for primoridials (the Lambdas) we put 0
      Tau2 = 0;
      //put in the average mass of the resonances (again from SHM or TF)
      //this is the value for protons
      fM1 = 1362;
      //generate a random path length for the propagation of the resonances
      //nothing to change!
      RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
      //adds a single entry into the PDF for the angular distribution to be used
      MagicSource.AddBGT_RP(RanVal1,cos(AngleRcP1));
  }
  delete F_EposDisto_pReso_Lam;

// 2. Primordials p and Reso Lambdas
  //open the magic file from dimi with the angular distributions.
  TFile* F_EposDisto_p_LamReso = new TFile("/Users/sartozza/cernbox/SourceStudies/EPOS_AngDistrib/EposDisto_p_LamReso.root");
   //set up the ntuple, do not change anything unless told so by dimi
  TNtuple* T_EposDisto_p_LamReso = (TNtuple*)F_EposDisto_p_LamReso->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_p_LamReso = T_EposDisto_p_LamReso->GetEntries();
  T_EposDisto_p_LamReso->SetBranchAddress("k_D",&k_D);
  T_EposDisto_p_LamReso->SetBranchAddress("P1",&fP1);
  T_EposDisto_p_LamReso->SetBranchAddress("P2",&fP2);
  T_EposDisto_p_LamReso->SetBranchAddress("M1",&fM1);
  T_EposDisto_p_LamReso->SetBranchAddress("M2",&fM2);
  T_EposDisto_p_LamReso->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_p_LamReso->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_p_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
  //iterate over the ntuple
  for(unsigned uEntry=0; uEntry<N_EposDisto_p_LamReso; uEntry++){
      //get each entry
      T_EposDisto_p_LamReso->GetEntry(uEntry);
      //disregard the entry of you are outside the desired k*
      if(k_D>k_CutOff) continue;
      //overwrite the value for the lifetime. This is computed from the
      //stat. hadronization model (Vale) or thermal fist (Max)
      //this is the value for the primordials protons
      Tau1 = 0;
      //for secondaries (the Lambdas) we put 0
      Tau2 = 4.69;
      //put in the average mass of the resonances (again from SHM or TF)
      //this is the value for protons
      fM2 = 1462;
      //generate a random path length for the propagation of the resonances
      //nothing to change!
      RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
      //adds a single entry into the PDF for the angular distribution to be used
      MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
  }
  delete F_EposDisto_p_LamReso;

// 3. Reso p and Reso Lambdas
  //open the magic file from dimi with the angular distributions.
  TFile* F_EposDisto_pReso_LamReso = new TFile("/Users/sartozza/cernbox/SourceStudies/EPOS_AngDistrib/EposDisto_pReso_LamReso.root");
   //set up the ntuple, do not change anything unless told so by dimi
  TNtuple* T_EposDisto_pReso_LamReso = (TNtuple*)F_EposDisto_pReso_LamReso->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_pReso_LamReso = T_EposDisto_pReso_LamReso->GetEntries();
  T_EposDisto_pReso_LamReso->SetBranchAddress("k_D",&k_D);
  T_EposDisto_pReso_LamReso->SetBranchAddress("P1",&fP1);
  T_EposDisto_pReso_LamReso->SetBranchAddress("P2",&fP2);
  T_EposDisto_pReso_LamReso->SetBranchAddress("M1",&fM1);
  T_EposDisto_pReso_LamReso->SetBranchAddress("M2",&fM2);
  T_EposDisto_pReso_LamReso->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_pReso_LamReso->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_pReso_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
  //iterate over the ntuple
  for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_LamReso; uEntry++){
      //get each entry
      T_EposDisto_pReso_LamReso->GetEntry(uEntry);
      //disregard the entry of you are outside the desired k*
      if(k_D>k_CutOff) continue;
      //overwrite the value for the lifetime. This is computed from the
      //stat. hadronization model (Vale) or thermal fist (Max)
      //this is the value for the secondary protons
      Tau1 = 1.65;
      //for primoridials (the Lambdas) we put 0
      Tau2 = 4.69;
      //put in the average mass of the resonances (again from SHM or TF)
      //this is the value for protons
      fM1 = 1362;
      fM2 = 1462;
      //generate a random path length for the propagation of the resonances
      //nothing to change!
      RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
      RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
      //adds a single entry into the PDF for the angular distribution to be used
      MagicSource.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
  }
  delete F_EposDisto_pReso_LamReso;
  //if you have resonances contributing to both particles, we need to repeat the above procedure
  //for the prim-reso (AddBGT_PR) and reso-reso (AddBGT_RR) cases

  const unsigned NumSourceBins = 128;
  const double rMin = 0;
  const double rMax = 16;
  TFile* fOutput = new TFile("fOutputpAL.root","recreate");
  TH1F* hSource = new TH1F("hSource","hSource",NumSourceBins,rMin,rMax);

  //fill the histo fro the source
  for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
    //get the x-axis (r value) of the current bin
    double xaxis = hSource->GetBinCenter(uBin+1);
    //an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
    double parameters[2];
    parameters[0] = CoreSize;
    parameters[1] = 2.0;
    double SourceValue = MagicSource.RootEval(&xaxis, parameters);
    hSource->SetBinContent(uBin+1,SourceValue);
    //infinite errors for now
    hSource->SetBinError(uBin+1,1000.);
  }

  //idea: fit the source distribution only in a range around its peak
  //to do this: silly idea: put very large uncertainties in the bins outside of this range
  //we can get this range automatically, by evaluating the central (median) integral of the source distribution
  //with this set up, we fit the 68% most central yield of the source distribution
  double lowerlimit;
  double upperlimit;
  GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
  unsigned lowerbin = hSource->FindBin(lowerlimit);
  unsigned upperbin = hSource->FindBin(upperlimit);
  for(unsigned uBin=lowerbin; uBin<=upperbin; uBin++){
    hSource->SetBinError(uBin+1,0.01);
  }

  printf("Core size of %.3f fm\n",CoreSize);
  printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
  //fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
  TF1* fSource = new TF1("fSource",GaussSourceTF1,rMin,rMax,1);
  fSource->SetParameter(0,CoreSize);
  fSource->SetParLimits(0,CoreSize*0.5,CoreSize*2.0);
  hSource->Fit(fSource,"S, N, R, M");
  printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));

  //get rid of weird plotting
  for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
    hSource->SetBinError(uBin+1,0.01);
  }
  hSource->Write();
  fSource->Write();

  delete hSource;
  delete fSource;
  delete fOutput;
}


void EffectiveGaussianLAL(){

//value taken from pp scaling (mT_scaling forlder)
  const double CoreSize = 0.92683;

  //DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
  DLM_CleverMcLevyResoTM MagicSource;

  //DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
  MagicSource.InitStability(1,2-1e-6,2+1e-6);
  MagicSource.InitScale(38,0.15,2.0);
  MagicSource.InitRad(257*2,0,64);
  MagicSource.InitType(2);
  ///////////////////

  //for p-antiL, set up the amount of secondaries
  //first for the lambdas (64.38%)
  MagicSource.SetUpReso(0,0.6438);
  //than for the Lambdas, here its 64.38%
  MagicSource.SetUpReso(1,0.6438);

  //the cut off scale in k*, for which the angular distributions from EPOS
  //are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
  const double k_CutOff = 200;

  //to be used for the NTuple later on
  Float_t k_D;
  Float_t fP1;
  Float_t fP2;
  Float_t fM1;
  Float_t fM2;
  Float_t Tau1;
  Float_t Tau2;
  Float_t AngleRcP1;
  Float_t AngleRcP2;
  Float_t AngleP1P2;
  //random generator dimi style. The input is incompatible with the ROOT random generator,
  //do not mix and match, do not ask me how I know this. Ask Bernie.
  //11 is the seed, you can change that to you favorite number
  DLM_Random RanGen(11);
  //dummies to save random shit
  double RanVal1;
  double RanVal2;
  double RanVal3;

// 1. Primordials Lambdas and Reso Lambdas
  //open the magic file from dimi with the angular distributions.
  TFile* F_EposDisto_Lam_LamReso = new TFile("/Users/sartozza/cernbox/SourceStudies/EPOS_AngDistrib/EposDisto_Lam_LamReso.root");
   //set up the ntuple, do not change anything unless told so by dimi
  TNtuple* T_EposDisto_Lam_LamReso = (TNtuple*)F_EposDisto_Lam_LamReso->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_Lam_LamReso = T_EposDisto_Lam_LamReso->GetEntries();
  T_EposDisto_Lam_LamReso->SetBranchAddress("k_D",&k_D);
  T_EposDisto_Lam_LamReso->SetBranchAddress("P1",&fP1);
  T_EposDisto_Lam_LamReso->SetBranchAddress("P2",&fP2);
  T_EposDisto_Lam_LamReso->SetBranchAddress("M1",&fM1);
  T_EposDisto_Lam_LamReso->SetBranchAddress("M2",&fM2);
  T_EposDisto_Lam_LamReso->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_Lam_LamReso->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_Lam_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_Lam_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_Lam_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
  //iterate over the ntuple
  for(unsigned uEntry=0; uEntry<N_EposDisto_Lam_LamReso; uEntry++){
      //get each entry
      T_EposDisto_Lam_LamReso->GetEntry(uEntry);
      //disregard the entry of you are outside the desired k*
      if(k_D>k_CutOff) continue;
      //overwrite the value for the lifetime. This is computed from the
      //stat. hadronization model (Vale) or thermal fist (Max)
      //this is the value for the secondary protons
      Tau1 = 0.;
      //for primoridials (the Lambdas) we put 0
      Tau2 = 4.69;
      //put in the average mass of the resonances (again from SHM or TF)
      //this is the value for protons
      fM2 = 1462;
      //generate a random path length for the propagation of the resonances
      //nothing to change!
      RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
      //adds a single entry into the PDF for the angular distribution to be used
      MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
      MagicSource.AddBGT_RP(RanVal2,-cos(AngleRcP2));
  }
  delete F_EposDisto_Lam_LamReso;


// 3. Reso Lambda and Reso Lambdas
  //open the magic file from dimi with the angular distributions.
  TFile* F_EposDisto_LamReso_LamReso = new TFile("/Users/sartozza/cernbox/SourceStudies/EPOS_AngDistrib/EposDisto_LamReso_LamReso.root");
   //set up the ntuple, do not change anything unless told so by dimi
  TNtuple* T_EposDisto_LamReso_LamReso = (TNtuple*)F_EposDisto_LamReso_LamReso->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_LamReso_LamReso = T_EposDisto_LamReso_LamReso->GetEntries();
  T_EposDisto_LamReso_LamReso->SetBranchAddress("k_D",&k_D);
  T_EposDisto_LamReso_LamReso->SetBranchAddress("P1",&fP1);
  T_EposDisto_LamReso_LamReso->SetBranchAddress("P2",&fP2);
  T_EposDisto_LamReso_LamReso->SetBranchAddress("M1",&fM1);
  T_EposDisto_LamReso_LamReso->SetBranchAddress("M2",&fM2);
  T_EposDisto_LamReso_LamReso->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_LamReso_LamReso->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_LamReso_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_LamReso_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_LamReso_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
  //iterate over the ntuple
  for(unsigned uEntry=0; uEntry<N_EposDisto_LamReso_LamReso; uEntry++){
      //get each entry
      T_EposDisto_LamReso_LamReso->GetEntry(uEntry);
      //disregard the entry of you are outside the desired k*
      if(k_D>k_CutOff) continue;
      //overwrite the value for the lifetime. This is computed from the
      //stat. hadronization model (Vale) or thermal fist (Max)
      //this is the value for the secondary protons
      Tau1 = 4.69;
      //for primoridials (the Lambdas) we put 0
      Tau2 = 4.69;
      //put in the average mass of the resonances (again from SHM or TF)
      //this is the value for protons
      fM1 = 1462;
      fM2 = 1462;
      //generate a random path length for the propagation of the resonances
      //nothing to change!
      RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
      RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
      //adds a single entry into the PDF for the angular distribution to be used
      MagicSource.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
  }
  delete F_EposDisto_LamReso_LamReso;
  //if you have resonances contributing to both particles, we need to repeat the above procedure
  //for the prim-reso (AddBGT_PR) and reso-reso (AddBGT_RR) cases

  const unsigned NumSourceBins = 128;
  const double rMin = 0;
  const double rMax = 16;
  TFile* fOutput = new TFile("fOutputLAL.root","recreate");
  TH1F* hSource = new TH1F("hSource","hSource",NumSourceBins,rMin,rMax);

  //fill the histo fro the source
  for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
    //get the x-axis (r value) of the current bin
    double xaxis = hSource->GetBinCenter(uBin+1);
    //an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
    double parameters[2];
    parameters[0] = CoreSize;
    parameters[1] = 2.0;
    double SourceValue = MagicSource.RootEval(&xaxis, parameters);
    hSource->SetBinContent(uBin+1,SourceValue);
    //infinite errors for now
    hSource->SetBinError(uBin+1,1000.);
  }

  //idea: fit the source distribution only in a range around its peak
  //to do this: silly idea: put very large uncertainties in the bins outside of this range
  //we can get this range automatically, by evaluating the central (median) integral of the source distribution
  //with this set up, we fit the 68% most central yield of the source distribution
  double lowerlimit;
  double upperlimit;
  GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
  unsigned lowerbin = hSource->FindBin(lowerlimit);
  unsigned upperbin = hSource->FindBin(upperlimit);
  for(unsigned uBin=lowerbin; uBin<=upperbin; uBin++){
    hSource->SetBinError(uBin+1,0.01);
  }

  printf("Core size of %.3f fm\n",CoreSize);
  printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
  //fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
  TF1* fSource = new TF1("fSource",GaussSourceTF1,rMin,rMax,1);
  fSource->SetParameter(0,CoreSize);
  fSource->SetParLimits(0,CoreSize*0.5,CoreSize*2.0);
  hSource->Fit(fSource,"S, N, R, M");
  printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));

  //get rid of weird plotting
  for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
    hSource->SetBinError(uBin+1,0.01);
  }
  hSource->Write();
  fSource->Write();

  delete hSource;
  delete fSource;
  delete fOutput;
}


int CalculateSource(int argc, char *argv[]){
    printf("Suck my cock\n");
    DLM_Timer TIMER;

    // EffectiveGaussianXiCock();

    // EffectiveGaussianpAL();

    EffectiveGaussianLAL();


    long long ExeTime = TIMER.Stop()/1000.;
    char* strtime = new char [128];
    ShowTime(ExeTime,strtime,0,true,6);
    printf("The script terminated after: %s\n",strtime);
    delete [] strtime;
    return 0;
}
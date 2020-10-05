#include "CATS.h"
#include "CommonAnaFunctions.h"
#include "DLM_Potentials.h"
#include "DLM_Random.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_WfModel.h"
#include "DLM_CkModels.h"
#include "TGraph.h"
#include "TGraphErrors.h"
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
#include "TFractionFitter.h"
#include "TGenPhaseSpace.h"
#include "TString.h"
#include "DLM_Fitters.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLatex.h"
#include "DLM_Integration.h"
#include "CATSconstants.h"
#include "DLM_SubPads.h"
#include <math.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TGaxis.h"

void SetStyle(bool graypalette, bool title)
{
      const int NCont = 255;
      gStyle->Reset("Plain");
      gStyle->SetNumberContours(NCont);
      gStyle->SetOptTitle(title);
      gStyle->SetTitleBorderSize(0);
      gStyle->SetOptStat(0);
      if(graypalette) gStyle->SetPalette(8,0);
      else gStyle->SetPalette(1);
      gStyle->SetCanvasColor(10);
      gStyle->SetCanvasBorderMode(0);
      gStyle->SetFrameLineWidth(1);
      gStyle->SetFrameFillColor(kWhite);
      gStyle->SetPadColor(10);
      gStyle->SetPadTickX(1);
      // gStyle->SetPadTickY(1);
      gStyle->SetPadBottomMargin(0.15);
      gStyle->SetPadLeftMargin(0.15);
      gStyle->SetHistLineWidth(1);
      gStyle->SetHistLineColor(kRed);
      gStyle->SetFuncWidth(2);
      gStyle->SetFuncColor(kGreen);
      gStyle->SetLineWidth(2);
      gStyle->SetLabelSize(0.045,"x");
      gStyle->SetLabelSize(0.045,"y");
      gStyle->SetLabelOffset(0.01,"y");
      gStyle->SetLabelOffset(0.01,"x");
      gStyle->SetLabelColor(kBlack,"xyz");
      gStyle->SetTitleSize(0.05,"xyz");
      gStyle->SetTitleOffset(1.4,"y");
      gStyle->SetTitleOffset(1.2,"x");
      gStyle->SetTitleFillColor(kWhite);
      gStyle->SetTextSizePixels(26);
      gStyle->SetTextFont(42);
      gStyle->SetLegendBorderSize(0);
      gStyle->SetLegendFillColor(kWhite);
      gStyle->SetLegendFont(42);
      gStyle->SetLegendBorderSize(0);

}

void PlotGraph(TGraph& graph, TString name, bool wf, int color, int linew, int lines){
  graph.SetName(name);
  // if(!wf) graph.SetTitle(name+";k* (MeV/c); C(k*)");
  // else graph.SetTitle(name +";r (fm); wf(r)");
  graph.SetLineColor(color);
  graph.SetLineWidth(linew);
  graph.SetLineStyle(lines);
  graph.SetMarkerColor(color);
}

void PlotHisto(TH1F& histo, TString name, bool wf, int color, int linew, int lines){
  histo.SetName(name);
  // if(!wf) histo.SetTitle(name+";k* (MeV/c); C(k*)");
  // else histo.SetTitle(name +";r (fm); wf(r)");
  histo.SetLineColor(color);
  histo.SetLineWidth(linew);
  histo.SetLineStyle(lines);
}


void PlotTF1(TF1& fun, TString name, bool wf, int color, int linew, int lines){
  fun.SetName(name);
  // if(!wf) histo.SetTitle(name+";k* (MeV/c); C(k*)");
  // else histo.SetTitle(name +";r (fm); wf(r)");
  fun.SetLineColor(color);
  fun.SetLineWidth(linew);
  fun.SetLineStyle(lines);
}


void Theo_pAntip(){

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);

    TFile* OutputFile; TFile* OutputFileWF;
    OutputFile = new TFile(TString::Format("pAp_Theo_Ck_2007_RadpAL.root"),"recreate");
    OutputFileWF = new TFile(TString::Format("pAp_Theo_WF2007_RadpAL.root"),"recreate");
    // SET THE GAUSSIAN SOURCE
    double radiuspAL = 0.;
    // radiuspAL = 1.146;//run2 HM based on avg mT
    radiuspAL = 1.70613;//run2 HM based on avg mT from pAL scaling
    printf("Radius for pAL = %.2f\n", radiuspAL);

    // SET THE BINNING
    double kmin_pAL = 0.;
    double kmax_pAL = 400.;
    unsigned NumMomBins_pAL = int((kmax_pAL)/2.);
    printf("........Setting the theoretical correlation functions...........\n");
    TString DataSample = "pp13TeV_HM_BBar";
    DLM_CommonAnaFunctions AnalysisObject;

//ch 0: 1S0 +1P1 -> 1/4
//ch 1: 3S1 + 3P0 -> 1/12
//ch 2: 3S1 + 3P1 -> 1/4
//ch 3: 3S1 + 3P2 -> 5/12

//nAn coupled-channel
//ch 4: 1S0 +1P1 -> 1/4
//ch 5: 3S1 + 3P0 -> 1/12
//ch 6: 3S1 + 3P1 -> 1/4
//ch 7: 3S1 + 3P2 -> 5/12

 // Setting Up the pAp CATS
  //With all channels to compare wfs
    DLM_Ck* Ck_pAp_All;//only pAp
    CATS AB_pAp_All;
    DLM_Ck* Ck_pAp_All1;// pAp+nAn
    CATS AB_pAp_All1;
  //With only S waves: Type 0-> only pAp, Type 1 pAp+nAn
    DLM_Ck* Ck_pAp0_S;//only pAp
    CATS AB_pAp0_S;
    DLM_Ck* Ck_pAp1_S;//+n-antin
    CATS AB_pAp1_S;
  // only 1S0
    DLM_Ck* Ck_pAp0_1S0;//only pAp
    CATS AB_pAp0_1S0;
  //only 3S1
    DLM_Ck* Ck_pAp0_3S1;//only pAp
    CATS AB_pAp0_3S1;

  //With 3P0 contributions  Type 0-> only pAp, Type 1 pAp+nAn
    DLM_Ck* Ck_pAp0_P;//only pAp
    CATS AB_pAp0_P;
    DLM_Ck* Ck_pAp1_P;//+n-antin
    CATS AB_pAp1_P;
  //only 3P0
    DLM_Ck* Ck_pAp0_3P0;//only pAp
    CATS AB_pAp0_3P0;

  //only 3S1+3P0
    DLM_Ck* Ck_pAp0_3S13P0;//only pAp
    CATS AB_pAp0_3S13P0;
  //only 3S1+3P1
    DLM_Ck* Ck_pAp0_3S13P1;//only pAp
    CATS AB_pAp0_3S13P1;
  //only 3S1+3P2
    DLM_Ck* Ck_pAp0_3S13P2;//only pAp
    CATS AB_pAp0_3S13P2;

    AB_pAp_All.SetMomBins(NumMomBins_pAL,kmin_pAL,kmax_pAL);
    AB_pAp_All1.SetMomBins(NumMomBins_pAL,kmin_pAL,kmax_pAL);

    AB_pAp0_S.SetMomBins(NumMomBins_pAL,kmin_pAL,kmax_pAL);
    AB_pAp1_S.SetMomBins(NumMomBins_pAL,kmin_pAL,kmax_pAL);

    AB_pAp0_1S0.SetMomBins(NumMomBins_pAL,kmin_pAL,kmax_pAL);
    AB_pAp0_3S1.SetMomBins(NumMomBins_pAL,kmin_pAL,kmax_pAL);
    AB_pAp0_P.SetMomBins(NumMomBins_pAL,kmin_pAL,kmax_pAL);
    AB_pAp0_3P0.SetMomBins(NumMomBins_pAL,kmin_pAL,kmax_pAL);
    AB_pAp0_3S13P0.SetMomBins(NumMomBins_pAL,kmin_pAL,kmax_pAL);
    AB_pAp0_3S13P1.SetMomBins(NumMomBins_pAL,kmin_pAL,kmax_pAL);
    AB_pAp0_3S13P2.SetMomBins(NumMomBins_pAL,kmin_pAL,kmax_pAL);

    printf("........Haidenbauer p-antip theoretical correlation from CATS...........\n");
//for wfs
    AnalysisObject.SetUpCats_pApHaide(AB_pAp_All,"HAIDE_1","Gauss", DataSample);
    AB_pAp_All.SetNotifications(CATS::nAll);
    AB_pAp_All.KillTheCat();
    Ck_pAp_All = new DLM_Ck(AB_pAp_All.GetNumSourcePars(),0,AB_pAp_All);//#source pars,#pot pars
// pAp+nAn
    AnalysisObject.SetUpCats_pApHaide(AB_pAp_All1,"HAIDE_2","Gauss", DataSample);
    AB_pAp_All1.SetNotifications(CATS::nAll);
    AB_pAp_All1.KillTheCat();
    Ck_pAp_All1 = new DLM_Ck(AB_pAp_All1.GetNumSourcePars(),0,AB_pAp_All1);//#source pars,#pot pars
    // Ck_pAp_All1->SetCutOff(400.,600.);//kc matching point at 500
// only S states pAp
    AnalysisObject.SetUpCats_pApHaide(AB_pAp0_S,"HAIDE_1","Gauss", DataSample);
    AB_pAp0_S.SetNotifications(CATS::nAll);
    AB_pAp0_S.RemoveExternalWaveFunction(0, 1);
    AB_pAp0_S.RemoveExternalWaveFunction(1, 1);
    AB_pAp0_S.RemoveExternalWaveFunction(2, 1);
    AB_pAp0_S.RemoveExternalWaveFunction(3, 1);

    AB_pAp0_S.SetChannelWeight(0,1./4.);
    AB_pAp0_S.SetChannelWeight(1,3./4.);
    AB_pAp0_S.SetChannelWeight(2,0.);
    AB_pAp0_S.SetChannelWeight(3,0.);

    AB_pAp0_S.KillTheCat();
    Ck_pAp0_S = new DLM_Ck(AB_pAp0_S.GetNumSourcePars(),0,AB_pAp0_S);//#source pars,#pot pars
// pAp+nAn only S states - now we have 6 channels!!
    AnalysisObject.SetUpCats_pApHaide(AB_pAp1_S,"HAIDE_2","Gauss", DataSample);
    AB_pAp1_S.SetNotifications(CATS::nAll);
    AB_pAp1_S.RemoveExternalWaveFunction(0, 1);
    AB_pAp1_S.RemoveExternalWaveFunction(1, 1);
    AB_pAp1_S.RemoveExternalWaveFunction(2, 1);
    AB_pAp1_S.RemoveExternalWaveFunction(3, 1);
    AB_pAp1_S.RemoveExternalWaveFunction(4, 1);
    AB_pAp1_S.RemoveExternalWaveFunction(5, 1);
    AB_pAp1_S.RemoveExternalWaveFunction(6, 1);
    AB_pAp1_S.RemoveExternalWaveFunction(7, 1);
    AB_pAp1_S.SetChannelWeight(0,1./4.);
    AB_pAp1_S.SetChannelWeight(1,3./4.);
    AB_pAp1_S.SetChannelWeight(2,0.);
    AB_pAp1_S.SetChannelWeight(3,0.);

    AB_pAp1_S.SetChannelWeight(4,1./4.);
    AB_pAp1_S.SetChannelWeight(5,3./4.);
    AB_pAp1_S.SetChannelWeight(6,0.);
    AB_pAp1_S.SetChannelWeight(7,0.);
    AB_pAp1_S.KillTheCat();
    Ck_pAp1_S = new DLM_Ck(AB_pAp1_S.GetNumSourcePars(),0,AB_pAp1_S);//#source pars,#pot pars
    // Ck_pAp1_S->SetCutOff(400.,500.);//kc matching point at 500

// only 1S0
    AnalysisObject.SetUpCats_pApHaide(AB_pAp0_1S0,"HAIDE_1","Gauss", DataSample);
    AB_pAp0_1S0.SetNotifications(CATS::nAll);
    AB_pAp0_1S0.RemoveExternalWaveFunction(0, 1);//removing 1P1
    AB_pAp0_1S0.SetChannelWeight(0,1.);
    AB_pAp0_1S0.SetChannelWeight(1,0.);
    AB_pAp0_1S0.SetChannelWeight(2,0.);
    AB_pAp0_1S0.SetChannelWeight(3,0.);
    AB_pAp0_1S0.KillTheCat();
    Ck_pAp0_1S0 = new DLM_Ck(AB_pAp0_1S0.GetNumSourcePars(),0,AB_pAp0_1S0);//#source pars,#pot pars
// only 3S1
    AnalysisObject.SetUpCats_pApHaide(AB_pAp0_3S1,"HAIDE_1","Gauss", DataSample);
    AB_pAp0_3S1.SetNotifications(CATS::nAll);
    AB_pAp0_3S1.RemoveExternalWaveFunction(1, 1);//removing 3P0
    AB_pAp0_3S1.SetChannelWeight(0,0.);
    AB_pAp0_3S1.SetChannelWeight(1,1.);
    AB_pAp0_3S1.SetChannelWeight(2,0.);
    AB_pAp0_3S1.SetChannelWeight(3,0.);
    AB_pAp0_3S1.KillTheCat();
    Ck_pAp0_3S1 = new DLM_Ck(AB_pAp0_3S1.GetNumSourcePars(),0,AB_pAp0_3S1);//#source pars,#pot pars
    // AnalysisObject.SetUpCats_pApHaide(AB_pAp1_S,"HAIDE_2","Gauss", DataSample);
    // AB_pAp1_S.SetNotifications(CATS::nAll);
    // AB_pAp1_S.SetChannelWeight(0,1./4.);
    // AB_pAp1_S.SetChannelWeight(1,0.);
    // AB_pAp1_S.SetChannelWeight(2,3./4.);
    // AB_pAp1_S.KillTheCat();
    // Ck_pAp1_S = new DLM_Ck(AB_pAp1_S.GetNumSourcePars(),0,AB_pAp1_S);//#source pars,#pot pars

// //with all P waves
//     AnalysisObject.SetUpCats_pApHaide(AB_pAp0_P,"HAIDE_1","Gauss", DataSample);
//     AB_pAp0_P.SetNotifications(CATS::nAll);
//     AB_pAp0_P.KillTheCat();
//     Ck_pAp0_P = new DLM_Ck(AB_pAp0_P.GetNumSourcePars(),0,AB_pAp0_P);//#source pars,#pot pars
//   //only 3P0
//     AnalysisObject.SetUpCats_pApHaide(AB_pAp0_3P0,"HAIDE_1","Gauss", DataSample);
//     AB_pAp0_3P0.SetNotifications(CATS::nAll);
//     AB_pAp0_3P0.RemoveExternalWaveFunction(1, 0);
//     AB_pAp0_3P0.SetChannelWeight(0,0.);
//     AB_pAp0_3P0.SetChannelWeight(1,1.);
//     AB_pAp0_3P0.SetChannelWeight(2,0.);
//     AB_pAp0_3P0.KillTheCat();
//     Ck_pAp0_3P0 = new DLM_Ck(AB_pAp0_3P0.GetNumSourcePars(),0,AB_pAp0_3P0);//#source pars,#pot pars

  //3S1+3P0
    AnalysisObject.SetUpCats_pApHaide(AB_pAp0_3S13P0,"HAIDE_1","Gauss", DataSample);
    AB_pAp0_3S13P0.SetNotifications(CATS::nAll);
    AB_pAp0_3S13P0.SetChannelWeight(0,0.);
    AB_pAp0_3S13P0.SetChannelWeight(1,1.);//C_3S1+C_3P0
    AB_pAp0_3S13P0.SetChannelWeight(2,0);
    AB_pAp0_3S13P0.SetChannelWeight(3,0);
    AB_pAp0_3S13P0.KillTheCat();
    Ck_pAp0_3S13P0 = new DLM_Ck(AB_pAp0_3S13P0.GetNumSourcePars(),0,AB_pAp0_3S13P0);//#source pars,#pot pars
  //3S1+3P1
    AnalysisObject.SetUpCats_pApHaide(AB_pAp0_3S13P1,"HAIDE_1","Gauss", DataSample);
    AB_pAp0_3S13P1.SetNotifications(CATS::nAll);
    AB_pAp0_3S13P1.SetChannelWeight(0,0.);
    AB_pAp0_3S13P1.SetChannelWeight(1,0.);//C_3S1+C_3P0
    AB_pAp0_3S13P1.SetChannelWeight(2,1.);
    AB_pAp0_3S13P1.SetChannelWeight(3,0.);
    AB_pAp0_3S13P1.KillTheCat();
    Ck_pAp0_3S13P1 = new DLM_Ck(AB_pAp0_3S13P1.GetNumSourcePars(),0,AB_pAp0_3S13P1);//#source pars,#pot pars
  //3S1+3P2
    AnalysisObject.SetUpCats_pApHaide(AB_pAp0_3S13P2,"HAIDE_1","Gauss", DataSample);
    AB_pAp0_3S13P2.SetNotifications(CATS::nAll);
    AB_pAp0_3S13P2.SetChannelWeight(0,0.);
    AB_pAp0_3S13P2.SetChannelWeight(1,0.);//C_3S1+C_3P0
    AB_pAp0_3S13P2.SetChannelWeight(2,0);
    AB_pAp0_3S13P2.SetChannelWeight(3,1.);
    AB_pAp0_3S13P2.KillTheCat();
    Ck_pAp0_3S13P2 = new DLM_Ck(AB_pAp0_3S13P2.GetNumSourcePars(),0,AB_pAp0_3S13P2);//#source pars,#pot pars
    // AnalysisObject.SetUpCats_pApHaide(AB_pAp1,"HAIDE_2","Gauss", DataSample);
    // AB_pAp1.SetNotifications(CATS::nAll);
    // AB_pAp1.KillTheCat();
    // Ck_pAp1 = new DLM_Ck(AB_pAp1.GetNumSourcePars(),0,AB_pAp1);//#source pars,#pot pars

     Ck_pAp_All->SetSourcePar(0,radiuspAL);
     Ck_pAp_All1->SetSourcePar(0,radiuspAL);

     Ck_pAp0_S->SetSourcePar(0,radiuspAL);
     Ck_pAp1_S->SetSourcePar(0,radiuspAL);

     Ck_pAp0_1S0->SetSourcePar(0,radiuspAL);
     Ck_pAp0_3S1->SetSourcePar(0,radiuspAL);

    //  Ck_pAp0_P->SetSourcePar(0,radiuspAL);
    //  Ck_pAp0_3P0->SetSourcePar(0,radiuspAL);
     Ck_pAp0_3S13P0->SetSourcePar(0,radiuspAL);
     Ck_pAp0_3S13P1->SetSourcePar(0,radiuspAL);
     Ck_pAp0_3S13P2->SetSourcePar(0,radiuspAL);

     Ck_pAp_All->Update();
     Ck_pAp_All1->Update();

     Ck_pAp0_S->Update();
     Ck_pAp1_S->Update();

     Ck_pAp0_1S0->Update();
     Ck_pAp0_3S1->Update();

    //  Ck_pAp0_P->Update();
    //  Ck_pAp0_3P0->Update();

     Ck_pAp0_3S13P0->Update();
     Ck_pAp0_3S13P1->Update();
     Ck_pAp0_3S13P2->Update();

    TGraph Graph_Ck_pAp_All;
    PlotGraph(Graph_Ck_pAp_All,"Graph_Ck_pAp_All",false, kRed,2,1);
    TGraph Graph_Ck_pAp_All1;
    PlotGraph(Graph_Ck_pAp_All1,"Graph_Ck_pAp_All1",false, kCyan+2,2,1);

    TGraph Graph_Ck_pAp_Type0_S;
    PlotGraph(Graph_Ck_pAp_Type0_S,"Graph_Ck_pAp_Type0_S",false, kRed+1,2,1);

    TGraph Graph_Ck_pAp_Type1_S;
    PlotGraph(Graph_Ck_pAp_Type1_S,"Graph_Ck_pAp_Type1_S",false, kCyan+2,2,1);

    TGraph Graph_Ck_pAp_Type0_1S0;
    PlotGraph(Graph_Ck_pAp_Type0_1S0,"Graph_Ck_pAp_Type0_1S0",false, kMagenta,2,1);
    TGraph Graph_Ck_pAp_Type0_3S1;
    PlotGraph(Graph_Ck_pAp_Type0_3S1,"Graph_Ck_pAp_Type0_3S1",false, kAzure,2,1);

    // TGraph Graph_Ck_pAp_Type0_P;
    // PlotGraph(Graph_Ck_pAp_Type0_P,"Graph_Ck_pAp_Type0_P",false, kBlue+1,2,1);
    // TGraph Graph_Ck_pAp_Type0_3P0;
    // PlotGraph(Graph_Ck_pAp_Type0_3P0,"Graph_Ck_pAp_Type0_3P0",false, kGreen+2,2,1);

    TGraph Graph_Ck_pAp_Type0_3S13P0;
    PlotGraph(Graph_Ck_pAp_Type0_3S13P0,"Graph_Ck_pAp_Type0_3S13P0",false, kOrange+7,2,1);
    TGraph Graph_Ck_pAp_Type0_3S13P1;
    PlotGraph(Graph_Ck_pAp_Type0_3S13P1,"Graph_Ck_pAp_Type0_3S13P1",false, kViolet+1,2,1);
    TGraph Graph_Ck_pAp_Type0_3S13P2;
    PlotGraph(Graph_Ck_pAp_Type0_3S13P2,"Graph_Ck_pAp_Type0_3S13P2",false, kCyan+2,2,1);

    TGraph Graph_Ratio_SP;
    PlotGraph(Graph_Ratio_SP,"Graph_Ratio_SP",false, kBlack+1,2,1);
    TGraph Graph_Ratio_SAll;
    PlotGraph(Graph_Ratio_SAll,"Graph_Ratio_SAll",false, kGray+1,2,1);



    for(unsigned k=0; k<NumMomBins_pAL; k++){
     Graph_Ck_pAp_All.SetPoint(k,Ck_pAp_All->GetBinCenter(0,k),Ck_pAp_All->GetBinContent(k));
     Graph_Ck_pAp_All1.SetPoint(k,Ck_pAp_All1->GetBinCenter(0,k),Ck_pAp_All1->GetBinContent(k));

     Graph_Ck_pAp_Type0_S.SetPoint(k,Ck_pAp0_S->GetBinCenter(0,k),Ck_pAp0_S->GetBinContent(k));
     Graph_Ck_pAp_Type1_S.SetPoint(k,Ck_pAp1_S->GetBinCenter(0,k),Ck_pAp1_S->GetBinContent(k));

     Graph_Ck_pAp_Type0_1S0.SetPoint(k,Ck_pAp0_1S0->GetBinCenter(0,k),Ck_pAp0_1S0->GetBinContent(k));
     Graph_Ck_pAp_Type0_3S1.SetPoint(k,Ck_pAp0_3S1->GetBinCenter(0,k),Ck_pAp0_3S1->GetBinContent(k));

    //  Graph_Ck_pAp_Type0_P.SetPoint(k,Ck_pAp0_P->GetBinCenter(0,k),Ck_pAp0_P->GetBinContent(k));
    //  Graph_Ck_pAp_Type0_3P0.SetPoint(k,Ck_pAp0_3P0->GetBinCenter(0,k),Ck_pAp0_3P0->GetBinContent(k));

     Graph_Ck_pAp_Type0_3S13P0.SetPoint(k,Ck_pAp0_3S13P0->GetBinCenter(0,k),Ck_pAp0_3S13P0->GetBinContent(k));
     Graph_Ck_pAp_Type0_3S13P1.SetPoint(k,Ck_pAp0_3S13P1->GetBinCenter(0,k),Ck_pAp0_3S13P1->GetBinContent(k));
     Graph_Ck_pAp_Type0_3S13P2.SetPoint(k,Ck_pAp0_3S13P2->GetBinCenter(0,k),Ck_pAp0_3S13P2->GetBinContent(k));

     Graph_Ratio_SAll.SetPoint(k,Ck_pAp0_S->GetBinCenter(0,k),Ck_pAp0_S->GetBinContent(k)/Ck_pAp_All->GetBinContent(k));
     }
    printf("Momentum = %.2f ---- C = %.4f \n", Ck_pAp_All1->GetBinCenter(0,199),Ck_pAp_All1->GetBinContent(199));


 //1. Checking the WFS
  double rmin = 0.05;
  double rmax = 10;
  unsigned count = 500;
  double step = (rmax-rmin)/double(count);
  double Radius;

  double ReWF2_1S0_pp [6]; double ImWF2_1S0_pp [6];
  double ReWF2_3S1_pp [6]; double ImWF2_3S1_pp [6];
  double ReWF2_3P0_pp [6]; double ImWF2_3P0_pp [6];

  double WF2_1S0_pp [6];
  double WF2_3S1_pp [6];
  double WF2_3P0_pp [6];

  double ReWF2_1P1_pp [6]; double ImWF2_1P1_pp [6];
  double ReWF2_3P1_pp [6]; double ImWF2_3P1_pp [6];
  double ReWF2_3P2_pp [6]; double ImWF2_3P2_pp [6];

  double WF2_1P1_pp [6];
  double WF2_3P1_pp [6];
  double WF2_3P2_pp [6];

  TGraph Graph_ReWF2_1S0_pp [6];
  TGraph Graph_ReWF2_3S1_pp [6];
  TGraph Graph_ReWF2_3P0_pp [6];
  TGraph Graph_ReWF2_1P1_pp [6];
  TGraph Graph_ReWF2_3P1_pp [6];
  TGraph Graph_ReWF2_3P2_pp [6];

  TGraph Graph_ImWF2_1S0_pp [6];
  TGraph Graph_ImWF2_3S1_pp [6];
  TGraph Graph_ImWF2_3P0_pp [6];
  TGraph Graph_ImWF2_1P1_pp [6];
  TGraph Graph_ImWF2_3P1_pp [6];
  TGraph Graph_ImWF2_3P2_pp [6];

  TGraph Graph_WF2_1S0_pp [6];
  TGraph Graph_WF2_3S1_pp [6];
  TGraph Graph_WF2_3P0_pp [6];
  TGraph Graph_WF2_1P1_pp [6];
  TGraph Graph_WF2_3P1_pp [6];
  TGraph Graph_WF2_3P2_pp [6];

  // std::vector<int> momvec={2,10,20,30,38,39,40,41,42,43,44,45,50};
  std::vector<int> momvec={40,41,42,43,44,45};

  OutputFileWF->cd();
  for (unsigned k = 0; k<momvec.size(); k++){
    int mombinsel=momvec[k];
    printf("Momentum Selected = %.2f\n",Ck_pAp0_S->GetBinCenter(0,mombinsel));
    PlotGraph(Graph_ReWF2_1S0_pp[k],TString::Format("Graph_ReWF2_1S0_pp_%i",mombinsel),true,kRed+k-1,4,1);
    PlotGraph(Graph_ReWF2_3S1_pp[k],TString::Format("Graph_ReWF2_3S1_pp_%i",mombinsel),true,kBlue+k-1,4,1);
    PlotGraph(Graph_ReWF2_3P0_pp[k],TString::Format("Graph_ReWF2_3P0_pp_%i",mombinsel),true,kGreen+k-1,4,1);
    PlotGraph(Graph_ReWF2_1P1_pp[k],TString::Format("Graph_ReWF2_1P1_pp_%i",mombinsel),true,kViolet+k-1,4,1);
    PlotGraph(Graph_ReWF2_3P1_pp[k],TString::Format("Graph_ReWF2_3P1_pp_%i",mombinsel),true,kAzure+k-1,4,1);
    PlotGraph(Graph_ReWF2_3P2_pp[k],TString::Format("Graph_ReWF2_3P2_pp_%i",mombinsel),true,kOrange+k-1,4,1);

    PlotGraph(Graph_ImWF2_1S0_pp[k],TString::Format("Graph_ImWF2_1S0_pp_%i",mombinsel),true,kRed+k-1,4,7);
    PlotGraph(Graph_ImWF2_3S1_pp[k],TString::Format("Graph_ImWF2_3S1_pp_%i",mombinsel),true,kBlue+k-1,4,7);
    PlotGraph(Graph_ImWF2_3P0_pp[k],TString::Format("Graph_ImWF2_3P0_pp_%i",mombinsel),true,kGreen+k-1,4,7);
    PlotGraph(Graph_ImWF2_1P1_pp[k],TString::Format("Graph_ImWF2_1P1_pp_%i",mombinsel),true,kViolet+k-1,4,1);
    PlotGraph(Graph_ImWF2_3P1_pp[k],TString::Format("Graph_ImWF2_3P1_pp_%i",mombinsel),true,kAzure+k-1,4,1);
    PlotGraph(Graph_ImWF2_3P2_pp[k],TString::Format("Graph_ImWF2_3P2_pp_%i",mombinsel),true,kOrange+k-1,4,1);

    PlotGraph(Graph_WF2_1S0_pp[k],TString::Format("Graph_WF2_1S0_pp_%i",mombinsel),true,kRed+k-1,4,2);
    PlotGraph(Graph_WF2_3S1_pp[k],TString::Format("Graph_WF2_3S1_pp_%i",mombinsel),true,kBlue+k-1,4,2);
    PlotGraph(Graph_WF2_3P0_pp[k],TString::Format("Graph_WF2_3P0_pp_%i",mombinsel),true,kGreen+k-1,4,2);
    PlotGraph(Graph_WF2_1P1_pp[k],TString::Format("Graph_WF2_1P1_pp_%i",mombinsel),true,kViolet+k-1,4,1);
    PlotGraph(Graph_WF2_3P1_pp[k],TString::Format("Graph_WF2_3P1_pp_%i",mombinsel),true,kAzure+k-1,4,1);
    PlotGraph(Graph_WF2_3P2_pp[k],TString::Format("Graph_WF2_3P2_pp_%i",mombinsel),true,kOrange+k-1,4,1);

  for (unsigned i = 0; i < count; i++) {

  Radius = rmin+step*double(i);
  ReWF2_1S0_pp[k] = std::real(AB_pAp_All.EvalRadialWaveFunction(mombinsel,0,0,Radius,true));
  ReWF2_3S1_pp[k] = std::real(AB_pAp_All.EvalRadialWaveFunction(mombinsel,1,0,Radius,true));
  ReWF2_3P0_pp[k] = std::real(AB_pAp_All.EvalRadialWaveFunction(mombinsel,1,1,Radius,true));
  ReWF2_1P1_pp[k] = std::real(AB_pAp_All.EvalRadialWaveFunction(mombinsel,0,1,Radius,true));
  ReWF2_3P1_pp[k] = std::real(AB_pAp_All.EvalRadialWaveFunction(mombinsel,2,1,Radius,true));
  ReWF2_3P2_pp[k] = std::real(AB_pAp_All.EvalRadialWaveFunction(mombinsel,3,1,Radius,true));

  ImWF2_1S0_pp[k] = std::imag(AB_pAp_All.EvalRadialWaveFunction(mombinsel,0,0,Radius,true));
  ImWF2_3S1_pp[k] = std::imag(AB_pAp_All.EvalRadialWaveFunction(mombinsel,1,0,Radius,true));
  ImWF2_3P0_pp[k] = std::imag(AB_pAp_All.EvalRadialWaveFunction(mombinsel,1,1,Radius,true));
  ImWF2_1P1_pp[k] = std::imag(AB_pAp_All.EvalRadialWaveFunction(mombinsel,0,1,Radius,true));
  ImWF2_3P1_pp[k] = std::imag(AB_pAp_All.EvalRadialWaveFunction(mombinsel,2,1,Radius,true));
  ImWF2_3P2_pp[k] = std::imag(AB_pAp_All.EvalRadialWaveFunction(mombinsel,3,1,Radius,true));

  WF2_1S0_pp[k] = pow(ReWF2_1S0_pp[k],2)+pow(ImWF2_1S0_pp[k],2);
  WF2_3S1_pp[k] = pow(ReWF2_3S1_pp[k],2)+pow(ImWF2_3S1_pp[k],2);
  WF2_3P0_pp[k] = pow(ReWF2_3P0_pp[k],2)+pow(ImWF2_3P0_pp[k],2);
  WF2_1P1_pp[k] = pow(ReWF2_1P1_pp[k],2)+pow(ImWF2_1P1_pp[k],2);
  WF2_3P1_pp[k] = pow(ReWF2_3P1_pp[k],2)+pow(ImWF2_3P1_pp[k],2);
  WF2_3P2_pp[k] = pow(ReWF2_3P2_pp[k],2)+pow(ImWF2_3P2_pp[k],2);

  Graph_ReWF2_1S0_pp[k].SetPoint(i,Radius,ReWF2_1S0_pp[k]);
  Graph_ReWF2_3S1_pp[k].SetPoint(i,Radius,ReWF2_3S1_pp[k]);
  Graph_ReWF2_3P0_pp[k].SetPoint(i,Radius,ReWF2_3P0_pp[k]);
  Graph_ReWF2_1P1_pp[k].SetPoint(i,Radius,ReWF2_1P1_pp[k]);
  Graph_ReWF2_3P1_pp[k].SetPoint(i,Radius,ReWF2_3P1_pp[k]);
  Graph_ReWF2_3P2_pp[k].SetPoint(i,Radius,ReWF2_3P2_pp[k]);

  Graph_ImWF2_1S0_pp[k].SetPoint(i,Radius,ImWF2_1S0_pp[k]);
  Graph_ImWF2_3S1_pp[k].SetPoint(i,Radius,ImWF2_3S1_pp[k]);
  Graph_ImWF2_3P0_pp[k].SetPoint(i,Radius,ImWF2_3P0_pp[k]);
  Graph_ImWF2_1P1_pp[k].SetPoint(i,Radius,ImWF2_1P1_pp[k]);
  Graph_ImWF2_3P1_pp[k].SetPoint(i,Radius,ImWF2_3P1_pp[k]);
  Graph_ImWF2_3P2_pp[k].SetPoint(i,Radius,ImWF2_3P2_pp[k]);

  Graph_WF2_1S0_pp[k].SetPoint(i,Radius,WF2_1S0_pp[k]);
  Graph_WF2_3S1_pp[k].SetPoint(i,Radius,WF2_3S1_pp[k]);
  Graph_WF2_3P0_pp[k].SetPoint(i,Radius,WF2_3P0_pp[k]);
  Graph_WF2_1P1_pp[k].SetPoint(i,Radius,WF2_1P1_pp[k]);
  Graph_WF2_3P1_pp[k].SetPoint(i,Radius,WF2_3P1_pp[k]);
  Graph_WF2_3P2_pp[k].SetPoint(i,Radius,WF2_3P2_pp[k]);

  }
    Graph_ReWF2_1S0_pp[k].Write();
    Graph_ReWF2_3S1_pp[k].Write();
    Graph_ReWF2_3P0_pp[k].Write();
    Graph_ReWF2_1P1_pp[k].Write();
    Graph_ReWF2_3P1_pp[k].Write();
    Graph_ReWF2_3P2_pp[k].Write();

    Graph_ImWF2_1S0_pp[k].Write();
    Graph_ImWF2_3S1_pp[k].Write();
    Graph_ImWF2_3P0_pp[k].Write();
    Graph_ImWF2_1P1_pp[k].Write();
    Graph_ImWF2_3P1_pp[k].Write();
    Graph_ImWF2_3P2_pp[k].Write();

    Graph_WF2_1S0_pp[k].Write();
    Graph_WF2_3S1_pp[k].Write();
    Graph_WF2_3P0_pp[k].Write();
    Graph_WF2_1P1_pp[k].Write();
    Graph_WF2_3P1_pp[k].Write();
    Graph_WF2_3P2_pp[k].Write();

  }

  OutputFileWF->Close();

  OutputFile->cd();

  Graph_Ck_pAp_All.Write();
  Graph_Ck_pAp_All1.Write();

  Graph_Ck_pAp_Type0_S.Write();
  Graph_Ck_pAp_Type1_S.Write();

  Graph_Ck_pAp_Type0_1S0.Write();
  Graph_Ck_pAp_Type0_3S1.Write();
  // Graph_Ck_pAp_Type0_P.Write();
  // Graph_Ck_pAp_Type0_3P0.Write();
  Graph_Ck_pAp_Type0_3S13P0.Write();
  Graph_Ck_pAp_Type0_3S13P1.Write();
  Graph_Ck_pAp_Type0_3S13P2.Write();

  Graph_Ratio_SAll.Write();

  OutputFile->Close();
  delete OutputFileWF;
  delete OutputFile;


}


void ReviewFemto_Calculations(){

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);

    TFile* OutputFileCF; TFile* OutputFileWF; TFile* OutputFileSource;
    OutputFileCF = new TFile(TString::Format("fOutput_CF.root"),"recreate");
    OutputFileWF = new TFile(TString::Format("fOutput_WF.root"),"recreate");
    OutputFileSource = new TFile(TString::Format("fOutput_Source.root"),"recreate");

    // SET THE GAUSSIAN SOURCE in fm
    //For the moment I do not have the effective gaussian fitted to core + res for pΞ and pΩ
    //Using the value quoted in Nature of r_core, since resonances are only coming from protons
    double radius_pXi = 1.02;
    double radius_pOm = 0.95;
    //For K+ p we will use for now the results published in PRL where we use the radius from pp 13 TeV
    double radius_pKp = 1.18;

    // SET THE BINNING
    double kmin = 0.;
    double kmax = 400.;
    unsigned NumMomBins = int((kmax)/2.);//2 MeV binning

    printf("........Setting the theoretical correlation functions...........\n");
    TString DataSample = "FemtoReview";
    DLM_CommonAnaFunctions AnalysisObject;

    printf("........p-Ξ theoretical correlation from CATS...........\n");
    DLM_Ck* Ck_pXi;
    CATS AB_pXi;
    AB_pXi.SetMomBins(NumMomBins,kmin,kmax);
    //Choose between POT: pXim_HALQCDPaper2020, pXim_HALQCD1
    AnalysisObject.SetUpCats_pXim_Review(AB_pXi,"pXim_HALQCDPaper2020","Gauss");
    AB_pXi.SetNotifications(CATS::nAll);
    AB_pXi.KillTheCat();
    Ck_pXi = new DLM_Ck(AB_pXi.GetNumSourcePars(),0,AB_pXi);//#source pars,#pot pars
    Ck_pXi->SetSourcePar(0,radius_pXi);
    Ck_pXi->Update();

    printf("........p-Ω theoretical correlation from CATS...........\n");
    DLM_Ck* Ck_pOmega;
    CATS AB_pOmega;
    AB_pOmega.SetMomBins(NumMomBins,kmin,kmax);
    //Choose between POT: pXim_HALQCDPaper2020, pXim_HALQCD1
    AnalysisObject.SetUpCats_pOmega_Review(AB_pOmega,"pOmega_Lattice","Gauss");
    AB_pOmega.SetNotifications(CATS::nAll);
    AB_pOmega.KillTheCat();
    Ck_pOmega = new DLM_Ck(AB_pOmega.GetNumSourcePars(),0,AB_pOmega);//#source pars,#pot pars
    Ck_pOmega->SetSourcePar(0,radius_pOm);
    Ck_pOmega->Update();

    printf("........p-K+ theoretical correlation from CATS with I=1 only...........\n");
    DLM_Ck* Ck_pKpI1;
    CATS AB_pKpI1;
    AB_pKpI1.SetMomBins(NumMomBins,kmin,kmax);
    //Choose between POT: pXim_HALQCDPaper2020, pXim_HALQCD1
    AnalysisObject.SetUpCats_pKp_Review(AB_pKpI1,"Isospin1","Gauss");
    AB_pKpI1.SetNotifications(CATS::nAll);
    AB_pKpI1.KillTheCat();
    Ck_pKpI1 = new DLM_Ck(AB_pKpI1.GetNumSourcePars(),0,AB_pKpI1);//#source pars,#pot pars
    Ck_pKpI1->SetSourcePar(0,radius_pKp);
    Ck_pKpI1->Update();


    printf("........p-K+ theoretical correlation from CATS with I=1 and I=0 (RamonaLike)...........\n");
    DLM_Ck* Ck_pKpRamona;
    CATS AB_pKpRamona;
    AB_pKpRamona.SetMomBins(NumMomBins,kmin,kmax);
    //Choose between POT: pXim_HALQCDPaper2020, pXim_HALQCD1
    AnalysisObject.SetUpCats_pKp_Review(AB_pKpRamona,"RamonaLike","Gauss");
    AB_pKpRamona.SetNotifications(CATS::nAll);
    AB_pKpRamona.KillTheCat();
    Ck_pKpRamona = new DLM_Ck(AB_pKpRamona.GetNumSourcePars(),0,AB_pKpRamona);//#source pars,#pot pars
    Ck_pKpRamona->SetSourcePar(0,radius_pKp);
    Ck_pKpRamona->Update();

    TGraph Graph_Ck_pXi;
    PlotGraph(Graph_Ck_pXi,"Graph_Ck_pXi",false, kOrange+7,2,1);
    TGraph Graph_Ck_pOmega;
    PlotGraph(Graph_Ck_pOmega,"Graph_Ck_pOmega",false, kPink+5,2,1);
    TGraph Graph_Ck_pKpI1;
    PlotGraph(Graph_Ck_pKpI1,"Graph_Ck_pKpI1",false, kAzure+7,2,1);
    TGraph Graph_Ck_pKpRamona;
    PlotGraph(Graph_Ck_pKpRamona,"Graph_Ck_pKpRamona",false, kAzure+7,2,1);

    for(unsigned k=0; k<NumMomBins; k++){
     Graph_Ck_pXi.SetPoint(k,Ck_pXi->GetBinCenter(0,k),Ck_pXi->GetBinContent(k));
     Graph_Ck_pOmega.SetPoint(k,Ck_pOmega->GetBinCenter(0,k),Ck_pOmega->GetBinContent(k));
     Graph_Ck_pKpI1.SetPoint(k,Ck_pKpI1->GetBinCenter(0,k),Ck_pKpI1->GetBinContent(k));
     Graph_Ck_pKpRamona.SetPoint(k,Ck_pKpRamona->GetBinCenter(0,k),Ck_pKpRamona->GetBinContent(k));
    }
    // printf("Momentum = %.2f ---- C = %.4f \n", Ck_pXi->GetBinCenter(0,199),Ck_pXi->GetBinContent(199));


 //1. Checking the WFS
  double rmin = 0.05;
  double rmax = 10.05;
  unsigned count = 200;
  int nbins_rad = 200;
  double step = (rmax-rmin)/double(count);
  double Radius;

// we only wfs only for two values of k*
//k* = 10 and k*=100

//For pXi we have 4 channels (I=0,1; S=0,1)
//in the Nature paper we plot the potential for I=0,S=0 since it´s the most attractive
//since I predict the future I will plot them all
  double WF2_pXi_I0S0 [2];
  double WF2_pXi_I0S1 [2];
  double WF2_pXi_I1S0 [2];
  double WF2_pXi_I1S1 [2];

  double WF2_pOmega_5S2 [2];

  double WF2_pKp_I1 [2];
  double WF2_pKp_RamonaI0 [2];
  double WF2_pKp_RamonaI1 [2];


  double source_pXi;
  double source_pOmega;
  double source_pKp;

  TGraph Graph_WF2_pXi_I0S0 [2];  TGraph Graph_WF2_pXi_I0S0_renorm [2];
  TGraph Graph_WF2_pXi_I0S1 [2];  TGraph Graph_WF2_pXi_I0S1_renorm [2];
  TGraph Graph_WF2_pXi_I1S0 [2];  TGraph Graph_WF2_pXi_I1S0_renorm [2];
  TGraph Graph_WF2_pXi_I1S1 [2];  TGraph Graph_WF2_pXi_I1S1_renorm [2];

  TGraph Graph_WF2_pOmega_5S2 [2];  TGraph Graph_WF2_pOmega_5S2_renorm [2];

  TGraph Graph_WF2_pKp_I1 [2];  TGraph Graph_WF2_pKp_I1_renorm [2];
  TGraph Graph_WF2_pKp_RamonaI0 [2];
  TGraph Graph_WF2_pKp_RamonaI1 [2];

  // TH1F* h_WF2_pOmega_5S25 = new TH1F("h_WF2_pOmega_5S25","h_WF2_pOmega_5S25",nbins_rad,rmin,rmax);
  // TH1F* h_WF2_pOmega_5S250 = new TH1F("h_WF2_pOmega_5S250","h_WF2_pOmega_5S250",nbins_rad,rmin,rmax);

  // TH1F* h_WF2_pKp_I15 = new TH1F("h_WF2_pKp_I15","h_WF2_pKp_I15",nbins_rad,rmin,rmax);
  // TH1F* h_WF2_pKp_I150 = new TH1F("h_WF2_pKp_I150","h_WF2_pKp_I150",nbins_rad,rmin,rmax);

  TGraph Graph_Source_pXi;
  TGraph Graph_Source_pOmega;
  TGraph Graph_Source_pKp;

  TGraph Graph_Source_pXi_renorm;
  TGraph Graph_Source_pOmega_renorm;
  TGraph Graph_Source_pKp_renorm;

  TH1F* h_Source_pXi_renorm = new TH1F("h_Source_pXi_renorm","h_Source_pXi_renorm",nbins_rad,rmin,rmax);
  TH1F* h_Source_pOmega_renorm = new TH1F("h_Source_pOmega_renorm","h_Source_pOmega_renorm",nbins_rad,rmin,rmax);
  TH1F* h_Source_pKp_renorm = new TH1F("h_Source_pKp_renorm","h_Source_pKp_renorm",nbins_rad,rmin,rmax);

  PlotGraph(Graph_Source_pXi,TString::Format("Graph_Source_pXi"),true,kOrange+7,4,7);
  PlotGraph(Graph_Source_pOmega,TString::Format("Graph_Source_pOmega"),true,kPink+5,4,7);
  PlotGraph(Graph_Source_pKp,TString::Format("Graph_Source_pKp"),true,kAzure+7,4,7);
  PlotGraph(Graph_Source_pXi_renorm,TString::Format("Graph_Source_pXi_renorm"),true,kOrange+7,4,7);
  PlotGraph(Graph_Source_pOmega_renorm,TString::Format("Graph_Source_pOmega_renorm"),true,kPink+5,4,7);
  PlotGraph(Graph_Source_pKp_renorm,TString::Format("Graph_Source_pKp_renorm"),true,kAzure+7,4,7);

//Getting the source
  OutputFileSource->cd();
  for (unsigned i = 0; i < count; i++) {
  Radius = rmin+step*double(i);

  source_pXi = AB_pXi.EvaluateTheSource(0,Radius,0);
  Graph_Source_pXi.SetPoint(i,Radius,source_pXi);
  h_Source_pXi_renorm->SetBinContent(i+1,source_pXi);

  source_pOmega = AB_pOmega.EvaluateTheSource(0,Radius,0);
  Graph_Source_pOmega.SetPoint(i,Radius,source_pOmega);
  h_Source_pOmega_renorm->SetBinContent(i+1,source_pOmega);

  source_pKp = AB_pKpI1.EvaluateTheSource(0,Radius,0);
  Graph_Source_pKp.SetPoint(i,Radius,source_pKp);
  h_Source_pKp_renorm->SetBinContent(i+1,source_pKp);
  }
  double norm_source_pXi = h_Source_pXi_renorm->Integral(rmin,rmax);
  double norm_source_pOmega = h_Source_pOmega_renorm->Integral(rmin,rmax);
  double norm_source_pKp = h_Source_pKp_renorm->Integral(rmin,rmax);
  h_Source_pXi_renorm->Sumw2("off");
  h_Source_pXi_renorm->Scale(1./norm_source_pXi);
  h_Source_pOmega_renorm->Sumw2("off");
  h_Source_pOmega_renorm->Scale(1./norm_source_pOmega);
  h_Source_pKp_renorm->Sumw2("off");
  h_Source_pKp_renorm->Scale(1./norm_source_pKp);
//Filling the renormalized TGraphs
  for (unsigned i = 0; i < count; i++) {
    Radius = rmin+step*double(i);
    source_pXi = AB_pXi.EvaluateTheSource(0,Radius,0);
    Graph_Source_pXi_renorm.SetPoint(i,Radius,source_pXi/norm_source_pOmega);
    source_pOmega = AB_pOmega.EvaluateTheSource(0,Radius,0);
    Graph_Source_pOmega_renorm.SetPoint(i,Radius,source_pOmega/norm_source_pOmega);
    source_pKp = AB_pKpI1.EvaluateTheSource(0,Radius,0);
    Graph_Source_pKp_renorm.SetPoint(i,Radius,source_pKp/norm_source_pOmega);
  }

  Graph_Source_pXi.Write();
  Graph_Source_pXi_renorm.Write();

  Graph_Source_pOmega.Write();
  Graph_Source_pOmega_renorm.Write();

  Graph_Source_pKp.Write();
  Graph_Source_pKp_renorm.Write();

  OutputFileSource->Close();


// bins
// k*=10 -> bin=5
// k*=100 -> bin = 50
  std::vector<int> momvec={5,50};

  double dummy_pXi1 = 0.;
  double dummy_pXi2 = 0.;
  double dummy_pOmega1 = 0.;
  double dummy_pOmega2 = 0.;
  double dummy_pKp1 = 0.;
  double dummy_pKp2 = 0.;

  TH1F* histotestpXi1 = new TH1F("histotestpXi1","histotestpXi1",200,0.05,10.05);
  TH1F* histotestpXi2 = new TH1F("histotestpXi2","histotestpXi2",200,0.05,10.05);

  TH1F* histotestpOmega1 = new TH1F("histotestpOmega1","histotestpOmega1",200,0.05,10.05);
  TH1F* histotestpOmega2 = new TH1F("histotestpOmega2","histotestpOmega2",200,0.05,10.05);

  TH1F* histotestpKp1 = new TH1F("histotestpKp1","histotestpKp1",200,0.05,10.05);
  TH1F* histotestpKp2 = new TH1F("histotestpKp2","histotestpKp2",200,0.05,10.05);

//FIlling histograms for normalization
    for (int i = 0; i < count; i++) {
    Radius = rmin+step*double(i);
    dummy_pXi1 = AB_pXi.EvalWaveFun2(5,Radius,0);//momBin, radius channel
    dummy_pXi2 = AB_pXi.EvalWaveFun2(50,Radius,0);//momBin, radius channel

    dummy_pOmega1 = AB_pOmega.EvalWaveFun2(5,Radius,0);//momBin, radius channel
    dummy_pOmega2 = AB_pOmega.EvalWaveFun2(50,Radius,0);//momBin, radius channel

    dummy_pKp1 = AB_pKpI1.EvalWaveFun2(5,Radius,0);//momBin, radius channel
    dummy_pKp2 = AB_pKpI1.EvalWaveFun2(50,Radius,0);//momBin, radius channel

    histotestpXi1->SetBinContent(i+1,dummy_pXi1);
    histotestpXi2->SetBinContent(i+1,dummy_pXi2);
    histotestpOmega1->SetBinContent(i+1,dummy_pOmega1);
    histotestpOmega2->SetBinContent(i+1,dummy_pOmega2);
    histotestpKp1->SetBinContent(i+1,dummy_pKp1);
    histotestpKp2->SetBinContent(i+1,dummy_pKp2);
    }

  double norm_WF2_pXi_5 = (histotestpXi1->Integral(rmin,rmax));
  double norm_WF2_pXi_50 = (histotestpXi2->Integral(rmin,rmax));
  double norm_WF2_pOmega_50 = (histotestpOmega2->Integral(rmin,rmax));
  double norm_WF2_pKp_5 = (histotestpKp1->Integral(rmin,rmax));



  OutputFileWF->cd();

  for (int k = 0; k<momvec.size(); k++){
    int mombinsel=momvec[k];
    printf("Momentum Selected = %.2f\n",Ck_pXi->GetBinCenter(0,mombinsel));
    PlotGraph(Graph_WF2_pXi_I0S0[k],TString::Format("Graph_WF2_pXi_I0S0%i",mombinsel),true,kOrange+7,4,1);
    PlotGraph(Graph_WF2_pXi_I0S1[k],TString::Format("Graph_WF2_pXi_I0S1%i",mombinsel),true,kOrange+7,4,1);
    PlotGraph(Graph_WF2_pXi_I1S0[k],TString::Format("Graph_WF2_pXi_I1S0%i",mombinsel),true,kOrange+7,4,1);
    PlotGraph(Graph_WF2_pXi_I1S1[k],TString::Format("Graph_WF2_pXi_I1S1%i",mombinsel),true,kOrange+7,4,1);

    PlotGraph(Graph_WF2_pXi_I0S0_renorm[k],TString::Format("Graph_WF2_pXi_I0S0_renorm%i",mombinsel),true,kOrange+7,4,1);

    PlotGraph(Graph_WF2_pOmega_5S2[k],TString::Format("Graph_WF2_pOmega_5S2%i",mombinsel),true,kPink+5,4,1);

    PlotGraph(Graph_WF2_pOmega_5S2_renorm[k],TString::Format("Graph_WF2_pOmega_5S2_renorm%i",mombinsel),true,kPink+5,4,1);

    PlotGraph(Graph_WF2_pKp_I1[k],TString::Format("Graph_WF2_pKp_I1%i",mombinsel),true,kAzure+7,4,1);
    PlotGraph(Graph_WF2_pKp_RamonaI0[k],TString::Format("Graph_WF2_pKp_RamonaI0%i",mombinsel),true,kAzure+7,4,1);
    PlotGraph(Graph_WF2_pKp_RamonaI1[k],TString::Format("Graph_WF2_pKp_RamonaI1%i",mombinsel),true,kAzure+7,4,1);

    PlotGraph(Graph_WF2_pKp_I1_renorm[k],TString::Format("Graph_WF2_pKp_I1_renorm%i",mombinsel),true,kAzure+7,4,1);

  for (int i = 0; i < count; i++) {
  Radius = rmin+step*double(i);
  WF2_pXi_I0S0[k] = AB_pXi.EvalWaveFun2(mombinsel,Radius,0);//momBin, radius channel
  WF2_pXi_I0S1[k] = AB_pXi.EvalWaveFun2(mombinsel,Radius,1);//momBin, radius channel
  WF2_pXi_I1S0[k] = AB_pXi.EvalWaveFun2(mombinsel,Radius,2);//momBin, radius channel
  WF2_pXi_I1S1[k] = AB_pXi.EvalWaveFun2(mombinsel,Radius,3);//momBin, radius channel

  WF2_pOmega_5S2[k] = AB_pOmega.EvalWaveFun2(mombinsel,Radius,0);//momBin, radius channel

  WF2_pKp_I1[k] = AB_pKpI1.EvalWaveFun2(mombinsel,Radius,0);//momBin, radius channel
  WF2_pKp_RamonaI0[k] = AB_pKpRamona.EvalWaveFun2(mombinsel,Radius,0);//momBin, radius channel
  WF2_pKp_RamonaI1[k] = AB_pKpRamona.EvalWaveFun2(mombinsel,Radius,1);//momBin, radius channel

  Graph_WF2_pXi_I0S0[k].SetPoint(i,Radius,WF2_pXi_I0S0[k]);
  Graph_WF2_pXi_I0S1[k].SetPoint(i,Radius,WF2_pXi_I0S1[k]);
  Graph_WF2_pXi_I1S0[k].SetPoint(i,Radius,WF2_pXi_I1S0[k]);
  Graph_WF2_pXi_I1S1[k].SetPoint(i,Radius,WF2_pXi_I1S1[k]);
  Graph_WF2_pXi_I0S0_renorm[k].SetPoint(i,Radius,WF2_pXi_I0S0[k]);
  Graph_WF2_pXi_I0S1_renorm[k].SetPoint(i,Radius,WF2_pXi_I0S1[k]);
  Graph_WF2_pXi_I1S0_renorm[k].SetPoint(i,Radius,WF2_pXi_I1S0[k]);
  Graph_WF2_pXi_I1S1_renorm[k].SetPoint(i,Radius,WF2_pXi_I1S1[k]);

  Graph_WF2_pOmega_5S2[k].SetPoint(i,Radius,WF2_pOmega_5S2[k]);
  Graph_WF2_pOmega_5S2_renorm[k].SetPoint(i,Radius,WF2_pOmega_5S2[k]);

  Graph_WF2_pKp_I1[k].SetPoint(i,Radius,WF2_pKp_I1[k]);
  Graph_WF2_pKp_I1_renorm[k].SetPoint(i,Radius,WF2_pKp_I1[k]);
  Graph_WF2_pKp_RamonaI0[k].SetPoint(i,Radius,WF2_pKp_RamonaI0[k]);
  Graph_WF2_pKp_RamonaI1[k].SetPoint(i,Radius,WF2_pKp_RamonaI1[k]);


  }

  if(k==0){
    for (int i=0;i<Graph_WF2_pXi_I0S0_renorm[0].GetN();i++) Graph_WF2_pXi_I0S0_renorm[0].GetY()[i] *= (1./norm_WF2_pXi_5);
    for (int i=0;i<Graph_WF2_pOmega_5S2_renorm[0].GetN();i++) Graph_WF2_pOmega_5S2_renorm[0].GetY()[i] *= (1./norm_WF2_pXi_5);
    for (int i=0;i<Graph_WF2_pKp_I1_renorm[0].GetN();i++) Graph_WF2_pKp_I1_renorm[0].GetY()[i] *= (1./norm_WF2_pXi_5);
  } else if(k==1){
    for (int i=0;i<Graph_WF2_pXi_I0S0_renorm[1].GetN();i++) Graph_WF2_pXi_I0S0_renorm[1].GetY()[i] *= (1./norm_WF2_pXi_50);
    for (int i=0;i<Graph_WF2_pOmega_5S2_renorm[1].GetN();i++) Graph_WF2_pOmega_5S2_renorm[1].GetY()[i] *= (1./norm_WF2_pXi_50);
    for (int i=0;i<Graph_WF2_pKp_I1_renorm[1].GetN();i++) Graph_WF2_pKp_I1_renorm[1].GetY()[i] *= (1./norm_WF2_pXi_50);
  }

    Graph_WF2_pXi_I0S0[k].Write();
    Graph_WF2_pXi_I0S1[k].Write();
    Graph_WF2_pXi_I1S0[k].Write();
    Graph_WF2_pXi_I1S1[k].Write();

    Graph_WF2_pXi_I0S0_renorm[k].Write();

    Graph_WF2_pOmega_5S2[k].Write();

    Graph_WF2_pOmega_5S2_renorm[k].Write();

    Graph_WF2_pKp_I1[k].Write();
    Graph_WF2_pKp_RamonaI0[k].Write();
    Graph_WF2_pKp_RamonaI1[k].Write();

    Graph_WF2_pKp_I1_renorm[k].Write();
  }



  OutputFileWF->Close();

  OutputFileCF->cd();
  Graph_Ck_pXi.Write();
  Graph_Ck_pOmega.Write();
  Graph_Ck_pKpI1.Write();
  Graph_Ck_pKpRamona.Write();
  OutputFileCF->Close();

  delete OutputFileWF;
  delete OutputFileCF;
  delete OutputFileSource;


}




void ReviewFemto_Plots(){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);

  SetStyle(false,false);


  TFile* file_WF = TFile::Open("fOutput_WF.root");
  TFile* file_Source = TFile::Open("fOutput_Source.root");
  TFile* file_CF = TFile::Open("fOutput_CF.root");

  TGraph* Graph_Source_pXi = (TGraph*)(file_Source->Get("Graph_Source_pXi"));
  TGraph* Graph_Source_pOmega = (TGraph*)(file_Source->Get("Graph_Source_pOmega"));
  TGraph* Graph_Source_pKp = (TGraph*)(file_Source->Get("Graph_Source_pKp"));

  //|WF|^2 per k*=10
  TGraph* Graph_WF2_pXi_I0S05 = (TGraph*)(file_WF->Get("Graph_WF2_pXi_I0S0_renorm5"));
  TGraph* Graph_WF2_pOmega_5S25 = (TGraph*)(file_WF->Get("Graph_WF2_pOmega_5S2_renorm5"));
  TGraph* Graph_WF2_pKp_I15 = (TGraph*)(file_WF->Get("Graph_WF2_pKp_I1_renorm5"));

  //|WF|^2 per k*=100
  TGraph* Graph_WF2_pXi_I0S050 = (TGraph*)(file_WF->Get("Graph_WF2_pXi_I0S0_renorm50"));
  TGraph* Graph_WF2_pOmega_5S250 = (TGraph*)(file_WF->Get("Graph_WF2_pOmega_5S2_renorm50"));
  TGraph* Graph_WF2_pKp_I150 = (TGraph*)(file_WF->Get("Graph_WF2_pKp_I1_renorm50"));


  //Test trial plotting the wfs for pXi and Source for pXi
  auto* canvas_1 = new TCanvas("canvas_1","canvas_1",1200,675);
  canvas_1->Divide(2,1);
  canvas_1->cd(1);
  // canvas_1->SetRightMargin(2.0);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.0);
  //Drawing the wfs for pXi, k*=10
  Graph_WF2_pXi_I0S05->GetYaxis()->SetTitleOffset(1.7);
  Graph_WF2_pXi_I0S05->GetYaxis()->SetTitle("|#Psi_{k*}(r)|^{2} (a.u.)");
  Graph_WF2_pXi_I0S05->GetXaxis()->SetTitle("r (fm)");
  Graph_WF2_pXi_I0S05->GetXaxis()->SetRangeUser(0.,7.9);
  // for (int i=0;i<Graph_WF2_pKp_I15->GetN();i++) Graph_WF2_pXi_I0S05->GetY()[i] *= 1./10.;
  // for (int i=0;i<Graph_WF2_pKp_I15->GetN();i++) Graph_WF2_pOmega_5S25->GetY()[i] *= 1./10.;
  for (int i=0;i<Graph_WF2_pKp_I15->GetN();i++) Graph_WF2_pKp_I15->GetY()[i] *= 100.;
  Graph_WF2_pXi_I0S05->GetYaxis()->SetRangeUser(0.,1.1*TMath::MaxElement(Graph_WF2_pXi_I0S05->GetN(),Graph_WF2_pXi_I0S05->GetY()));
  Graph_WF2_pXi_I0S05->Draw();
  Graph_WF2_pOmega_5S25->Draw("same");
  Graph_WF2_pKp_I15->Draw("same");

  TH1F* hdummy1 = new TH1F("","",1,0,1);
  hdummy1->SetLineColor(kBlack);
  hdummy1->SetLineWidth(2);
  hdummy1->SetLineStyle(2);

  TH1F* hdummy2 = new TH1F("","",1,0,1);
  hdummy2->SetLineColor(kBlack);
  hdummy2->SetLineWidth(2);
  hdummy2->SetLineStyle(1);

  TLegend* legend = new TLegend(0.5, 0.7,0.99, 0.9);
  legend->SetFillStyle(0);
  legend->SetLineWidth(0);
  legend->SetNColumns(2);
  legend->SetColumnSeparation(0.01);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(gStyle->GetTextSize()*0.90);
  legend->AddEntry(hdummy1, "4#pi r^{2}S(r)", "l");
  legend->AddEntry(hdummy2, "| #Psi_{k*}(r)|^{ 2}", "l");
  legend->Draw("same");


  TLegend* legend1 = new TLegend(0.8, 0.55, 0.99, 0.73);
  legend1->SetFillStyle(0);
  legend1->SetLineWidth(0);
  legend1->SetBorderSize(0);
  legend1->SetTextFont(42);
  legend1->SetTextSize(gStyle->GetTextSize()*0.90);
  legend1->AddEntry(Graph_WF2_pXi_I0S05,"p-#Xi^{-}","lp");
  legend1->AddEntry(Graph_WF2_pOmega_5S25,"p-#Omega^{-}","lp");
  legend1->AddEntry(Graph_WF2_pKp_I15,"p-K^{+}","lp");
  legend1->Draw("same");

  TLatex BeamText;
  BeamText.SetTextFont(42);
  BeamText.SetTextSize(gStyle->GetTextSize()*0.85);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.7, 0.3, "k* = 10 MeV/c");


  canvas_1->Update();
  //scale the source to the pad coordinates
   double gr_max = TMath::MaxElement(Graph_Source_pXi->GetN(),Graph_Source_pXi->GetY());
   Float_t rightmax = 1.1*gr_max;
   Float_t scale = gPad->GetUymax()/rightmax;
   printf("gr_max = %.2f\n",gr_max);
   printf("scale = %.2f\n",scale);
   printf("rightmax = %.2f\n",rightmax);

   for (int i=0;i<Graph_Source_pXi->GetN();i++) Graph_Source_pXi->GetY()[i] *= scale;
   for (int i=0;i<Graph_Source_pOmega->GetN();i++) Graph_Source_pOmega->GetY()[i] *= scale;
   for (int i=0;i<Graph_Source_pKp->GetN();i++) Graph_Source_pKp->GetY()[i] *= scale;




   Graph_Source_pXi->Draw("same");
   Graph_Source_pOmega->Draw("same");
   Graph_Source_pKp->Draw("same");

   // draw an axis on the right side
   printf("gPad->GetUxmax() = %.2f\n",gPad->GetUxmax());
   printf("gPad->GetUymin() = %.2f\n",gPad->GetUymin());
   printf("gPad->GetUymax() = %.2f\n",gPad->GetUymax());


   TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
   gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
    axis->SetTextFont(42);
    axis->SetTitleSize(0.05);
    axis->SetTitleOffset(1.25);
    axis->SetLabelFont(42);
    axis->SetLabelSize(0.0);
    axis->SetLabelOffset(0.01);
    // axis->SetTitle("4#pi r^{2} S (r) (1/fm)");
    axis->Draw();

  // canvas_1->Print("canvas_1.pdf");


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Same stuff for k*=100 MeV
  // auto* canvas_2 = new TCanvas("canvas_2","canvas_2",1000,800);
  // canvas_2->SetRightMargin(0.2);
  canvas_1->cd(2);
  gPad->SetLeftMargin(0.04);
  gPad->SetRightMargin(0.16);
  // Graph_WF2_pXi_I0S050->GetYaxis()->SetTitle("|#Psi_{k*}(r)|^{2} (a.u.)");
  Graph_WF2_pXi_I0S050->GetXaxis()->SetTitle("r (fm)");
  Graph_WF2_pXi_I0S050->GetXaxis()->SetRangeUser(0.,7.9);
  Graph_WF2_pXi_I0S050->GetYaxis()->SetLabelSize(0.0);
  for (int i=0;i<Graph_WF2_pKp_I150->GetN();i++) Graph_WF2_pXi_I0S050->GetY()[i] *= 0.5;
  for (int i=0;i<Graph_WF2_pKp_I150->GetN();i++) Graph_WF2_pOmega_5S250->GetY()[i] *= 0.5;
  for (int i=0;i<Graph_WF2_pKp_I150->GetN();i++) Graph_WF2_pKp_I150->GetY()[i] *= 0.5;
  Graph_WF2_pXi_I0S050->GetYaxis()->SetRangeUser(0.,1.1*TMath::MaxElement(Graph_WF2_pXi_I0S05->GetN(),Graph_WF2_pXi_I0S05->GetY()));
  Graph_WF2_pXi_I0S050->Draw();
  Graph_WF2_pOmega_5S250->Draw("same");
  Graph_WF2_pKp_I150->Draw("same");
  TLatex BeamText2;
  BeamText2.SetTextFont(42);
  BeamText2.SetTextSize(gStyle->GetTextSize()*0.85);
  BeamText2.SetNDC(kTRUE);
  BeamText2.DrawLatex(0.53, 0.3, "k* = 100 MeV/c");
  canvas_1->Update();

  // TH1F* hdummy1 = new TH1F("","",1,0,1);
  // hdummy1->SetLineColor(kBlack);
  // hdummy1->SetLineWidth(4);
  // hdummy1->SetLineStyle(7);

  // TH1F* hdummy2 = new TH1F("","",1,0,1);
  // hdummy2->SetLineColor(kBlack);
  // hdummy2->SetLineWidth(4);
  // hdummy2->SetLineStyle(1);

  // TLegend* legend = new TLegend(0.4, 0.7, 0.77, 0.9);
  // legend->SetFillStyle(0);
  // legend->SetLineWidth(0);
  // legend->SetNColumns(2);
  // legend->SetBorderSize(0);
  // legend->SetTextFont(42);
  // legend->SetTextSize(gStyle->GetTextSize()*0.90);
  // legend->AddEntry(hdummy1, "4#pi r^{2}S(r)", "l");
  // legend->AddEntry(hdummy2, "| #Psi_{k*}(r)|^{ 2}", "l");
  // legend->Draw("same");


  // TLegend* legend1 = new TLegend(0.65, 0.55, 0.77, 0.73);
  // legend1->SetFillStyle(0);
  // legend1->SetLineWidth(0);
  // legend1->SetBorderSize(0);
  // legend1->SetTextFont(42);
  // legend1->SetTextSize(gStyle->GetTextSize()*0.90);
  // legend1->AddEntry(Graph_WF2_pXi_I0S05,"p-#Xi^{-}","lp");
  // legend1->AddEntry(Graph_WF2_pOmega_5S25,"p-#Omega^{-}","lp");
  // legend1->AddEntry(Graph_WF2_pKp_I15,"p-K^{+}","lp");
  // legend1->Draw("same");

  //scale the source to the pad coordinates
  // double
  gr_max = TMath::MaxElement(Graph_Source_pXi->GetN(),Graph_Source_pXi->GetY());
  // //  Float_t
  rightmax = 1.1*gr_max;
  // //  Float_t
  scale = gPad->GetUymax()/rightmax;

  for (int i=0;i<Graph_Source_pXi->GetN();i++) Graph_Source_pXi->GetY()[i] *= scale;
  for (int i=0;i<Graph_Source_pOmega->GetN();i++) Graph_Source_pOmega->GetY()[i] *= scale;
  for (int i=0;i<Graph_Source_pKp->GetN();i++) Graph_Source_pKp->GetY()[i] *= scale;

  Graph_Source_pXi->Draw("same");
  Graph_Source_pOmega->Draw("same");
  Graph_Source_pKp->Draw("same");

  axis->SetTextFont(42);
  axis->SetTitleSize(0.05);
  axis->SetTitleOffset(1.7);
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.045);
  axis->SetLabelOffset(0.01);
  axis->SetTitle("4#pi r^{2} S (r) (1/fm)");
  axis->Draw();
  // canvas_1->Print("canvas_1.pdf");



}


void SetStyleAxis_2(TH1F* hAxis){
    hAxis->SetStats(false);
    hAxis->SetTitle("");

    hAxis->GetXaxis()->SetTitleSize(0.06);
    hAxis->GetXaxis()->SetLabelSize(0.06);
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);

    hAxis->GetYaxis()->SetTitleSize(0.06);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    hAxis->GetYaxis()->SetTitleOffset(0.95);
    hAxis->GetYaxis()->SetNdivisions(510);
}


void ReviewFemto_PlotsDimi(){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);

  SetStyle(false,false);

  TString InputFolder = "/Users/sartozza/Desktop/ReviewFemto/MaterialForPlots/";

  TFile* file = TFile::Open(InputFolder+"fOutput_Dimi.root");

  TF1* PotRep = (TF1*)(file->Get("fPot_Toy1_X"));
  TF1* PotAttr = (TF1*)(file->Get("fPot_ND54_X"));
  TF1* PotBS = (TF1*)(file->Get("fPot_ND46_X"));

  TGraph* WF2Rep = (TGraph*)(file->Get("gTotWF_Toy1_X"));
  TGraph* WF2Attr = (TGraph*)(file->Get("gTotWF_ND54_X"));
  TGraph* WF2BS = (TGraph*)(file->Get("gTotWF_ND46_X"));

  TGraph* Source1 = (TGraph*)(file->Get("gSource_Toy1_X"));
  TGraph* Source15 = (TGraph*)(file->Get("gSource_ND54_X"));
  TGraph* Source4 = (TGraph*)(file->Get("gSource_ND46_X"));

  TGraph* CkAttrSource1 = (TGraph*)(file->Get("gCk_ND54_X_1.0"));
  TGraph* CkAttrSource15 = (TGraph*)(file->Get("gCk_ND54_X_1.5"));
  TGraph* CkAttrSource4 = (TGraph*)(file->Get("gCk_ND54_X_4.0"));

  TGraph* CkRepSource1 = (TGraph*)(file->Get("gCk_Toy1_X_1.0"));
  TGraph* CkRepSource15 = (TGraph*)(file->Get("gCk_Toy1_X_1.5"));
  TGraph* CkRepSource4 = (TGraph*)(file->Get("gCk_Toy1_X_4.0"));

  TGraph* CkBSSource1 = (TGraph*)(file->Get("gCk_ND46_X_1.0"));
  TGraph* CkBSSource15 = (TGraph*)(file->Get("gCk_ND46_X_1.5"));
  TGraph* CkBSSource4 = (TGraph*)(file->Get("gCk_ND46_X_4.0"));

  PlotTF1(*PotRep,"PotRep",true,kAzure+7,4,1);
  PlotTF1(*PotAttr,"PotRep",true,kOrange+7,4,1);
  PlotTF1(*PotBS,"PotRep",true,kPink+5,4,1);

  PlotGraph(*WF2Rep,"WF2Rep",true,kAzure+7,4,1);
  PlotGraph(*WF2Attr,"WF2Attr",true,kOrange+7,4,1);
  PlotGraph(*WF2BS,"WF2BS",true,kPink+5,4,1);

  PlotGraph(*Source1,"Source1",true,kGray+3,2,9);
  PlotGraph(*Source15,"Source15",true,kGray+3,2,3);
  PlotGraph(*Source4,"Source4",true,kGray+3,2,4);

  PlotGraph(*CkAttrSource1,"CkAttrSource1",true,kOrange+7,3,9);
  PlotGraph(*CkAttrSource15,"CkAttrSource15",true,kOrange+7,3,3);
  PlotGraph(*CkAttrSource4,"CkAttrSource4",true,kOrange+7,3,4);

  PlotGraph(*CkRepSource1,"CkRepSource1",true,kAzure+7,3,9);
  PlotGraph(*CkRepSource15,"CkRepSource15",true,kAzure+7,3,3);
  PlotGraph(*CkRepSource4,"CkRepSource4",true,kAzure+7,3,4);

  PlotGraph(*CkBSSource1,"CkBSSource1",true,kPink+5,3,9);
  PlotGraph(*CkBSSource15,"CkBSSource15",true,kPink+5,3,3);
  PlotGraph(*CkBSSource4,"CkBSSource4",true,kPink+5,3,4);



  double WFu_yMax = 1.1*TMath::MaxElement(WF2Attr->GetN(),WF2Attr->GetY());
  double WFtot_yMax = 1.1*TMath::MaxElement(WF2Attr->GetN(),WF2Attr->GetY());
  double WF_rMax = 6.;
  double Kmax = 150.;

  TH1F* hPot_Dummy = new TH1F("hPot_Dummy","hPot_Dummy",128,0,WF_rMax);
  SetStyleAxis_2(hPot_Dummy);
  hPot_Dummy->GetXaxis()->SetTitle("r (fm)");
  hPot_Dummy->GetYaxis()->SetTitle("V(r) (MeV)");
  hPot_Dummy->GetYaxis()->SetRangeUser(-120,400);
  TH1F* hWFtot_Dummy = new TH1F("hWFtot_Dummy","hWFtot_Dummy",128,0,WF_rMax);
  SetStyleAxis_2(hWFtot_Dummy);
  hWFtot_Dummy->GetXaxis()->SetTitle("r (fm)");
  hWFtot_Dummy->GetYaxis()->SetTitle("|#Psi_{k}(r)|^{2}");
  hWFtot_Dummy->GetYaxis()->SetRangeUser(-0.1,WFtot_yMax);
  TH1F* hSource_Dummy = new TH1F("hSource_Dummy","hSource_Dummy",128,0,WF_rMax);
  SetStyleAxis_2(hSource_Dummy);
  hSource_Dummy->GetXaxis()->SetTitle("r (fm)");
  hSource_Dummy->GetYaxis()->SetTitle("|#Psi_{k*}(r)|^{2}");
  hSource_Dummy->GetYaxis()->SetRangeUser(-0.1,WFtot_yMax);
  TH1F* hCk_Dummy = new TH1F("hCk_Dummy","hCk_Dummy",128,0,Kmax);
  SetStyleAxis_2(hCk_Dummy);
  hCk_Dummy->GetXaxis()->SetTitle("k* (MeV/c)");
  hCk_Dummy->GetYaxis()->SetTitle("C(k*)");
  hCk_Dummy->GetYaxis()->SetRangeUser(0.,3.2);

  TF1* fWFtot_base = new TF1("fWFtot_base","1",0,WF_rMax);
  fWFtot_base->SetLineColor(kBlack);
  fWFtot_base->SetLineWidth(1);
  fWFtot_base->SetLineStyle(2);

  TF1* fPot_base = new TF1("fPot_base","0",0,WF_rMax);
  fPot_base->SetLineColor(kBlack);
  fPot_base->SetLineWidth(1);
  fPot_base->SetLineStyle(2);

  TF1* fCk_base = new TF1("fCk_base","1",0,Kmax);
  fCk_base->SetLineColor(kBlack);
  fCk_base->SetLineWidth(1);
  fCk_base->SetLineStyle(2);

  //Starting the 3 Subpads mode
  DLM_SubPads Pad_WF(1080,1440);
  Pad_WF.AddSubPad(0,1,0.7,1);//lrbt
  Pad_WF.AddSubPad(0,1,0.4,0.7);
  Pad_WF.AddSubPad(0,1,0,0.37);

  Pad_WF.SetMargin(0,0.15,0.15,0.0,0.06);
  Pad_WF.SetMargin(1,0.15,0.15,0.05,0.);
  Pad_WF.SetMargin(2,0.15,0.15,0.06,0.01);

//First pad with potentials
  Pad_WF.cd(0);
  hPot_Dummy->Draw("axis");
  fPot_base->Draw("same");
  PotRep->Draw("same");
  PotAttr->Draw("same");
  PotBS->Draw("same");
  TLegend* legend1 = new TLegend(0.45, 0.25, 0.83, 0.73);
  legend1->SetFillStyle(0);
  legend1->SetLineWidth(0);
  legend1->SetBorderSize(0);
  legend1->SetTextFont(42);
  legend1->SetTextSize(gStyle->GetTextSize()*0.90);
  legend1->AddEntry(PotRep,"Repulsive","l");
  legend1->AddEntry(PotAttr,"Attractive","l");
  legend1->AddEntry(PotBS,"with a bound State","l");
  legend1->Draw("same");

//Second Pad with |wf|^2 e Sources
  Pad_WF.cd(1);
  Pad_WF.GetCanvas()->SetFillColor(0);
  Pad_WF.GetCanvas()->SetFillStyle(0);
  Pad_WF.GetCanvas()->SetFrameFillStyle(0);
  Pad_WF.GetCanvas()->Update();

  hWFtot_Dummy->Draw("axis");
  fWFtot_base->Draw("same");
  WF2Attr->Draw("same");
  WF2Rep->Draw("same");
  WF2BS->Draw("same");

  TH1F* hdummy1 = new TH1F("","",1,0,1);
  hdummy1->SetLineColor(kBlack);
  hdummy1->SetLineWidth(2);
  hdummy1->SetLineStyle(2);

  TH1F* hdummy2 = new TH1F("","",1,0,1);
  hdummy2->SetLineColor(kBlack);
  hdummy2->SetLineWidth(2);
  hdummy2->SetLineStyle(1);

  TLegend* legend = new TLegend(0.5, 0.8,0.83, 0.9);
  legend->SetFillStyle(0);
  legend->SetLineWidth(0);
  legend->SetNColumns(2);
  legend->SetColumnSeparation(0.);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(gStyle->GetTextSize()*0.90);
  legend->AddEntry(hdummy2, "| #Psi_{k*}(r)|^{ 2}", "l");
  legend->AddEntry(hdummy1, "4#pi r^{2}S(r)", "l");
  legend->Draw("same");


  //  Pad_WF.GetCanvas();
   Source1->Draw("same");
   Source15->Draw("same");
   Source4->Draw("same");

  TLegend* legendSource = new TLegend(0.57, 0.65, 0.96, 0.8);
  legendSource->SetFillStyle(0);
  legendSource->SetLineWidth(0);
  legendSource->SetBorderSize(0);
  legendSource->SetTextFont(42);
  legendSource->SetTextSize(gStyle->GetTextSize()*0.90);
  legendSource->AddEntry(Source1,"r_{G} = 1 fm","l");
  legendSource->AddEntry(Source15,"r_{G} = 1.5 fm","l");
  legendSource->AddEntry(Source4,"r_{G} = 4 fm","l");
  legendSource->Draw("same");

//scale the source to the pad coordinates
  double gr_max = TMath::MaxElement(Source1->GetN(),Source1->GetY());
  printf("gr_max = %.2f",gr_max);
  Float_t rightmax = 1.2*gr_max;
  printf("rightmax = %.2f",rightmax);
  Float_t scale = WFu_yMax/rightmax;

  for (int i=0;i<Source1->GetN();i++) Source1->GetY()[i] *= scale;
  for (int i=0;i<Source15->GetN();i++) Source15->GetY()[i] *= scale;
  for (int i=0;i<Source4->GetN();i++) Source4->GetY()[i] *= scale;

    // TGaxis* Arho = new TGaxis(0,WFu_yMax,WF_rMax,WFu_yMax,0,WF_rhoMax,505,"-");
    // Arho->SetTitle("#rho=kr");
    // Arho->SetLabelSize(hWFu_Dummy->GetXaxis()->GetLabelSize());
    // Arho->SetTitleSize(hWFu_Dummy->GetXaxis()->GetTitleSize());
    // Arho->SetLabelFont(hWFu_Dummy->GetXaxis()->GetLabelFont());
    // Arho->SetTitleFont(hWFu_Dummy->GetXaxis()->GetTitleFont());
    // Arho->Draw("same");
   // draw an axis on the right side
  TGaxis *axis = new TGaxis(WF_rMax,0,WF_rMax,WFu_yMax,0.,rightmax,505,"+L");
  // axis->SetTextFont(42);
  // axis->SetTitleSize(0.05);
  // axis->SetTitleOffset(0.7);
  // axis->SetLabelFont(42);
  // axis->SetLabelSize(0.045);
  // axis->SetLabelOffset(0.01);
  axis->SetLabelSize(hWFtot_Dummy->GetXaxis()->GetLabelSize());
  axis->SetTitleSize(hWFtot_Dummy->GetXaxis()->GetTitleSize());
  axis->SetLabelFont(hWFtot_Dummy->GetXaxis()->GetLabelFont());
  axis->SetTitleFont(hWFtot_Dummy->GetXaxis()->GetTitleFont());
  axis->SetTitle("4#pi r^{2} S (r) (1/fm)");
  axis->Draw();


//Subpad for correlations
  Pad_WF.cd(2);
  Pad_WF.GetCanvas()->SetFillColor(0);
  Pad_WF.GetCanvas()->SetFrameFillStyle(0);
  Pad_WF.GetCanvas()->Update();
  hCk_Dummy->Draw("axis");
  fCk_base->Draw("same");
  CkAttrSource1->Draw("same");
  CkAttrSource15->Draw("same");
  CkAttrSource4->Draw("same");
  CkRepSource1->Draw("same");
  CkRepSource15->Draw("same");
  CkRepSource4->Draw("same");
  CkBSSource1->Draw("same");
  CkBSSource15->Draw("same");
  CkBSSource4->Draw("same");

  TH1F* hdummy1C = new TH1F("","",1,0,1);
  hdummy1C->SetLineColor(kBlack);
  hdummy1C->SetLineWidth(2);
  hdummy1C->SetLineStyle(9);

  TH1F* hdummy15C = new TH1F("","",1,0,1);
  hdummy15C->SetLineColor(kBlack);
  hdummy15C->SetLineWidth(2);
  hdummy15C->SetLineStyle(3);

  TH1F* hdummy4C = new TH1F("","",1,0,1);
  hdummy4C->SetLineColor(kBlack);
  hdummy4C->SetLineWidth(2);
  hdummy4C->SetLineStyle(4);

  TLegend* legendCk = new TLegend(0.5, 0.8,0.83, 0.9);
  legendCk->SetFillStyle(0);
  legendCk->SetLineWidth(0);
  legendCk->SetBorderSize(0);
  legendCk->SetTextFont(42);
  legendCk->SetTextSize(gStyle->GetTextSize()*0.90);
  legendCk->AddEntry(hdummy1C, "r_{G} = 1 fm", "l");
  legendCk->AddEntry(hdummy15C, "r_{G} = 1.5 fm", "l");
  legendCk->AddEntry(hdummy4C, "r_{G} = 4 fm", "l");
  legendSource->Draw("same");



  TString Outputfolder = gSystem->pwd();
  Pad_WF.GetCanvas()->SaveAs(Outputfolder+"/Dimi_Canvas.pdf");

/*
  TGraph* Graph_Source_pXi = (TGraph*)(file_Source->Get("Graph_Source_pXi"));
  TGraph* Graph_Source_pOmega = (TGraph*)(file_Source->Get("Graph_Source_pOmega"));
  TGraph* Graph_Source_pKp = (TGraph*)(file_Source->Get("Graph_Source_pKp"));

  //|WF|^2 per k*=10
  TGraph* Graph_WF2_pXi_I0S05 = (TGraph*)(file_WF->Get("Graph_WF2_pXi_I0S0_renorm5"));
  TGraph* Graph_WF2_pOmega_5S25 = (TGraph*)(file_WF->Get("Graph_WF2_pOmega_5S2_renorm5"));
  TGraph* Graph_WF2_pKp_I15 = (TGraph*)(file_WF->Get("Graph_WF2_pKp_I1_renorm5"));

  //|WF|^2 per k*=100
  TGraph* Graph_WF2_pXi_I0S050 = (TGraph*)(file_WF->Get("Graph_WF2_pXi_I0S0_renorm50"));
  TGraph* Graph_WF2_pOmega_5S250 = (TGraph*)(file_WF->Get("Graph_WF2_pOmega_5S2_renorm50"));
  TGraph* Graph_WF2_pKp_I150 = (TGraph*)(file_WF->Get("Graph_WF2_pKp_I1_renorm50"));


  //Test trial plotting the wfs for pXi and Source for pXi
  auto* canvas_1 = new TCanvas("canvas_1","canvas_1",1200,675);
  canvas_1->Divide(2,1);
  canvas_1->cd(1);
  // canvas_1->SetRightMargin(2.0);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.0);
  //Drawing the wfs for pXi, k*=10
  Graph_WF2_pXi_I0S05->GetYaxis()->SetTitleOffset(1.7);
  Graph_WF2_pXi_I0S05->GetYaxis()->SetTitle("|#Psi_{k*}(r)|^{2} (a.u.)");
  Graph_WF2_pXi_I0S05->GetXaxis()->SetTitle("r (fm)");
  Graph_WF2_pXi_I0S05->GetXaxis()->SetRangeUser(0.,7.9);
  // for (int i=0;i<Graph_WF2_pKp_I15->GetN();i++) Graph_WF2_pXi_I0S05->GetY()[i] *= 1./10.;
  // for (int i=0;i<Graph_WF2_pKp_I15->GetN();i++) Graph_WF2_pOmega_5S25->GetY()[i] *= 1./10.;
  for (int i=0;i<Graph_WF2_pKp_I15->GetN();i++) Graph_WF2_pKp_I15->GetY()[i] *= 100.;
  Graph_WF2_pXi_I0S05->GetYaxis()->SetRangeUser(0.,1.1*TMath::MaxElement(Graph_WF2_pXi_I0S05->GetN(),Graph_WF2_pXi_I0S05->GetY()));
  Graph_WF2_pXi_I0S05->Draw();
  Graph_WF2_pOmega_5S25->Draw("same");
  Graph_WF2_pKp_I15->Draw("same");

  TH1F* hdummy1 = new TH1F("","",1,0,1);
  hdummy1->SetLineColor(kBlack);
  hdummy1->SetLineWidth(2);
  hdummy1->SetLineStyle(2);

  TH1F* hdummy2 = new TH1F("","",1,0,1);
  hdummy2->SetLineColor(kBlack);
  hdummy2->SetLineWidth(2);
  hdummy2->SetLineStyle(1);

  TLegend* legend = new TLegend(0.5, 0.7,0.99, 0.9);
  legend->SetFillStyle(0);
  legend->SetLineWidth(0);
  legend->SetNColumns(2);
  legend->SetColumnSeparation(0.01);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(gStyle->GetTextSize()*0.90);
  legend->AddEntry(hdummy1, "4#pi r^{2}S(r)", "l");
  legend->AddEntry(hdummy2, "| #Psi_{k*}(r)|^{ 2}", "l");
  legend->Draw("same");


  TLegend* legend1 = new TLegend(0.8, 0.55, 0.99, 0.73);
  legend1->SetFillStyle(0);
  legend1->SetLineWidth(0);
  legend1->SetBorderSize(0);
  legend1->SetTextFont(42);
  legend1->SetTextSize(gStyle->GetTextSize()*0.90);
  legend1->AddEntry(Graph_WF2_pXi_I0S05,"p-#Xi^{-}","lp");
  legend1->AddEntry(Graph_WF2_pOmega_5S25,"p-#Omega^{-}","lp");
  legend1->AddEntry(Graph_WF2_pKp_I15,"p-K^{+}","lp");
  legend1->Draw("same");

  TLatex BeamText;
  BeamText.SetTextFont(42);
  BeamText.SetTextSize(gStyle->GetTextSize()*0.85);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.7, 0.3, "k* = 10 MeV/c");


  canvas_1->Update();
  //scale the source to the pad coordinates
  double gr_max = TMath::MaxElement(Graph_Source_pXi->GetN(),Graph_Source_pXi->GetY());
   Float_t rightmax = 1.1*gr_max;
   Float_t scale = gPad->GetUymax()/rightmax;
   printf("gr_max = %.2f\n",gr_max);
   printf("scale = %.2f\n",scale);
   printf("rightmax = %.2f\n",rightmax);

   for (int i=0;i<Graph_Source_pXi->GetN();i++) Graph_Source_pXi->GetY()[i] *= scale;
   for (int i=0;i<Graph_Source_pOmega->GetN();i++) Graph_Source_pOmega->GetY()[i] *= scale;
   for (int i=0;i<Graph_Source_pKp->GetN();i++) Graph_Source_pKp->GetY()[i] *= scale;




   Graph_Source_pXi->Draw("same");
   Graph_Source_pOmega->Draw("same");
   Graph_Source_pKp->Draw("same");

   // draw an axis on the right side
   printf("gPad->GetUxmax() = %.2f\n",gPad->GetUxmax());
   printf("gPad->GetUymin() = %.2f\n",gPad->GetUymin());
   printf("gPad->GetUymax() = %.2f\n",gPad->GetUymax());


   TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
   gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
    axis->SetTextFont(42);
    axis->SetTitleSize(0.05);
    axis->SetTitleOffset(1.25);
    axis->SetLabelFont(42);
    axis->SetLabelSize(0.0);
    axis->SetLabelOffset(0.01);
    // axis->SetTitle("4#pi r^{2} S (r) (1/fm)");
    axis->Draw();

  // canvas_1->Print("canvas_1.pdf");


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Same stuff for k*=100 MeV
  // auto* canvas_2 = new TCanvas("canvas_2","canvas_2",1000,800);
  // canvas_2->SetRightMargin(0.2);
  canvas_1->cd(2);
  gPad->SetLeftMargin(0.04);
  gPad->SetRightMargin(0.16);
  // Graph_WF2_pXi_I0S050->GetYaxis()->SetTitle("|#Psi_{k*}(r)|^{2} (a.u.)");
  Graph_WF2_pXi_I0S050->GetXaxis()->SetTitle("r (fm)");
  Graph_WF2_pXi_I0S050->GetXaxis()->SetRangeUser(0.,7.9);
  Graph_WF2_pXi_I0S050->GetYaxis()->SetLabelSize(0.0);
  for (int i=0;i<Graph_WF2_pKp_I150->GetN();i++) Graph_WF2_pXi_I0S050->GetY()[i] *= 0.5;
  for (int i=0;i<Graph_WF2_pKp_I150->GetN();i++) Graph_WF2_pOmega_5S250->GetY()[i] *= 0.5;
  for (int i=0;i<Graph_WF2_pKp_I150->GetN();i++) Graph_WF2_pKp_I150->GetY()[i] *= 0.5;
  Graph_WF2_pXi_I0S050->GetYaxis()->SetRangeUser(0.,1.1*TMath::MaxElement(Graph_WF2_pXi_I0S05->GetN(),Graph_WF2_pXi_I0S05->GetY()));
  Graph_WF2_pXi_I0S050->Draw();
  Graph_WF2_pOmega_5S250->Draw("same");
  Graph_WF2_pKp_I150->Draw("same");
  TLatex BeamText2;
  BeamText2.SetTextFont(42);
  BeamText2.SetTextSize(gStyle->GetTextSize()*0.85);
  BeamText2.SetNDC(kTRUE);
  BeamText2.DrawLatex(0.53, 0.3, "k* = 100 MeV/c");
  canvas_1->Update();

  // TH1F* hdummy1 = new TH1F("","",1,0,1);
  // hdummy1->SetLineColor(kBlack);
  // hdummy1->SetLineWidth(4);
  // hdummy1->SetLineStyle(7);

  // TH1F* hdummy2 = new TH1F("","",1,0,1);
  // hdummy2->SetLineColor(kBlack);
  // hdummy2->SetLineWidth(4);
  // hdummy2->SetLineStyle(1);

  // TLegend* legend = new TLegend(0.4, 0.7, 0.77, 0.9);
  // legend->SetFillStyle(0);
  // legend->SetLineWidth(0);
  // legend->SetNColumns(2);
  // legend->SetBorderSize(0);
  // legend->SetTextFont(42);
  // legend->SetTextSize(gStyle->GetTextSize()*0.90);
  // legend->AddEntry(hdummy1, "4#pi r^{2}S(r)", "l");
  // legend->AddEntry(hdummy2, "| #Psi_{k*}(r)|^{ 2}", "l");
  // legend->Draw("same");


  // TLegend* legend1 = new TLegend(0.65, 0.55, 0.77, 0.73);
  // legend1->SetFillStyle(0);
  // legend1->SetLineWidth(0);
  // legend1->SetBorderSize(0);
  // legend1->SetTextFont(42);
  // legend1->SetTextSize(gStyle->GetTextSize()*0.90);
  // legend1->AddEntry(Graph_WF2_pXi_I0S05,"p-#Xi^{-}","lp");
  // legend1->AddEntry(Graph_WF2_pOmega_5S25,"p-#Omega^{-}","lp");
  // legend1->AddEntry(Graph_WF2_pKp_I15,"p-K^{+}","lp");
  // legend1->Draw("same");

  //scale the source to the pad coordinates
  // double
  gr_max = TMath::MaxElement(Graph_Source_pXi->GetN(),Graph_Source_pXi->GetY());
  // //  Float_t
  rightmax = 1.1*gr_max;
  // //  Float_t
  scale = gPad->GetUymax()/rightmax;

  for (int i=0;i<Graph_Source_pXi->GetN();i++) Graph_Source_pXi->GetY()[i] *= scale;
  for (int i=0;i<Graph_Source_pOmega->GetN();i++) Graph_Source_pOmega->GetY()[i] *= scale;
  for (int i=0;i<Graph_Source_pKp->GetN();i++) Graph_Source_pKp->GetY()[i] *= scale;

  Graph_Source_pXi->Draw("same");
  Graph_Source_pOmega->Draw("same");
  Graph_Source_pKp->Draw("same");

  axis->SetTextFont(42);
  axis->SetTitleSize(0.05);
  axis->SetTitleOffset(1.7);
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.045);
  axis->SetLabelOffset(0.01);
  axis->SetTitle("4#pi r^{2} S (r) (1/fm)");
  axis->Draw();
  // canvas_1->Print("canvas_1.pdf");
  */
}

//########################################
int pAntiLambda(int argc, char *argv[]){
    printf("\033[0;33m mT differential Fits \033[0m\n");

  //Theo_pAntip();
  // ReviewFemto_Calculations();
  // ReviewFemto_Plots();
  ReviewFemto_PlotsDimi();

  return 0;
}

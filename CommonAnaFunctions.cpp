
#include "CommonAnaFunctions.h"
#include "CATS.h"
#include "CATStools.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_WfModel.h"
#include <iostream>
#include "TString.h"
#include "TH2F.h"
#include "TFile.h"
#include "TROOT.h"
#include "TGraph.h"


DLM_CommonAnaFunctions::DLM_CommonAnaFunctions():NumCleverLevyObjects(5){
    //Simple_Reso = NULL;
    //Simple_Reso = new MS_GaussExp_mT_Simple [NumCleverLevyObjects];
    CleverLevy = NULL;
    CleverLevy = new DLM_CleverLevy [NumCleverLevyObjects];
    CleverMcLevyReso = NULL;
    CleverMcLevyReso = new DLM_CleverMcLevyReso [NumCleverLevyObjects];
    CatsFilesFolder = new TString();
}

DLM_CommonAnaFunctions::~DLM_CommonAnaFunctions(){
    //if(Simple_Reso){delete[]Simple_Reso;Simple_Reso=NULL;}
    if(CleverLevy){delete[]CleverLevy;CleverLevy=NULL;}
    if(CleverMcLevyReso){delete[]CleverMcLevyReso;CleverMcLevyReso=NULL;}
    delete CatsFilesFolder;
}

//Implement Lednicky
void DLM_CommonAnaFunctions::SetUpCats_pAp(CATS& Kitty, const TString& POT, const TString& SOURCE){
  std::cout<<"SetUpCats_pAp dummy \n"<<std::endl;
}
// Proton-AntiProton ONLY COULOMB
void DLM_CommonAnaFunctions::SetUpCats_pApCoulomb(CATS& Kitty, const TString& POT, const TString& SOURCE,
		 const TString& DataSample){

    CATSparameters* cPars = NULL;
    double radius;
    if(DataSample=="pp13TeV_MB_BBar") radius=1.188;
    if(DataSample=="pp13TeV_HM_BBar") radius=1.28;

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,radius);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="McGauss_Reso"){
        CleverMcLevyReso[1].InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyReso[1].InitScale(35,0.25,2.0);
        CleverMcLevyReso[1].InitRad(256,0,64);
        CleverMcLevyReso[1].InitType(2);
        CleverMcLevyReso[1].InitReso(0,1);
        CleverMcLevyReso[1].InitReso(1,1);
        CleverMcLevyReso[1].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso[1].SetUpReso(1,0,1.-0.3562,1462.93,4.69,Mass_L,Mass_pic);
        CleverMcLevyReso[1].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[1], 2);
        Kitty.SetAnaSource(0,0.87);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        goto CLEAN_SetUpCats_pApCoulomb;
    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, -2212);
    Kitty.SetRedMass( (Mass_p*Mass_p)/(Mass_p+Mass_p) );

    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);
    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);
    Kitty.SetChannelWeight(0, 1./4.);
    Kitty.SetChannelWeight(1, 3./4.);

    CLEAN_SetUpCats_pApCoulomb: ;

}

void DLM_CommonAnaFunctions::SetUpCats_OmegaOmegaCoulomb(CATS& Kitty, const TString& POT, const TString& SOURCE,
		 const TString& DataSample){

    CATSparameters* cPars = NULL;
    double radius;
    if(DataSample=="pp13TeV_MB_BBar") radius=1.188;
    if(DataSample=="pp13TeV_HM_BBar") radius=1.28;

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,radius);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="McGauss_Reso"){
        CleverMcLevyReso[1].InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyReso[1].InitScale(35,0.25,2.0);
        CleverMcLevyReso[1].InitRad(256,0,64);
        CleverMcLevyReso[1].InitType(2);
        CleverMcLevyReso[1].InitReso(0,1);
        CleverMcLevyReso[1].InitReso(1,1);
        CleverMcLevyReso[1].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso[1].SetUpReso(1,0,1.-0.3562,1462.93,4.69,Mass_L,Mass_pic);
        CleverMcLevyReso[1].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[1], 2);
        Kitty.SetAnaSource(0,0.87);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        goto CLEAN_SetUpCats_pApCoulomb;
    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(+1);
    Kitty.SetPdgId(3334, 3334);
    Kitty.SetRedMass( (Mass_Omega*Mass_Omega)/(Mass_Omega+Mass_Omega) );

    Kitty.SetNumChannels(4);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);
    Kitty.SetNumPW(2,1);
    Kitty.SetNumPW(3,1);

    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);
    Kitty.SetSpin(2,2);
    Kitty.SetSpin(3,3);

    Kitty.SetChannelWeight(0, 1./16.);
    Kitty.SetChannelWeight(1, 3./16.);
    Kitty.SetChannelWeight(2, 5./16.);
    Kitty.SetChannelWeight(3, 7./16.);

    CLEAN_SetUpCats_pApCoulomb: ;

}
//POT:
//  "AV18"
void DLM_CommonAnaFunctions::SetUpCats_pp(CATS& Kitty, const TString& POT, const TString& SOURCE){

    CATSparameters* cPars = NULL;

    CATSparameters* cPotPars1S0 = NULL;
    CATSparameters* cPotPars3P0 = NULL;
    CATSparameters* cPotPars3P1 = NULL;
    CATSparameters* cPotPars3P2 = NULL;

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.306);//only gaussian
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Cauchy"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Nolan"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,1.2);
        cPars->SetParameter(1,1.6);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Single"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,sqrt(1.6)*1.2);
        cPars->SetParameter(1,1.6);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Diff"){
        cPars->SetParameter(0,0.5*1.6*1.2);
        cPars->SetParameter(1,1.6);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Nolan"){
        CleverLevy[0].InitStability(20,1,2);
        CleverLevy[0].InitScale(35,0.25,2.0);
        CleverLevy[0].InitRad(256,0,64);
        CleverLevy[0].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[0], 2);
        Kitty.SetAnaSource(0,1.2);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Single"){
        CleverLevy[0].InitStability(20,1,2);
        CleverLevy[0].InitScale(35,0.25,2.0);
        CleverLevy[0].InitRad(256,0,64);
        CleverLevy[0].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[0], 2);
        Kitty.SetAnaSource(0,sqrt(1.6)*1.2);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Diff"){
        CleverLevy[0].InitStability(20,1,2);
        CleverLevy[0].InitScale(35,0.25,2.0);
        CleverLevy[0].InitRad(256,0,64);
        CleverLevy[0].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[0], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0,0.5*1.6*1.2);
        Kitty.SetAnaSource(1,1.6);
    }
    else if(SOURCE=="GaussExpTotSimple_2body"){
        //printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Gauss_mT_Reso)\n");
        cPars = new CATSparameters(CATSparameters::tSource,11,true);
        cPars->SetParameter(0,1.2);
        cPars->SetParameter(1,1.65);//tau
        cPars->SetParameter(2,1.-0.3578);//prim
        cPars->SetParameter(3,1361.52);//reso mass
        cPars->SetParameter(4,Mass_p);
        cPars->SetParameter(5,Mass_pic);
        cPars->SetParameter(6,1.65);
        cPars->SetParameter(7,1.-0.3578);
        cPars->SetParameter(8,1361.52);
        cPars->SetParameter(9,Mass_p);
        cPars->SetParameter(10,Mass_pic);
        Kitty.SetAnaSource(GaussExpTotSimple_2body, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="McLevyNolan_Reso"){
        CleverMcLevyReso[0].InitStability(21,1,2);
        CleverMcLevyReso[0].InitScale(38,0.15,2.0);
        CleverMcLevyReso[0].InitRad(257,0,64);
        CleverMcLevyReso[0].InitType(2);
        CleverMcLevyReso[0].InitReso(0,1);
        CleverMcLevyReso[0].InitReso(1,1);
        CleverMcLevyReso[0].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso[0].SetUpReso(1,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso[0].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[0], 2);
        Kitty.SetAnaSource(0,1.2);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="McGauss_Reso"){
        CleverMcLevyReso[0].InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyReso[0].InitScale(38,0.15,2.0);
        CleverMcLevyReso[0].InitRad(1000,0,64);
        CleverMcLevyReso[0].InitType(2);
        CleverMcLevyReso[0].InitReso(0,1);
        CleverMcLevyReso[0].InitReso(1,1);
        CleverMcLevyReso[0].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso[0].SetUpReso(1,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso[0].InitNumMcIter(5000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[0], 2);
        Kitty.SetAnaSource(0,0.96);//put the core
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="EPOS"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pp;
    }
    else if(SOURCE=="EPOSrescaled"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pp;
    }
    else if(SOURCE=="Levy_mT_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pp;
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        goto CLEAN_SetUpCats_pp;
    }

    if(POT=="AV18"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotPars1S0[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
        double PotPars3P0[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
        double PotPars3P1[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
        double PotPars3P2[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,2};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3P0 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P0->SetParameters(PotPars3P0);
        cPotPars3P1 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P1->SetParameters(PotPars3P1);
        cPotPars3P2 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P2->SetParameters(PotPars3P2);
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing pp potential '%s'\n",POT.Data());
        goto CLEAN_SetUpCats_pp;
    }
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(2212, 2212);
    Kitty.SetRedMass( 0.5*Mass_p );

    Kitty.SetNumChannels(4);
    Kitty.SetNumPW(0,2);
    Kitty.SetNumPW(1,2);
    Kitty.SetNumPW(2,2);
    Kitty.SetNumPW(3,2);
    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);
    Kitty.SetSpin(2,1);
    Kitty.SetSpin(3,1);
    Kitty.SetChannelWeight(0, 3./12.);
    Kitty.SetChannelWeight(1, 1./12.);
    Kitty.SetChannelWeight(2, 3./12.);
    Kitty.SetChannelWeight(3, 5./12.);

    if(cPotPars1S0) Kitty.SetShortRangePotential(0,0,fDlmPot,*cPotPars1S0);
    if(cPotPars3P0) Kitty.SetShortRangePotential(1,1,fDlmPot,*cPotPars3P0);
    if(cPotPars3P1) Kitty.SetShortRangePotential(2,1,fDlmPot,*cPotPars3P1);
    if(cPotPars3P2) Kitty.SetShortRangePotential(3,1,fDlmPot,*cPotPars3P2);

    CLEAN_SetUpCats_pp: ;
    if(cPars){delete cPars; cPars=NULL;}
    //if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if(cPotPars1S0){delete cPotPars1S0; cPotPars1S0=NULL;}
    if(cPotPars3P0){delete cPotPars3P0; cPotPars3P0=NULL;}
    if(cPotPars3P1){delete cPotPars3P1; cPotPars3P1=NULL;}
    if(cPotPars3P2){delete cPotPars3P2; cPotPars3P2=NULL;}

}

//Implement Lednicky
void DLM_CommonAnaFunctions::SetUpCats_pAL(CATS& Kitty, const TString& POT, const TString& SOURCE){
std::cout<<"SetUpCats_pAL dummy \n"<<std::endl;
}

//POT:
//  "LO"
//  "LO_Coupled_S"
//  "NLO"
//  "NLO_Coupled_S"
//  "Usmani"
void DLM_CommonAnaFunctions::SetUpCats_pL(CATS& Kitty, const TString& POT, const TString& SOURCE){
    CATSparameters* cPars = NULL;
    CATSparameters* pPars = NULL;

    CATSparameters* cPotPars1S0 = NULL;
    CATSparameters* cPotPars3S1 = NULL;

    DLM_Histo<complex<double>>*** ExternalWF=NULL;
    unsigned NumChannels=0;

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.306);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Cauchy"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Nolan"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,1.2);
        cPars->SetParameter(1,1.2);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Single"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,sqrt(1.2)*1.2);
        cPars->SetParameter(1,1.2);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Diff"){
        cPars->SetParameter(0,0.5*1.2*1.2);
        cPars->SetParameter(1,1.2);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Nolan"){
        CleverLevy[1].InitStability(20,1,2);
        CleverLevy[1].InitScale(35,0.25,2.0);
        CleverLevy[1].InitRad(256,0,64);
        CleverLevy[1].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[1], 2);
        Kitty.SetAnaSource(0,1.2);
        Kitty.SetAnaSource(1,1.2);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Single"){
        CleverLevy[1].InitStability(20,1,2);
        CleverLevy[1].InitScale(35,0.25,2.0);
        CleverLevy[1].InitRad(256,0,64);
        CleverLevy[1].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[1], 2);
        Kitty.SetAnaSource(0,sqrt(1.2)*1.2);
        Kitty.SetAnaSource(1,1.2);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Diff"){
        CleverLevy[1].InitStability(20,1,2);
        CleverLevy[1].InitScale(35,0.25,2.0);
        CleverLevy[1].InitRad(256,0,64);
        CleverLevy[1].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[1], 2);
        Kitty.SetAnaSource(0,0.5*1.2*1.2);
        Kitty.SetAnaSource(1,1.2);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="GaussExpTotSimple_2body"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (GaussExpTotSimple_2body)\n");
        goto CLEAN_SetUpCats_pL;
    }
    else if(SOURCE=="McLevyNolan_Reso"){
        CleverMcLevyReso[1].InitStability(21,1,2);
        CleverMcLevyReso[1].InitScale(38,0.15,2.0);
        CleverMcLevyReso[1].InitRad(257,0,64);
        CleverMcLevyReso[1].InitType(2);
        CleverMcLevyReso[1].InitReso(0,1);
        CleverMcLevyReso[1].InitReso(1,1);
        CleverMcLevyReso[1].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso[1].SetUpReso(1,0,1.-0.3562,1462.93,4.69,Mass_L,Mass_pic);
        CleverMcLevyReso[1].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[1], 2);
        Kitty.SetAnaSource(0,1.2);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="McGauss_Reso"){
        CleverMcLevyReso[1].InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyReso[1].InitScale(35,0.25,2.0);
        CleverMcLevyReso[1].InitRad(256,0,64);
        CleverMcLevyReso[1].InitType(2);
        CleverMcLevyReso[1].InitReso(0,1);
        CleverMcLevyReso[1].InitReso(1,1);
        CleverMcLevyReso[1].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso[1].SetUpReso(1,0,1.-0.3562,1462.93,4.69,Mass_L,Mass_pic);
        CleverMcLevyReso[1].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[1], 2);
        Kitty.SetAnaSource(0,0.87);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="EPOS"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pL;
    }
    else if(SOURCE=="EPOSrescaled"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pL;
    }
    else if(SOURCE=="Levy_mT_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pL;
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        goto CLEAN_SetUpCats_pL;
    }

    if(POT=="LO"){
//        ExternalWF = Init_pL_Haidenbauer(   "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaLO_600/",
//                                Kitty, 0, 600);
//        ExternalWF = Init_pL_Haidenbauer(CatsFilesFolder[0]+"/Interaction/Haidenbauer/pLambdaLO_600/",

        ExternalWF = Init_pL_Haidenbauer(   "/Users/sartozza/FemtoStuff/Haidenbauer/pLambdaLO_600/",
                                                        Kitty, 0, 600);
        NumChannels=2;
    }
    else if(POT=="LO_Coupled_S"){
        ExternalWF = Init_pL_Haidenbauer(   "/Users/sartozza/FemtoStuff/Haidenbauer/pLambdaLO_Coupling/",
                                Kitty, 1, 600);
        NumChannels=4;
    }
    else if(POT=="NLO"){
        ExternalWF = Init_pL_Haidenbauer(   "/Users/sartozza/FemtoStuff/Haidenbauer/pLambdaNLO/",
                                Kitty, 10, 600);
        NumChannels=2;
    }
    //s and p waves
    else if(POT=="NLO_sp"){
        ExternalWF = Init_pL_Haidenbauer(   "/Users/sartozza/FemtoStuff/Haidenbauer/pLambdaNLO/",
                                Kitty, 12, 600);
        NumChannels=4;
    }
    else if(POT=="NLO_Coupled_S"){
        ExternalWF = Init_pL_Haidenbauer(   "/Users/sartozza/FemtoStuff/Haidenbauer/pLambdaNLO_Coupling/",
                                Kitty, 11, 600);
        NumChannels=4;
    }
    else if(POT=="Usmani"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotPars1S0[8]={pL_UsmaniOli,0,0,0,0,0,0,0};
        double PotPars3S1[8]={pL_UsmaniOli,0,0,0,0,1,0,1};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3S1 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3S1->SetParameters(PotPars3S1);
        NumChannels=2;
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing pp potential '%s'\n",POT.Data());
        goto CLEAN_SetUpCats_pL;
    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3122);
    Kitty.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );

    Kitty.SetNumChannels(NumChannels);
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        if(!ExternalWF){
            Kitty.SetNumPW(uCh,1);
            Kitty.SetSpin(uCh, uCh%2==0?0:1);
            Kitty.SetChannelWeight(uCh, uCh%2==0?0.25:0.75);
        }


        if(cPotPars1S0&&uCh==0)Kitty.SetShortRangePotential(uCh,0,fDlmPot,*cPotPars1S0);
        else if(cPotPars3S1&&uCh==1) Kitty.SetShortRangePotential(uCh,0,fDlmPot,*cPotPars3S1);
        else if(ExternalWF){
            //for(unsigned uMomBin=0; uMomBin<Kitty.GetNumMomBins(); uMomBin++){
                //Kitty.UseExternalWaveFunction(uMomBin,uCh,0,WaveFunctionU[uMomBin][uCh][0], NumRadBins, RadBins, PhaseShifts[uMomBin][uCh][0]);
//printf("Look at that view (%u)!\n",uCh);
                Kitty.SetExternalWaveFunction(uCh,0,ExternalWF[0][uCh][0],ExternalWF[1][uCh][0]);
                if(Kitty.GetNumPW(uCh)==2){
                    Kitty.SetExternalWaveFunction(uCh,1,ExternalWF[0][uCh][1],ExternalWF[1][uCh][1]);
                }
//printf(" --Look at that view (%u)!\n",uCh);
            //}
        }
        else{
            printf("\033[1;31mERROR:\033[0m SetUpCats_pL says that you should NEVER see this message! BIG BUG!\n");
            goto CLEAN_SetUpCats_pL;
        }

    }
//Kitty.KillTheCat();
//printf("------------------------");
    CLEAN_SetUpCats_pL: ;

    if(cPars){delete cPars; cPars=NULL;}
    if(pPars){delete pPars; pPars=NULL;}
    //if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if(cPotPars1S0){delete cPotPars1S0; cPotPars1S0=NULL;}
    if(cPotPars3S1){delete cPotPars3S1; cPotPars3S1=NULL;}
    CleanUpWfHisto(Kitty,ExternalWF);

}

//Implement Lednicky
void DLM_CommonAnaFunctions::SetUpCats_LAL(CATS& Kitty, const TString& POT, const TString& SOURCE){
std::cout<<"SetUpCats_LAL dummy \n"<<std::endl;
}



//POT:
//  "pXim_Lattice" (the first version)
//  "pXim_HALQCD1" (the second version)
//  "pXim_ESC16_IS"

void DLM_CommonAnaFunctions::SetUpCats_pXim(CATS& Kitty, const TString& POT, const TString& SOURCE){
    CATSparameters* cPars = NULL;
    CATSparameters* pPars = NULL;

    CATSparameters* cPotParsI0S0 = NULL;
    CATSparameters* cPotParsI0S1 = NULL;
    CATSparameters* cPotParsI1S0 = NULL;
    CATSparameters* cPotParsI1S1 = NULL;

    DLM_Histo<complex<double>>*** ExternalWF=NULL;
    unsigned NumChannels=0;

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.306);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Cauchy"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Nolan"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,1.2);
        cPars->SetParameter(1,1.8);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Single"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,sqrt(1.8)*1.2);
        cPars->SetParameter(1,1.8);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Diff"){
        cPars->SetParameter(0,0.5*1.8*1.2);
        cPars->SetParameter(1,1.8);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Nolan"){
        CleverLevy[2].InitStability(20,1,2);
        CleverLevy[2].InitScale(35,0.25,2.0);
        CleverLevy[2].InitRad(256,0,64);
        CleverLevy[2].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[2], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0,1.2);
        Kitty.SetAnaSource(1,1.8);
    }
    else if(SOURCE=="CleverLevy_Single"){
        CleverLevy[2].InitStability(20,1,2);
        CleverLevy[2].InitScale(35,0.25,2.0);
        CleverLevy[2].InitRad(256,0,64);
        CleverLevy[2].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[2], 2);
        Kitty.SetAnaSource(0,sqrt(1.8)*1.2);
        Kitty.SetAnaSource(1,1.8);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Diff"){
        CleverLevy[2].InitStability(20,1,2);
        CleverLevy[2].InitScale(35,0.25,2.0);
        CleverLevy[2].InitRad(256,0,64);
        CleverLevy[2].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[2], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0,0.5*1.8*1.2);
        Kitty.SetAnaSource(1,1.8);
    }
    else if(SOURCE=="GaussExpTotSimple_2body"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (GaussExpTotSimple_2body)\n");
        goto CLEAN_SetUpCats_pp;
    }
    else if(SOURCE=="McLevyNolan_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (McLevyNolan_Reso)\n");
        goto CLEAN_SetUpCats_pp;
    }
    else if(SOURCE=="EPOS"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pp;
    }
    else if(SOURCE=="EPOSrescaled"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pp;
    }
    else if(SOURCE=="Levy_mT_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pp;
    }
    else if(SOURCE=="McGauss_Reso"){
        CleverMcLevyReso[2].InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyReso[2].InitScale(35,0.25,2.0);
        CleverMcLevyReso[2].InitRad(256,0,64);
        CleverMcLevyReso[2].InitType(2);
        CleverMcLevyReso[2].InitReso(0,1);
        CleverMcLevyReso[2].InitReso(1,0);
        CleverMcLevyReso[2].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso[2].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[2], 2);
        Kitty.SetAnaSource(0,0.77);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
      }
      else{
          printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
          goto CLEAN_SetUpCats_pp;
      }

    if(POT=="pXim_Lattice"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI0S0[9]={pXim_Lattice,12,0,-1,1,0,0,0,0};
        double PotParsI0S1[9]={pXim_Lattice,12,0,-1,1,1,0,1,0};
        double PotParsI1S0[9]={pXim_Lattice,6,1,1,1,0,0,0,0};
        double PotParsI1S1[9]={pXim_Lattice,6,1,1,1,1,0,1,0};
        cPotParsI0S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S0->SetParameters(PotParsI0S0);
        cPotParsI0S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S1->SetParameters(PotParsI0S1);
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S1->SetParameters(PotParsI1S1);
        NumChannels=4;
    }
    else if(POT=="pXim_HALQCD1"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI0S0[9]={pXim_HALQCD1,12,0,-1,1,0,0,0,0};
        double PotParsI0S1[9]={pXim_HALQCD1,12,0,-1,1,1,0,1,0};
        double PotParsI1S0[9]={pXim_HALQCD1,12,1,1,1,0,0,0,0};
        double PotParsI1S1[9]={pXim_HALQCD1,12,1,1,1,1,0,1,0};
        cPotParsI0S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S0->SetParameters(PotParsI0S0);
        cPotParsI0S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S1->SetParameters(PotParsI0S1);
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S1->SetParameters(PotParsI1S1);
        NumChannels=4;
    }
//     else if(POT=="pXim_ESC16_IS"){
// //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
// // here POT_FLAG = 12 is dummy
// double PotParsI0S0[9]={pXim_ESC16_IS,12,0,-1,1,0,0,0,0};
// double PotParsI0S1[9]={pXim_ESC16_IS,12,0,-1,1,1,0,1,0};
// double PotParsI1S0[9]={pXim_ESC16_IS,12,1,1,1,0,0,0,0};
// double PotParsI1S1[9]={pXim_ESC16_IS,12,1,1,1,1,0,1,0};
// cPotParsI0S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S0->SetParameters(PotParsI0S0);
// cPotParsI0S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S1->SetParameters(PotParsI0S1);
// cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S0->SetParameters(PotParsI1S0);
// cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S1->SetParameters(PotParsI1S1);
// NumChannels=4;
//     }

    else{
        printf("\033[1;31mERROR:\033[0m Non-existing pXim potential '%s'\n",POT.Data());
        goto CLEAN_SetUpCats_pp;
//        goto CLEAN_SetUpCats_pL;

    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, 3122);
    Kitty.SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );
    Kitty.SetNumChannels(NumChannels);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);
    Kitty.SetNumPW(2,1);
    Kitty.SetNumPW(3,1);
    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);
    Kitty.SetSpin(2,0);
    Kitty.SetSpin(3,1);
    Kitty.SetChannelWeight(0, 1./8.);
    Kitty.SetChannelWeight(1, 3./8.);
    Kitty.SetChannelWeight(2, 1./8.);
    Kitty.SetChannelWeight(3, 3./8.);

    if(cPotParsI0S0) Kitty.SetShortRangePotential(0,0,fDlmPot,*cPotParsI0S0);
    if(cPotParsI0S1) Kitty.SetShortRangePotential(1,0,fDlmPot,*cPotParsI0S1);
    if(cPotParsI1S0) Kitty.SetShortRangePotential(2,0,fDlmPot,*cPotParsI1S0);
    if(cPotParsI1S1) Kitty.SetShortRangePotential(3,0,fDlmPot,*cPotParsI1S1);


    CLEAN_SetUpCats_pp: ;
//    CLEAN_SetUpCats_pL: ;
    if(cPars){delete cPars; cPars=NULL;}
    if(pPars){delete pPars; pPars=NULL;}
    //if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if(cPotParsI0S0){delete cPotParsI0S0; cPotParsI0S0=NULL;}
    if(cPotParsI0S1){delete cPotParsI0S1; cPotParsI0S1=NULL;}
    if(cPotParsI1S0){delete cPotParsI1S0; cPotParsI1S0=NULL;}
    if(cPotParsI1S1){delete cPotParsI1S1; cPotParsI1S1=NULL;}

}

void DLM_CommonAnaFunctions::SetUpCats_pApHaide(CATS& Kitty, const TString& POT, const TString& SOURCE, const TString& DataSample){
    CATSparameters* cPars = NULL;
    CATSparameters* pPars = NULL;

    CATSparameters* cPotPars1S0 = NULL;
    CATSparameters* cPotPars3S1 = NULL;


    DLM_Histo<complex<double>>*** ExternalWF=NULL;
    unsigned NumChannels=0;

    double radius;
    if(DataSample=="pp13TeV_MB_BBar") radius=1.188;
    if(DataSample=="pp13TeV_HM_BBar") radius=1.25;

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,radius);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="McGauss_Reso"){
        CleverMcLevyReso[1].InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyReso[1].InitScale(35,0.25,2.0);
        CleverMcLevyReso[1].InitRad(256,0,64);
        CleverMcLevyReso[1].InitType(2);
        CleverMcLevyReso[1].InitReso(0,1);
        CleverMcLevyReso[1].InitReso(1,1);
        CleverMcLevyReso[1].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso[1].SetUpReso(1,0,1.-0.3562,1462.93,4.69,Mass_L,Mass_pic);
        CleverMcLevyReso[1].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[1], 2);
        Kitty.SetAnaSource(0,0.87);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        goto CLEAN_SetUpCats_pApHaide;
    }

    if(POT=="HAIDE_1"){//only ppbar->ppbar
//        ExternalWF = Init_pL_Haidenbauer(   "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaLO_600/",
//                                Kitty, 0, 600);
        ExternalWF=Init_pantip_Haidenbauer("/Users/sartozza/cernbox/Analysis/BBbar/Wavefunctions/Haidenbauer/p_antip_wCoulomb/wf_18092019/",Kitty,0);

        NumChannels=4;
    }
    else if(POT=="HAIDE_2"){//ppbar->ppbar + nnbar->ppbar
      ExternalWF=Init_pantip_Haidenbauer("/Users/sartozza/cernbox/Analysis/BBbar/Wavefunctions/Haidenbauer/p_antip_wCoulomb/wf_18092019/",Kitty,1);
      NumChannels=8;
    } else if(POT=="HAIDE_3"){//ppbar->ppbar + nnbar->ppbar + "2pi"->ppbar in 1S0
      ExternalWF=Init_pantip_Haidenbauer("/Users/sartozza/cernbox/Analysis/BBbar/Wavefunctions/Haidenbauer/p_antip_wCoulomb/wf_18092019/",Kitty,2);
      NumChannels=9;
    } else if(POT=="HAIDE_4"){//ppbar->ppbar + nnbar->ppbar + "2pi"->ppbar in 3S1
      ExternalWF=Init_pantip_Haidenbauer("/Users/sartozza/cernbox/Analysis/BBbar/Wavefunctions/Haidenbauer/p_antip_wCoulomb/wf_18092019/",Kitty,3);
      NumChannels=9;
    } else if(POT=="HAIDE_5"){//ppbar->ppbar + nnbar->ppbar + "2pi"->ppbar in all PWs
      ExternalWF=Init_pantip_Haidenbauer("/Users/sartozza/cernbox/Analysis/BBbar/Wavefunctions/Haidenbauer/p_antip_wCoulomb/wf_18092019/",Kitty,4);
      NumChannels=12;
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing pp potential '%s'\n",POT.Data());
        goto CLEAN_SetUpCats_pApHaide;
    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, -2212);
    Kitty.SetRedMass( (Mass_p*Mass_p)/(Mass_p+Mass_p) );

    Kitty.SetNumChannels(NumChannels);
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
      if(ExternalWF){
      //Setting the wfs for the different channels
      if(POT=="HAIDE_1"){//only ppbar->ppbar           ch pw               ch pw
        Kitty.SetExternalWaveFunction(0,0,ExternalWF[0][0][0],ExternalWF[1][0][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(0,1,ExternalWF[0][0][1],ExternalWF[1][0][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(1,0,ExternalWF[0][1][0],ExternalWF[1][1][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(1,1,ExternalWF[0][1][1],ExternalWF[1][1][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(2,0,ExternalWF[0][2][0],ExternalWF[1][2][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(2,1,ExternalWF[0][2][1],ExternalWF[1][2][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(3,0,ExternalWF[0][3][0],ExternalWF[1][3][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(3,1,ExternalWF[0][3][1],ExternalWF[1][3][1]);//ExternalWF_pantip[1] is phase shifts
      } else if(POT=="HAIDE_2"){//ppbar->ppbar + nnbar->ppbar
        Kitty.SetExternalWaveFunction(0,0,ExternalWF[0][0][0],ExternalWF[1][0][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(0,1,ExternalWF[0][0][1],ExternalWF[1][0][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(1,0,ExternalWF[0][1][0],ExternalWF[1][1][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(1,1,ExternalWF[0][1][1],ExternalWF[1][1][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(2,0,ExternalWF[0][2][0],ExternalWF[1][2][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(2,1,ExternalWF[0][2][1],ExternalWF[1][2][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(3,0,ExternalWF[0][3][0],ExternalWF[1][3][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(3,1,ExternalWF[0][3][1],ExternalWF[1][3][1]);//ExternalWF_pantip[1] is phase shifts

// nAn part
        Kitty.SetExternalWaveFunction(4,0,ExternalWF[0][4][0],ExternalWF[1][4][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(4,1,ExternalWF[0][4][1],ExternalWF[1][4][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(5,0,ExternalWF[0][5][0],ExternalWF[1][5][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(5,1,ExternalWF[0][5][1],ExternalWF[1][5][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(6,0,ExternalWF[0][6][0],ExternalWF[1][6][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(6,1,ExternalWF[0][6][1],ExternalWF[1][6][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(7,0,ExternalWF[0][7][0],ExternalWF[1][7][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(7,1,ExternalWF[0][7][1],ExternalWF[1][7][1]);//ExternalWF_pantip[1] is phase shifts
      }
else if(POT=="HAIDE_3" || POT=="HAIDE_4"){//ppbar->ppbar + nnbar->ppbar + "2pi"->ppbar
        Kitty.SetExternalWaveFunction(0,0,ExternalWF[0][0][0],ExternalWF[1][0][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(0,1,ExternalWF[0][0][1],ExternalWF[1][0][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(1,0,ExternalWF[0][1][0],ExternalWF[1][1][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(1,1,ExternalWF[0][1][1],ExternalWF[1][1][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(2,0,ExternalWF[0][2][0],ExternalWF[1][2][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(2,1,ExternalWF[0][2][1],ExternalWF[1][2][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(3,0,ExternalWF[0][3][0],ExternalWF[1][3][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(3,1,ExternalWF[0][3][1],ExternalWF[1][3][1]);//ExternalWF_pantip[1] is phase shifts

// nAn part
        Kitty.SetExternalWaveFunction(4,0,ExternalWF[0][4][0],ExternalWF[1][4][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(4,1,ExternalWF[0][4][1],ExternalWF[1][4][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(5,0,ExternalWF[0][5][0],ExternalWF[1][5][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(5,1,ExternalWF[0][5][1],ExternalWF[1][5][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(6,0,ExternalWF[0][6][0],ExternalWF[1][6][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(6,1,ExternalWF[0][6][1],ExternalWF[1][6][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(7,0,ExternalWF[0][7][0],ExternalWF[1][7][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(7,1,ExternalWF[0][7][1],ExternalWF[1][7][1]);//ExternalWF_pantip[1] is phase shifts
// "2pi" part
        Kitty.SetExternalWaveFunction(8,0,ExternalWF[0][8][0],ExternalWF[1][8][0]);//ExternalWF_pantip[1] is phase shifts
      }
else if(POT=="HAIDE_5"){//ppbar->ppbar + nnbar->ppbar + "2pi"->ppbar in all PW
        Kitty.SetExternalWaveFunction(0,0,ExternalWF[0][0][0],ExternalWF[1][0][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(0,1,ExternalWF[0][0][1],ExternalWF[1][0][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(1,0,ExternalWF[0][1][0],ExternalWF[1][1][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(1,1,ExternalWF[0][1][1],ExternalWF[1][1][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(2,0,ExternalWF[0][2][0],ExternalWF[1][2][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(2,1,ExternalWF[0][2][1],ExternalWF[1][2][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(3,0,ExternalWF[0][3][0],ExternalWF[1][3][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(3,1,ExternalWF[0][3][1],ExternalWF[1][3][1]);//ExternalWF_pantip[1] is phase shifts

// nAn part
        Kitty.SetExternalWaveFunction(4,0,ExternalWF[0][4][0],ExternalWF[1][4][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(4,1,ExternalWF[0][4][1],ExternalWF[1][4][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(5,0,ExternalWF[0][5][0],ExternalWF[1][5][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(5,1,ExternalWF[0][5][1],ExternalWF[1][5][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(6,0,ExternalWF[0][6][0],ExternalWF[1][6][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(6,1,ExternalWF[0][6][1],ExternalWF[1][6][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(7,0,ExternalWF[0][7][0],ExternalWF[1][7][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(7,1,ExternalWF[0][7][1],ExternalWF[1][7][1]);//ExternalWF_pantip[1] is phase shifts
// "2pi" part
        Kitty.SetExternalWaveFunction(8,0,ExternalWF[0][8][0],ExternalWF[1][8][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(8,1,ExternalWF[0][8][1],ExternalWF[1][8][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(9,0,ExternalWF[0][9][0],ExternalWF[1][9][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(9,1,ExternalWF[0][9][1],ExternalWF[1][9][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(10,0,ExternalWF[0][10][0],ExternalWF[1][10][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(10,1,ExternalWF[0][10][1],ExternalWF[1][10][1]);//ExternalWF_pantip[1] is phase shifts

        Kitty.SetExternalWaveFunction(11,0,ExternalWF[0][11][0],ExternalWF[1][11][0]);//ExternalWF_pantip[1] is phase shifts
        Kitty.SetExternalWaveFunction(11,1,ExternalWF[0][11][1],ExternalWF[1][11][1]);//ExternalWF_pantip[1] is phase shifts
      }


    }
      else{
            printf("\033[1;31mERROR:\033[0m SetUpCats_pApHaide says that you should NEVER see this message! BIG BUG!\n");
            goto CLEAN_SetUpCats_pApHaide;
        }

    }//end of for
//Kitty.KillTheCat();

    CLEAN_SetUpCats_pApHaide: ;

    if(cPars){delete cPars; cPars=NULL;}
    if(pPars){delete pPars; pPars=NULL;}
    //if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if(cPotPars1S0){delete cPotPars1S0; cPotPars1S0=NULL;}
    if(cPotPars3S1){delete cPotPars3S1; cPotPars3S1=NULL;}
    CleanUpWfHisto(Kitty,ExternalWF);
}

//----------------------SETTING UP BINNING---------------------------
void DLM_CommonAnaFunctions::SetUpBinning_pAp(const TString& DataSample, unsigned& NumMomBins, double*& MomBins, double*& FitRegion){
  if(DataSample=="pp13TeV_MB_Run2paper"){
      const double kMin=0;
      const double kStep=4;
      NumMomBins=94;//(i.e. max=376 MeV)
      if(MomBins) delete [] MomBins;
      MomBins = new double [NumMomBins+1];
      MomBins[0] = kMin;
      for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
          MomBins[uBin] = MomBins[uBin-1]+kStep;
      }
      if(FitRegion) delete [] FitRegion;
      FitRegion = new double [4];
      FitRegion[0] = MomBins[0];
      FitRegion[1] = MomBins[NumMomBins];
      FitRegion[2] = MomBins[NumMomBins]+kStep;
      FitRegion[3] = MomBins[NumMomBins]+kStep*31.;//till 500
  }
  else if (DataSample=="pp13TeV_MB_BBar" || DataSample=="pp13TeV_HM_BBar") {
        const double kMin=0;
        const double kStep=4;
        // NumMomBins=75;//94(376)//500->125 300->75//400->100
        if(MomBins) delete [] MomBins;
        MomBins = new double [NumMomBins+1];
        MomBins[0] = kMin;
        for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
            MomBins[uBin] = MomBins[uBin-1]+kStep;
        }
        if(FitRegion) delete [] FitRegion;
        FitRegion = new double [4];
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins]+kStep;
        FitRegion[3] = MomBins[NumMomBins]+kStep*31.;//till 624
    }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        NumMomBins=0;
        return;
    }
}

void DLM_CommonAnaFunctions::SetUpBinning_pp(const TString& DataSample, unsigned& NumMomBins, double*& MomBins, double*& FitRegion){
    if(DataSample=="pp13TeV_MB_Run2paper"){
        const double kMin=0;
        const double kStep=4;
        NumMomBins=94;//(i.e. max=376 MeV)
        if(MomBins) delete [] MomBins;
        MomBins = new double [NumMomBins+1];
        MomBins[0] = kMin;
        for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
            MomBins[uBin] = MomBins[uBin-1]+kStep;
        }
        if(FitRegion) delete [] FitRegion;
        FitRegion = new double [4];
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins]+kStep;
        FitRegion[3] = MomBins[NumMomBins]+kStep*31.;//till 500
    }
    else if(DataSample=="pp13TeV_HM_March19"){
        const double kMin=4;
        const double kStep=4;
        NumMomBins=94;//(i.e. max=376 MeV)
        if(MomBins) delete [] MomBins;
        MomBins = new double [NumMomBins+1];
        MomBins[0] = kMin;
        for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
            MomBins[uBin] = MomBins[uBin-1]+kStep;
        }
        if(FitRegion) delete [] FitRegion;
        FitRegion = new double [4];
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins]+kStep;
        FitRegion[3] = MomBins[NumMomBins]+kStep*31.;//till 500
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        const double kMin=4;
        const double kStep=4;
        NumMomBins=93;//(i.e. max=376 MeV)
        if(MomBins) delete [] MomBins;
        MomBins = new double [NumMomBins+1];
        MomBins[0] = kMin;
        for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
            MomBins[uBin] = MomBins[uBin-1]+kStep;
        }
        if(FitRegion) delete [] FitRegion;
        FitRegion = new double [4];
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins]+kStep;
        FitRegion[3] = MomBins[NumMomBins]+kStep*31.;//till 500
    }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        NumMomBins=0;
        return;
    }
}
//bin width 16 MeV (_0)
void DLM_CommonAnaFunctions::SetUpBinning_pAL(const TString& DataSample, unsigned& NumMomBins, double*& MomBins, double*& FitRegion, double& kfitrange){
  if(DataSample=="pp13TeV_MB_Run2paper"){
      const double kMin=0;
      const double kFineMin=336;//272//216
      const double kFineMax=336;//304
      const double kMax=336;//upper value of last train
      const double kCoarseStep=12;
      const double kFineStep=12;

      //the number of coarse bins below kFineMin
      //floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
      //and if there is rounding needed it is done so that we make our region larger, not smaller!
      unsigned NumCoarseBinsBelow = floor((kFineMin-kMin)/kCoarseStep);
      unsigned NumFineBins = ceil((kFineMax-double(NumCoarseBinsBelow)*kCoarseStep)/kFineStep);
      //we floor the highest point, to make sure we do not run out of the range provided by experimental data
      unsigned NumCoarseBinsAbove = floor((kMax-double(NumCoarseBinsBelow)*kCoarseStep-double(NumFineBins)*kFineStep)/kCoarseStep);

      NumMomBins=NumCoarseBinsBelow+NumFineBins+NumCoarseBinsAbove;

      if(MomBins) delete [] MomBins;
      MomBins = new double [NumMomBins+1];
      MomBins[0] = kMin;
      for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
          if(uBin<=NumCoarseBinsBelow||uBin>NumCoarseBinsBelow+NumFineBins){
              MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;
          }
          else{
              MomBins[uBin] = MomBins[uBin-1]+kFineStep;
          }
      }
      if(FitRegion) delete [] FitRegion;
      FitRegion = new double [4];
      FitRegion[0] = MomBins[0];
      FitRegion[1] = MomBins[NumMomBins];
      FitRegion[2] = MomBins[NumMomBins]+kCoarseStep;
      FitRegion[3] = MomBins[NumMomBins]+kCoarseStep*10.;//till 496
  }
  else if(DataSample=="pp13TeV_MB_BBar" || DataSample=="pp13TeV_HM_BBar"){
      const double kMin=0;
      // const double kFineMin=304;//608//544//448//400//352//336//272//216//496
      // const double kFineMax=304;//304//496
      // const double kMax=304;//336
      double kFineMin=kfitrange;//608//544//448//400//352//336//272//216//496
      double kFineMax=kfitrange;//304//496
      double kMax=kfitrange;//336
      const double kCoarseStep=16;
      const double kFineStep=16;

      //the number of coarse bins below kFineMin
      //floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
      //and if there is rounding needed it is done so that we make our region larger, not smaller!
      unsigned NumCoarseBinsBelow = floor((kFineMin-kMin)/kCoarseStep);
      unsigned NumFineBins = ceil((kFineMax-double(NumCoarseBinsBelow)*kCoarseStep)/kFineStep);
      //we floor the highest point, to make sure we do not run out of the range provided by experimental data
      unsigned NumCoarseBinsAbove = floor((kMax-double(NumCoarseBinsBelow)*kCoarseStep-double(NumFineBins)*kFineStep)/kCoarseStep);

      NumMomBins=NumCoarseBinsBelow+NumFineBins+NumCoarseBinsAbove;

      if(MomBins) delete [] MomBins;
      MomBins = new double [NumMomBins+1];
      MomBins[0] = kMin;
      for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
          if(uBin<=NumCoarseBinsBelow||uBin>NumCoarseBinsBelow+NumFineBins){
              MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;
          }
          else{
              MomBins[uBin] = MomBins[uBin-1]+kFineStep;
          }
      }
      if(FitRegion) delete [] FitRegion;
      FitRegion = new double [4];
      FitRegion[0] = MomBins[0];
      FitRegion[1] = MomBins[NumMomBins];
      FitRegion[2] = MomBins[NumMomBins]+kCoarseStep;
      FitRegion[3] = MomBins[NumMomBins]+kCoarseStep*10.;//till 496
  }

}

void DLM_CommonAnaFunctions::SetUpBinning_pL(const TString& DataSample, unsigned& NumMomBins, double*& MomBins, double*& FitRegion){

    if(DataSample=="pp13TeV_MB_Run2paper"){
        const double kMin=0;
        const double kFineMin=336;//272//216
        const double kFineMax=336;//304
        const double kMax=336;//336
        const double kCoarseStep=12;
        const double kFineStep=12;

        //the number of coarse bins below kFineMin
        //floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
        //and if there is rounding needed it is done so that we make our region larger, not smaller!
        unsigned NumCoarseBinsBelow = floor((kFineMin-kMin)/kCoarseStep);
        unsigned NumFineBins = ceil((kFineMax-double(NumCoarseBinsBelow)*kCoarseStep)/kFineStep);
        //we floor the highest point, to make sure we do not run out of the range provided by experimental data
        unsigned NumCoarseBinsAbove = floor((kMax-double(NumCoarseBinsBelow)*kCoarseStep-double(NumFineBins)*kFineStep)/kCoarseStep);

        NumMomBins=NumCoarseBinsBelow+NumFineBins+NumCoarseBinsAbove;

        if(MomBins) delete [] MomBins;
        MomBins = new double [NumMomBins+1];
        MomBins[0] = kMin;
        for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
            if(uBin<=NumCoarseBinsBelow||uBin>NumCoarseBinsBelow+NumFineBins){
                MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;
            }
            else{
                MomBins[uBin] = MomBins[uBin-1]+kFineStep;
            }
        }
        if(FitRegion) delete [] FitRegion;
        FitRegion = new double [4];
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins]+kCoarseStep;
        FitRegion[3] = MomBins[NumMomBins]+kCoarseStep*10.;//till 496
    }
    else if(DataSample=="pp13TeV_HM_March19"){
        const double kMin=0;
        const double kFineMin=336;
        const double kFineMax=336;
        const double kMax=336;
        const double kCoarseStep=12;
        const double kFineStep=12;

        //the number of coarse bins below kFineMin
        //floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
        //and if there is rounding needed it is done so that we make our region larger, not smaller!
        unsigned NumCoarseBinsBelow = floor((kFineMin-kMin)/kCoarseStep);
        unsigned NumFineBins = ceil((kFineMax-double(NumCoarseBinsBelow)*kCoarseStep)/kFineStep);
        //we floor the highest point, to make sure we do not run out of the range provided by experimental data
        unsigned NumCoarseBinsAbove = floor((kMax-double(NumCoarseBinsBelow)*kCoarseStep-double(NumFineBins)*kFineStep)/kCoarseStep);

        NumMomBins=NumCoarseBinsBelow+NumFineBins+NumCoarseBinsAbove;

        if(MomBins) delete [] MomBins;
        MomBins = new double [NumMomBins+1];
        MomBins[0] = kMin;
        for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
            if(uBin<=NumCoarseBinsBelow||uBin>NumCoarseBinsBelow+NumFineBins){
                MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;
            }
            else{
                MomBins[uBin] = MomBins[uBin-1]+kFineStep;
            }
        }
        if(FitRegion) delete [] FitRegion;
        FitRegion = new double [4];
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins]+kCoarseStep;
        FitRegion[3] = MomBins[NumMomBins]+kCoarseStep*10.;//till 496
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        const double kMin=0;
        const double kFineMin=336;
        const double kFineMax=336;
        const double kMax=336;
        const double kCoarseStep=12;
        const double kFineStep=12;

        //the number of coarse bins below kFineMin
        //floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
        //and if there is rounding needed it is done so that we make our region larger, not smaller!
        unsigned NumCoarseBinsBelow = floor((kFineMin-kMin)/kCoarseStep);
        unsigned NumFineBins = ceil((kFineMax-double(NumCoarseBinsBelow)*kCoarseStep)/kFineStep);
        //we floor the highest point, to make sure we do not run out of the range provided by experimental data
        unsigned NumCoarseBinsAbove = floor((kMax-double(NumCoarseBinsBelow)*kCoarseStep-double(NumFineBins)*kFineStep)/kCoarseStep);

        NumMomBins=NumCoarseBinsBelow+NumFineBins+NumCoarseBinsAbove;

        if(MomBins) delete [] MomBins;
        MomBins = new double [NumMomBins+1];
        MomBins[0] = kMin;
        for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
            if(uBin<=NumCoarseBinsBelow||uBin>NumCoarseBinsBelow+NumFineBins){
                MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;
            }
            else{
                MomBins[uBin] = MomBins[uBin-1]+kFineStep;
            }
        }
        if(FitRegion) delete [] FitRegion;
        FitRegion = new double [4];
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins]+kCoarseStep;
        FitRegion[3] = MomBins[NumMomBins]+kCoarseStep*10.;//till 496
    }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        NumMomBins=0;
        return;
    }

}
//bin width 16 MeV
void DLM_CommonAnaFunctions::SetUpBinning_LAL(const TString& DataSample, unsigned& NumMomBins, double*& MomBins, double*& FitRegion, double& kfitrange){
  if(DataSample=="pp13TeV_MB_Run2paper"){
      const double kMin=4;
      const double kFineMin=304;//272//216
      const double kFineMax=304;//304//500
      const double kMax=304;//336
      const double kCoarseStep=16;
      const double kFineStep=16;
      //the number of coarse bins below kFineMin
      //floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
      //and if there is rounding needed it is done so that we make our region larger, not smaller!
      unsigned NumCoarseBinsBelow = floor((kFineMin-kMin)/kCoarseStep);
      unsigned NumFineBins = ceil((kFineMax-double(NumCoarseBinsBelow)*kCoarseStep)/kFineStep);
      //we floor the highest point, to make sure we do not run out of the range provided by experimental data
      unsigned NumCoarseBinsAbove = floor((kMax-double(NumCoarseBinsBelow)*kCoarseStep-double(NumFineBins)*kFineStep)/kCoarseStep);

      NumMomBins=NumCoarseBinsBelow+NumFineBins+NumCoarseBinsAbove;

      if(MomBins) delete [] MomBins;
      MomBins = new double [NumMomBins+1];
      MomBins[0] = kMin;
      for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
          if(uBin<=NumCoarseBinsBelow||uBin>NumCoarseBinsBelow+NumFineBins){
              MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;
          }
          else{
              MomBins[uBin] = MomBins[uBin-1]+kFineStep;
          }
      }
      if(FitRegion) delete [] FitRegion;
      FitRegion = new double [4];
      FitRegion[0] = MomBins[0];
      FitRegion[1] = MomBins[NumMomBins];
      FitRegion[2] = MomBins[NumMomBins]+kCoarseStep;
      FitRegion[3] = MomBins[NumMomBins]+kCoarseStep*10.;//till 496
  } else if(DataSample=="pp13TeV_MB_BBar" || DataSample=="pp13TeV_HM_BBar"){
      const double kMin=4;
      // const double kFineMin=308;//548//404//356//500//272//222
      // const double kFineMax=308;//308
      // const double kMax=308;//336//348
      double kFineMin=kfitrange;//548//404//356//500//272//222
      double kFineMax=kfitrange;//308
      double kMax=kfitrange;//336//348
      const double kCoarseStep=16;
      const double kFineStep=16;

      //the number of coarse bins below kFineMin
      //floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
      //and if there is rounding needed it is done so that we make our region larger, not smaller!
      unsigned NumCoarseBinsBelow = floor((kFineMin-kMin)/kCoarseStep);
      unsigned NumFineBins = ceil((kFineMax-double(NumCoarseBinsBelow)*kCoarseStep)/kFineStep);
      //we floor the highest point, to make sure we do not run out of the range provided by experimental data
      unsigned NumCoarseBinsAbove = floor((kMax-double(NumCoarseBinsBelow)*kCoarseStep-double(NumFineBins)*kFineStep)/kCoarseStep);

      NumMomBins=NumCoarseBinsBelow+NumFineBins+NumCoarseBinsAbove;

      if(MomBins) delete [] MomBins;
      MomBins = new double [NumMomBins+1];
      MomBins[0] = kMin;
      for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
          if(uBin<=NumCoarseBinsBelow||uBin>NumCoarseBinsBelow+NumFineBins){
              MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;
          }
          else{
              MomBins[uBin] = MomBins[uBin-1]+kFineStep;
          }
      }
      if(FitRegion) delete [] FitRegion;
      FitRegion = new double [4];
      FitRegion[0] = MomBins[0];
      FitRegion[1] = MomBins[NumMomBins];
      FitRegion[2] = MomBins[NumMomBins]+kCoarseStep;
      FitRegion[3] = MomBins[NumMomBins]+kCoarseStep*10.;//till 496
  }
}

//------------------------------------

void DLM_CommonAnaFunctions::GetPurities_Ap(const TString& DataSample, const int& Variation, double* Purities){
  double PurityAntiProton;
  if(DataSample=="pp13TeV_MB_Run2paper"){
      PurityAntiProton = 0.989859;
  }else if(DataSample=="pp13TeV_MB_BBar"){
      PurityAntiProton = 0.989859;
  }
  else if(DataSample=="pp13TeV_HM_BBar"){
      PurityAntiProton = 0.989859;
  }
  else{
      printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
      PurityAntiProton = 1.0;
  }

  //following my lambda pars with the 3 possible modifications
  //for the proton:
  //0 = primary
  //1 = from Lambda
  //2 = other feeddown (flat)
  //3 = missidentified
  //const unsigned NumChannels_p = 4;
  //if(Purities){delete [] Purities; Purities = new double [NumChannels_p];}
  Purities[0] = PurityAntiProton;
  Purities[1] = PurityAntiProton;
  Purities[2] = PurityAntiProton;
  Purities[3] = 1.-PurityAntiProton;
}

void DLM_CommonAnaFunctions::GetPurities_p(const TString& DataSample, const int& Variation, double* Purities){
    double PurityProton;
    if(DataSample=="pp13TeV_MB_Run2paper"){
         PurityProton = 0.989859;
//        PurityProton = 0.9943;
    }
    else if(DataSample=="pp13TeV_HM_March19"){
        PurityProton = 0.9943;
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        PurityProton = 0.984265;
    }
    else if(DataSample=="pp13TeV_MB_BBar"){
         PurityProton = 0.989859;
        //PurityProton = 0.9943;
    }
    else if(DataSample=="pp13TeV_HM_BBar"){
        PurityProton = 0.9943;
    }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        PurityProton = 1.0;
    }

    //following my lambda pars with the 3 possible modifications
    //for the proton:
    //0 = primary
    //1 = from Lambda
    //2 = other feeddown (flat)
    //3 = missidentified
    //const unsigned NumChannels_p = 4;
    //if(Purities){delete [] Purities; Purities = new double [NumChannels_p];}
    Purities[0] = PurityProton;
    Purities[1] = PurityProton;
    Purities[2] = PurityProton;
    Purities[3] = 1.-PurityProton;
}

void DLM_CommonAnaFunctions::GetPurities_AL(const TString& DataSample, const int& Variation, double* Purities){
  double PurityAntiLambda;
  if(DataSample=="pp13TeV_MB_Run2paper"){
      PurityAntiLambda = 0.96768;
  }else if(DataSample=="pp13TeV_MB_BBar"){
        // PurityAntiLambda = 0.96768;
        PurityAntiLambda = 0.971;
    }
    else if(DataSample=="pp13TeV_HM_BBar"){
          // PurityAntiLambda = 0.96768;
          PurityAntiLambda = 0.971;
      }
  else{
      printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
      PurityAntiLambda = 1.0;
  }

  //for the Lambda:
  //0 = primary
  //1 = from Sigma0
  //2 = from Xim
  //3 = from Xi0
  //4 = missidentified
  Purities[0] = PurityAntiLambda;
  Purities[1] = PurityAntiLambda;
  Purities[2] = PurityAntiLambda;
  Purities[3] = PurityAntiLambda;
  Purities[4] = 1.-PurityAntiLambda;
}

void DLM_CommonAnaFunctions::GetPurities_L(const TString& DataSample, const int& Variation, double* Purities){
    double PurityLambda;
    if(DataSample=="pp13TeV_MB_Run2paper"){
        PurityLambda = 0.96768;
    }
    else if(DataSample=="pp13TeV_HM_March19"){
        printf("\033[1;33mWARNING:\033[0m pp13TeV_HM_March19 is not available yet!\n");
        PurityLambda = 0.96768;
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        PurityLambda = 0.937761;
    }
    else if(DataSample=="pp13TeV_MB_BBar"){
        PurityLambda = 0.96768;
    }
    else if(DataSample=="pp13TeV_HM_BBar"){
        PurityLambda = 0.96768;
    }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        PurityLambda = 1.0;
    }

    //for the Lambda:
    //0 = primary
    //1 = from Sigma0
    //2 = from Xim
    //3 = from Xi0
    //4 = missidentified
    Purities[0] = PurityLambda;
    Purities[1] = PurityLambda;
    Purities[2] = PurityLambda;
    Purities[3] = PurityLambda;
    Purities[4] = 1.-PurityLambda;
}

void DLM_CommonAnaFunctions::GetPurities_Xim(const TString& DataSample, const int& Variation, double* Purities){
    double PurityXim;
    if(DataSample=="pp13TeV_MB_Run2paper"){
        PurityXim = 0.956;
    }
    else if(DataSample=="pp13TeV_HM_March19"){
        printf("\033[1;33mWARNING:\033[0m pp13TeV_HM_March19 is not available yet!\n");
        PurityXim = 0.956;
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        PurityXim = 0.88;
    }
    else if(DataSample=="pp13TeV_MB_BBar"){
        PurityXim = 0.956;
    }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        PurityXim = 1.0;
    }

    //0 = primary
    //1 = from Xi-(1530)
    //2 = from Xi0(1530)
    //3 = from Omega
    //4 = missidentified
    Purities[0] = PurityXim;
    Purities[1] = PurityXim;
    Purities[2] = PurityXim;
    Purities[3] = PurityXim;
    Purities[4] = 1.-PurityXim;
}

//no variations are possible
void DLM_CommonAnaFunctions::GetFractions_Ap(const TString& DataSample, const int& Variation, double* Fractions){
  //following my lambda pars with the 3 possible modifications
  //for the proton:
  //0 = primary
  //1 = from Lambda
  //2 = other feeddown (flat)
  //3 = missidentified
  double Modify_pAp;
  switch(Variation){
      case 0 : Modify_pAp=1.; break;
      case 1 : Modify_pAp=0.8; break;
      case 2 : Modify_pAp=1.2; break;
      default : Modify_pAp=1; break;
  }
  double pAp_f0;//primary protons
  double pAp_f1;//fraction of Lambda
  if(DataSample=="pp13TeV_MB_BBar"){
      // pAp_f0 = 0.8802;
      // pAp_f1 = 0.0843;
      // pAp_f0 = 0.845975;//calculations from 24.06.2019
      // pAp_f1 = 0.100644;
      pAp_f0 = 0.833596;//calculations from 29.06.2019 TFractionFitter
      pAp_f1 = 0.116024;
  }else if(DataSample=="pp13TeV_HM_BBar"){
        // pAp_f0 = 0.8802;
        // pAp_f1 = 0.0843;
        // pAp_f0 = 0.845975;//calculations from 24.06.2019
        // pAp_f1 = 0.100644;
//        pAp_f0 = 0.833596;//calculations from 29.06.2019 TFractionFitter
//        pAp_f1 = 0.116024;
        pAp_f0 = 0.82266;//calculations from Andi update on MC 5% 27.09.2019 TFractionFitter
        pAp_f1 = 0.124465;
    }
  else{
      printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
      pAp_f0 = 1.0;
      pAp_f1 = 0.0;
  }
  double arrayPercLamAntiProton = pAp_f1/(1.-pAp_f0)*Modify_pAp;
  Fractions[0] = pAp_f0;
  Fractions[1] = (1.-pAp_f0)*(arrayPercLamAntiProton);
  Fractions[2] = (1.-pAp_f0)*(1.-arrayPercLamAntiProton);
  Fractions[3] = 1.;
  // printf("--------------------------------------------------------------\n");
  // printf("---------------DATASample analysed = %s-----------------------\n",DataSample.Data());
  printf("-----------------------ANTIPROTONS----------------------------\n");
      printf("Modify_pAp = %.1f ---- Modify_pAp = %.1f\n",Modify_pAp,Modify_pAp);
  printf("Fractions of primaries = %.3f\n", Fractions[0]);
  printf("Fractions of Sec.from Λ = %.3f\n", Fractions[1]);
  printf("Fractions of Sec.from Σ+ = %.3f\n", Fractions[2]);

}
void DLM_CommonAnaFunctions::GetFractions_p(const TString& DataSample, const int& Variation, double* Fractions){
    //following my lambda pars with the 3 possible modifications
    //for the proton:
    //0 = primary
    //1 = from Lambda
    //2 = other feeddown (flat)
    //3 = missidentified
    double Modify_pp;
    switch(Variation){
        case 0 : Modify_pp=1.; break;
        case 1 : Modify_pp=0.8; break;
        case 2 : Modify_pp=1.2; break;
        default : Modify_pp=1.; break;
    }
    double pp_f0;//primary protons
    double pp_f1;//fraction of Lambda
    if(DataSample=="pp13TeV_MB_Run2paper"){
        pp_f0 = 0.87397;
        pp_f1 = 0.0882211;
    }
    else if(DataSample=="pp13TeV_HM_March19"){
        printf("\033[1;33mWARNING:\033[0m pp13TeV_HM_March19 CROSS CHECK pp_f0 and pp_f1!\n");
        pp_f0 = 0.873;
        pp_f1 = 0.0898;
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        pp_f0 = 0.862814;
        pp_f1 = 0.09603;
    }
    else if(DataSample=="pp13TeV_MB_BBar"){
        // pp_f0 = 0.8797;
        // pp_f1 = 0.0848;
        // pp_f0 = 0.854938;//calculations from 24.06.2019
        // pp_f1 = 0.101178;
         pp_f0 = 0.83531;//calculations from 29.07.2019 from TFractionFitter
         pp_f1 = 0.114904;
    }
    else if(DataSample=="pp13TeV_HM_BBar"){
//         pp_f0 = 0.83531;//calculations from 29.07.2019 from TFractionFitter
//         pp_f1 = 0.114904;
         pp_f0 = 0.82164;//calculations from Andi update on MC 5% 27.09.2019 TFractionFitter
         pp_f1 = 0.125253;
    }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        pp_f0 = 1.0;
        pp_f1 = 0.0;
    }
    double arrayPercLamProton = pp_f1/(1.-pp_f0)*Modify_pp;
    Fractions[0] = pp_f0;
    Fractions[1] = (1.-pp_f0)*(arrayPercLamProton);
    Fractions[2] = (1.-pp_f0)*(1.-arrayPercLamProton);
    Fractions[3] = 1.;
    // printf("--------------------------------------------------------------\n");
    // printf("---------------DATASample analysed = %s-----------------------\n",DataSample.Data());
    printf("-----------------------PROTONS--------------------------------\n");
      printf("Modify_pp = %.1f ---- Modify_pp = %.1f\n",Modify_pp,Modify_pp);
     printf("Fractions of primaries = %.3f\n", Fractions[0]);
     printf("Fractions of Sec.from Λ = %.3f\n", Fractions[1]);
     printf("Fractions of Sec.from Σ+ = %.3f\n", Fractions[2]);

}
//Variation -> use the two digits
//first digit->Modify_SigL
//second digit->Modify_XiL
//0 -> default; 1 = -20%; 2 = +20%
void DLM_CommonAnaFunctions::GetFractions_AL(const TString& DataSample, const int& Variation, double* Fractions){
  double Modify_SigAL;
  double Modify_XiAL;
  switch(Variation%10){
      case 0 : Modify_SigAL=1.; break;
      case 1 : Modify_SigAL=0.8;break;
      case 2 : Modify_SigAL=1.2; break;
      default : Modify_SigAL=1; break;
  }
  switch(Variation/10){
      case 0 : Modify_XiAL=1.; break;
      case 1 : Modify_XiAL=0.8;break;
      case 2 : Modify_XiAL=1.2; break;
      default : Modify_XiAL=1; break;
  }
  double pAL_f0;//fraction of primary Lambdas
  double pAL_f1;//fraction of Sigma0
  double pAL_f2;//fractions of Xi0/m
  if(DataSample=="pp13TeV_MB_BBar"){
      // pAL_f0 = 0.5912;
      // pAL_f1 = 0.1971;
      // pAL_f2 = 0.1059;

      pAL_f0 = 0.577143;//23.07.2019 from FractionFitter
      pAL_f1 = 0.192381;
      pAL_f2 = 0.102758;
  }else if(DataSample=="pp13TeV_HM_BBar"){
        // pAL_f0 = 0.5912;
        // pAL_f1 = 0.1971;
        // pAL_f2 = 0.1059;

        pAL_f0 = 0.577143;//23.07.2019 from FractionFitter
        pAL_f1 = 0.192381;
        pAL_f2 = 0.102758;
    }
  else{
      printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
      pAL_f0 = 1.0;
      pAL_f1 = 0.0;
      pAL_f2 = 0.0;
  }

  double SigAntiLambdaPrimDir = pAL_f0+pAL_f1;
  double arrayPercSigAntiLambda=pAL_f1/pAL_f0*Modify_SigAL;
  double arrayPercXiAntiLambda=pAL_f2/(1.-pAL_f0-pAL_f1)*Modify_XiAL;
  double FracOfAntiLambda = 1./(1.+arrayPercSigAntiLambda);
  //0 is primary
  //1 is from Sigma0
  //2 is is from Xim
  //3 is is the flat feeddown
  //4 is for the missid
  Fractions[0] = SigAntiLambdaPrimDir*FracOfAntiLambda;
  Fractions[1] = SigAntiLambdaPrimDir*(1.-FracOfAntiLambda);
  Fractions[2] = (1.-SigAntiLambdaPrimDir)*(arrayPercXiAntiLambda);
  Fractions[3] = 1.-Fractions[0]-Fractions[1]-Fractions[2];
  Fractions[4] = 1.;
  // printf("--------------------------------------------------------------\n");
  // printf("---------------DATASample analysed = %s-----------------------\n",DataSample.Data());
  printf("-----------------------ANTILAMBDAS FRACTIONS----------------------------\n");
      printf("Modify_SigL = %.1f ---- Modify_XiL = %.1f\n",Modify_SigAL,Modify_XiAL);
  printf("Fractions of primaries = %.3f\n", Fractions[0]);
  printf("Fractions of Sec.from Σ0 = %.3f\n", Fractions[1]);
  printf("Fractions of Sec.from Ξm0 = %.3f\n", Fractions[2]);
  std::cout<<"Fractions for Lambdabar:\n"<<std::endl;
  std::cout<<"Primaries = "<<Fractions[0]<<std::endl;
  std::cout<<"From Sigma0 = "<<Fractions[1]<<std::endl;
  std::cout<<"From Xim = "<<Fractions[2]<<std::endl;
  std::cout<<"Flat feed = "<<Fractions[3]<<std::endl;
  std::cout<<"MisID = "<<Fractions[4]<<std::endl;
  std::cout<<"----------------------------\n"<<std::endl;

}

void DLM_CommonAnaFunctions::GetFractions_L(const TString& DataSample, const int& Variation, double* Fractions){
    double Modify_SigL;
    double Modify_XiL;
    switch(Variation%10){
        case 0 : Modify_SigL=1.; break;
        case 1 : Modify_SigL=0.8;break;
        case 2 : Modify_SigL=1.2; break;
        default : Modify_SigL=1.; break;
    }
    switch(Variation/10){
        case 0 : Modify_XiL=1.; break;
        case 1 : Modify_XiL=0.8;break;
        case 2 : Modify_XiL=1.2; break;
        default : Modify_XiL=1.; break;
    }
    double pL_f0;//fraction of primary Lambdas
    double pL_f1;//fraction of Sigma0
    double pL_f2;//fractions of Xi0/m
    if(DataSample=="pp13TeV_MB_Run2paper"){
        pL_f0 = 0.601008;
        pL_f1 = 0.200336;
        pL_f2 = 0.099328;
    }
    else if(DataSample=="pp13TeV_HM_March19"){
        pL_f0 = 0.576066;
        pL_f1 = 0.192022;
        pL_f2 = 0.115956;
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        pL_f0 = 0.521433;
        pL_f1 = 0.173811;
        pL_f2 = 0.152378;
    }
    else if(DataSample=="pp13TeV_MB_BBar"){
        // pL_f0 = 0.5941;
        // pL_f1 = 0.198;
        // pL_f2 = 0.1039;

        pL_f0 = 0.569624;// 23.07.2019 from FractionFitter
        pL_f1 = 0.189875;
        pL_f2 = 0.106327;
    }
    else if(DataSample=="pp13TeV_HM_BBar"){
            // pL_f0 = 0.5941;
            // pL_f1 = 0.198;
            // pL_f2 = 0.1039;

            pL_f0 = 0.569624;// 23.07.2019 from FractionFitter
            pL_f1 = 0.189875;
            pL_f2 = 0.106327;
        }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        pL_f0 = 1.0;
        pL_f1 = 0.0;
        pL_f2 = 0.0;
    }
    double SigLambdaPrimDir = pL_f0+pL_f1;
    double arrayPercSigLambda=pL_f1/pL_f0*Modify_SigL;
    double arrayPercXiLambda=pL_f2/(1.-pL_f0-pL_f1)*Modify_XiL;
    double FracOfLambda = 1./(1.+arrayPercSigLambda);
    //0 is primary
    //1 is from Sigma0
    //2 is is from Xim
    //3 is is the flat feeddown
    //4 is for the missid
    Fractions[0] = SigLambdaPrimDir*FracOfLambda;
    Fractions[1] = SigLambdaPrimDir*(1.-FracOfLambda);
    Fractions[2] = (1.-SigLambdaPrimDir)*(arrayPercXiLambda);
    Fractions[3] = 1.-Fractions[0]-Fractions[1]-Fractions[2];
    Fractions[4] = 1.;
    // printf("--------------------------------------------------------------\n");
    // printf("---------------DATASample analysed = %s-----------------------\n",DataSample.Data());
    printf("-----------------------LAMBDAS FRACTIONS--------------------------------\n");
    printf("Modify_SigL = %.1f ---- Modify_XiL = %.1f\n",Modify_SigL,Modify_XiL);
    printf("Fractions of primaries = %.3f\n", Fractions[0]);
    printf("Fractions of Sec.from Σ0 = %.3f\n", Fractions[1]);
    printf("Fractions of Sec.from Ξm0 = %.3f\n", Fractions[2]);
    std::cout<<"Fractions for Lambda:\n"<<std::endl;
    std::cout<<"Primaries = "<<Fractions[0]<<std::endl;
    std::cout<<"From Sigma0 = "<<Fractions[1]<<std::endl;
    std::cout<<"From Xim = "<<Fractions[2]<<std::endl;
    std::cout<<"Flat feed = "<<Fractions[3]<<std::endl;
    std::cout<<"MisID = "<<Fractions[4]<<std::endl;
    std::cout<<"----------------------------\n"<<std::endl;

}
void DLM_CommonAnaFunctions::GetFractions_Xim(const TString& DataSample, const int& Variation, double* Fractions){
    //0 = primary
    //1 = from Xi-(1530)
    //2 = from Omega
    //3 = flat
    //4 = missidentified

    //ratio Xi-(1530) to Xi-
    const double Xim1530_to_Xim = 0.32*(1./3.);
    //ratio Xi0(1530) to Xi0 (n=neutral)
    //const double Xin1530_to_Xim = 0.32*(2./3.);
    const double Omegam_to_Xim = 0.1;
    const double OmegamXim_BR = 0.086;
    //the ratios that we have for Xis are referred to the total number of Xi particles (which already include all contributions)
    //hence Xi1530_to_Xi indeed is simply the number of Xis that stem from a Xi1530
    Fractions[0] = 1.-3.*Xim1530_to_Xim-Omegam_to_Xim*OmegamXim_BR;
    Fractions[1] = Xim1530_to_Xim;
    Fractions[2] = Omegam_to_Xim*OmegamXim_BR;
    Fractions[3] = 1.-Fractions[2]-Fractions[1];
    Fractions[4] = 1.;
}

// Single contributions of Lambda parameters
//0 is primary
//1 is pL->pp
//2 is flat feed
//3 is missid

void DLM_CommonAnaFunctions::SetUpLambdaPars_pAp(const TString& DataSample, const int& Variation_p, const int& Variation_Ap, double* lambda_pars){
  double Purities_Ap[4];
  double Fraction_Ap[4];
  double Purities_p[4];
  double Fraction_p[4];
  GetPurities_Ap(DataSample,Variation_Ap,Purities_Ap);
  GetFractions_Ap(DataSample,Variation_Ap,Fraction_Ap);
  GetPurities_p(DataSample,Variation_p,Purities_p);
  GetFractions_p(DataSample,Variation_p,Fraction_p);
  lambda_pars[0] = Purities_Ap[0]*Fraction_Ap[0]*Purities_p[0]*Fraction_p[0];
  lambda_pars[1] = Purities_Ap[0]*Fraction_Ap[0]*Purities_p[1]*Fraction_p[1]+Purities_Ap[1]*Fraction_Ap[1]*Purities_p[0]*Fraction_p[0];
  lambda_pars[3] = Purities_Ap[0]*Fraction_Ap[0]*Purities_p[3]*Fraction_Ap[3]+Purities_Ap[3]*Fraction_Ap[3]*Purities_p[0]*Fraction_Ap[0];
  lambda_pars[2] = 1.-lambda_pars[3]-lambda_pars[1]-lambda_pars[0];
  std::cout<<"Lambda parameters for ppbar:\n"<<std::endl;
  std::cout<<"Primaries = "<<lambda_pars[0]<<std::endl;
  std::cout<<"From Lambda = "<<lambda_pars[1]<<std::endl;
  std::cout<<"Flat feed = "<<lambda_pars[2]<<std::endl;
  std::cout<<"MisID = "<<lambda_pars[3]<<std::endl;
  std::cout<<"-----------------------------\n"<<std::endl;
  std::cout<<"-----------------------------\n"<<std::endl;
  std::cout<<"TOT.ppbar= "<<lambda_pars[0]+lambda_pars[1]+lambda_pars[2]+lambda_pars[3]<<std::endl;

}

//0 is primary
//1 is pL->pp
//2 is flat feed
//3 is missid


void DLM_CommonAnaFunctions::SetUpLambdaPars_pp(const TString& DataSample, const int& Variation_p, double* lambda_pars){
    double Purities_p[4];
    double Fraction_p[4];
    GetPurities_p(DataSample,Variation_p,Purities_p);
    GetFractions_p(DataSample,Variation_p,Fraction_p);
    lambda_pars[0] = Purities_p[0]*Fraction_p[0]*Purities_p[0]*Fraction_p[0];
    lambda_pars[1] = Purities_p[0]*Fraction_p[0]*Purities_p[1]*Fraction_p[1]*2.;
    lambda_pars[3] = Purities_p[0]*Fraction_p[0]*Purities_p[3]*Purities_p[3]*2.;
//    lambda_pars[3] = (Purities_p[0]+Purities_p[0]+Purities_p[3])*Purities_p[3];
    lambda_pars[2] = 1.-lambda_pars[3]-lambda_pars[1]-lambda_pars[0];

    //double SUM=0;
    //for(unsigned uLam=0; uLam<4; uLam++){
    //    printf("λ(pp)_%u = %.1f\n",uLam,lambda_pars[uLam]*100.);
    //    SUM+=lambda_pars[uLam]*100.;
    //}
    //printf("SUM: %.1f\n------------\n",SUM);
}
//0 is primary
//1 is pSigma0->pL
//2 is pXim->pL
//3 is the flat feeddown
//4 is missid
void DLM_CommonAnaFunctions::SetUpLambdaPars_pAL(const TString& DataSample, const int& Variation_p, const int& Variation_AL, double* lambda_pars){
  double Purities_p[4];
  double Fraction_p[4];
  double Purities_AL[5];
  double Fraction_AL[5];
  GetPurities_p(DataSample,Variation_p,Purities_p);
  GetFractions_p(DataSample,Variation_p,Fraction_p);
  GetPurities_AL(DataSample,Variation_AL,Purities_AL);
  GetFractions_AL(DataSample,Variation_AL,Fraction_AL);



  lambda_pars[0] =    Purities_p[0]*Fraction_p[0]*Purities_AL[0]*Fraction_AL[0];//primaries
  lambda_pars[1] =    Purities_p[0]*Fraction_p[0]*Purities_AL[1]*Fraction_AL[1];//L from Sigma0
  lambda_pars[2] =    Purities_p[0]*Fraction_p[0]*Purities_AL[2]*Fraction_AL[2];//L from Xim
  lambda_pars[3] =    Purities_p[1]*Fraction_p[1]*Purities_AL[0]*Fraction_AL[0];//p from Lambda
  lambda_pars[5] =    Purities_p[0]*Fraction_p[0]*Purities_AL[4]*Fraction_AL[4]+ Purities_p[3]*Fraction_p[3]*Purities_AL[0]*Fraction_AL[0]+
                      2*(Purities_p[3]*Fraction_p[3]*Purities_AL[2]*Fraction_AL[2])+
                      Purities_p[3]*Fraction_p[3]*Purities_AL[1]*Fraction_AL[1]+
                      Purities_p[1]*Fraction_p[1]*Purities_AL[4]*Fraction_AL[4]+
                      Purities_p[2]*Fraction_p[2]*Purities_AL[4]*Fraction_AL[4]+
                      Purities_p[3]*Fraction_p[3]*Purities_AL[4]*Fraction_AL[4];
  lambda_pars[4] =    1.-lambda_pars[0]-lambda_pars[1]-lambda_pars[2]-lambda_pars[3]-lambda_pars[5];
  std::cout<<"Lambda parameters for pLbar:\n"<<std::endl;
  std::cout<<"Primaries = "<<lambda_pars[0]<<std::endl;
  std::cout<<"From Sigma0 = "<<lambda_pars[1]<<std::endl;
  std::cout<<"From Xim = "<<lambda_pars[2]<<std::endl;
  std::cout<<"From Lambda = "<<lambda_pars[3]<<std::endl;
  std::cout<<"Flat feed = "<<lambda_pars[4]<<std::endl;
  std::cout<<"MisID = "<<lambda_pars[5]<<std::endl;
  std::cout<<"-----------------------------\n"<<std::endl;
  std::cout<<"-----------------------------\n"<<std::endl;
  std::cout<<"TOT.pLbar= "<<lambda_pars[0]+lambda_pars[1]+lambda_pars[2]+lambda_pars[3]+lambda_pars[4]+lambda_pars[5]<<std::endl;

}

void DLM_CommonAnaFunctions::SetUpLambdaPars_pL(const TString& DataSample, const int& Variation_p, const int& Variation_L, double* lambda_pars){
    double Purities_p[4];
    double Fraction_p[4];
    double Purities_L[5];
    double Fraction_L[5];
    GetPurities_p(DataSample,Variation_p,Purities_p);
    GetFractions_p(DataSample,Variation_p,Fraction_p);
    GetPurities_L(DataSample,Variation_L,Purities_L);
    GetFractions_L(DataSample,Variation_L,Fraction_L);
    lambda_pars[0] =    Purities_p[0]*Fraction_p[0]*Purities_L[0]*Fraction_L[0];
    lambda_pars[1] =    Purities_p[0]*Fraction_p[0]*Purities_L[1]*Fraction_L[1];
    lambda_pars[2] =    Purities_p[0]*Fraction_p[0]*Purities_L[2]*Fraction_L[2];
    lambda_pars[4] =    Purities_p[0]*Fraction_p[0]*Purities_L[4]*Fraction_L[4]+ Purities_p[3]*Fraction_p[3]*Purities_L[0]*Fraction_L[0];
    lambda_pars[3] =    1.-lambda_pars[0]-lambda_pars[1]-lambda_pars[2]-lambda_pars[4];

    //double SUM=0;
    //for(unsigned uLam=0; uLam<5; uLam++){
    //    printf("λ(pΛ)_%u = %.1f\n",uLam,lambda_pars[uLam]*100.);
    //    SUM+=lambda_pars[uLam]*100.;
    //}
    //printf("SUM: %.1f\n------------\n",SUM);
}
//0 is primary
//1 is LSigma0->LL
//2 is LXim->LL
//3 is the flat feeddown
//4 is missid
void DLM_CommonAnaFunctions::SetUpLambdaPars_LAL(const TString& DataSample, const int& Variation_L, const int& Variation_AL, double* lambda_pars){
  double Purities_L[5];
  double Fraction_L[5];
  double Purities_AL[5];
  double Fraction_AL[5];
  GetPurities_L(DataSample,Variation_L,Purities_L);
  GetFractions_L(DataSample,Variation_L,Fraction_L);
  GetPurities_AL(DataSample,Variation_AL,Purities_AL);
  GetFractions_AL(DataSample,Variation_AL,Fraction_AL);
  lambda_pars[0] =    Purities_L[0]*Fraction_L[0]*Purities_AL[0]*Fraction_AL[0];
  lambda_pars[1] =    Purities_L[0]*Fraction_L[0]*Purities_AL[1]*Fraction_AL[1]+Purities_L[1]*Fraction_L[1]*Purities_AL[0]*Fraction_AL[0];
  lambda_pars[2] =    Purities_L[0]*Fraction_L[0]*Purities_AL[2]*Fraction_AL[2]+Purities_L[2]*Fraction_L[2]*Purities_AL[0]*Fraction_AL[0];
  lambda_pars[4] =    Purities_L[0]*Fraction_L[0]*Purities_AL[4]*Fraction_AL[4]+Purities_L[4]*Fraction_L[4]*Purities_AL[0]*Fraction_AL[0];
  lambda_pars[3] =    1.-lambda_pars[0]-lambda_pars[1]-lambda_pars[2]-lambda_pars[4];
  std::cout<<"-----------------------------\n"<<std::endl;
  std::cout<<"-----------------------------\n"<<std::endl;
  std::cout<<"Lambda parameters for LLbar:\n"<<std::endl;
  std::cout<<"Primaries = "<<lambda_pars[0]<<std::endl;
  std::cout<<"From Sigma0 = "<<lambda_pars[1]<<std::endl;
  std::cout<<"From Xim = "<<lambda_pars[2]<<std::endl;
  std::cout<<"Flat feed = "<<lambda_pars[3]<<std::endl;
  std::cout<<"TOTAL Flat feed = "<<lambda_pars[1]+lambda_pars[2]+lambda_pars[3]<<std::endl;
  std::cout<<"MisID = "<<lambda_pars[4]<<std::endl;
  std::cout<<"-----------------------------\n"<<std::endl;
  std::cout<<"TOT.LLbar = "<<lambda_pars[0]+lambda_pars[1]+lambda_pars[2]+lambda_pars[3]+lambda_pars[4]<<std::endl;
  //
  // printf("-------debug 0-------------------\n");


}
// Single contributions of Lambda parameters for crosscheck with RUN2 BB

// Single contributions of Lambda parameters
//0 is primary
//1 is pL->pp
//2 is flat feed
//3 is missid

void DLM_CommonAnaFunctions::SetUpLambdaParsContrib_pAp(const TString& DataSample, const int& Variation_Ap, double* lambda_pars){
  double Purities_Ap[4];
  double Fraction_Ap[4];
  double Purities_p[4];
  double Fraction_p[4];
  GetPurities_Ap(DataSample,Variation_Ap,Purities_Ap);
  GetFractions_Ap(DataSample,Variation_Ap,Fraction_Ap);
  GetPurities_p(DataSample,Variation_Ap,Purities_p);
  GetFractions_p(DataSample,Variation_Ap,Fraction_p);
  lambda_pars[0] = Purities_Ap[0]*Fraction_Ap[0]*Purities_p[0]*Fraction_p[0];//Primaries
  lambda_pars[1] = Purities_Ap[0]*Fraction_Ap[0]*Purities_p[1]*Fraction_p[1]
  +Purities_Ap[1]*Fraction_Ap[1]*Purities_p[0]*Fraction_p[0];//from Lambda/AntiLambda
  lambda_pars[2] = Purities_Ap[1]*Fraction_Ap[1]*Purities_p[1]*Fraction_p[1];//both from LAMBDA
  lambda_pars[3] = Purities_Ap[0]*Fraction_Ap[0]*Purities_p[2]*Fraction_p[2]
  +Purities_Ap[2]*Fraction_Ap[2]*Purities_p[0]*Fraction_p[0];//from Sigma+/AntiSigma+
  lambda_pars[4] = Purities_Ap[2]*Fraction_Ap[2]*Purities_p[2]*Fraction_p[2];//both from Sigma+
  lambda_pars[5] = Purities_Ap[1]*Fraction_Ap[1]*Purities_p[2]*Fraction_p[2]
  +Purities_Ap[2]*Fraction_Ap[2]*Purities_p[1]*Fraction_p[1];//From Lambda and Sigma+
  lambda_pars[6] = Purities_Ap[0]*Fraction_Ap[0]*Purities_p[3]*Fraction_p[3]
  +Purities_Ap[3]*Fraction_Ap[3]*Purities_p[0]*Fraction_p[0];//one misID and one primary
  lambda_pars[7] = Purities_Ap[1]*Fraction_Ap[1]*Purities_p[3]*Fraction_p[3]
  +Purities_Ap[3]*Fraction_Ap[3]*Purities_p[1]*Fraction_p[1];//one misID and one from Lambda
  lambda_pars[8] = Purities_Ap[2]*Fraction_Ap[2]*Purities_p[3]*Fraction_p[3]
  +Purities_Ap[3]*Fraction_Ap[3]*Purities_p[2]*Fraction_p[2];//one misID and one from Sigma+
  lambda_pars[9] = Purities_Ap[3]*Fraction_Ap[3]*Purities_p[3]*Fraction_p[3];//both misID

  // lambda_pars[3] = (Purities_Ap[0]+Purities_Ap[0]+Purities_Ap[3])*Purities_Ap[3];
  // lambda_pars[2] = 1.-lambda_pars[3]-lambda_pars[1]-lambda_pars[0];
}

void DLM_CommonAnaFunctions::SetUpLambdaParsContrib_pp(const TString& DataSample, const int& Variation_p, double* lambda_pars){
    double Purities_p[4];
    double Fraction_p[4];
    GetPurities_p(DataSample,Variation_p,Purities_p);
    GetFractions_p(DataSample,Variation_p,Fraction_p);
    lambda_pars[0] = Purities_p[0]*Fraction_p[0]*Purities_p[0]*Fraction_p[0];//Primaries
    lambda_pars[1] = Purities_p[0]*Fraction_p[0]*Purities_p[1]*Fraction_p[1]+Purities_p[1]*Fraction_p[1]*Purities_p[0]*Fraction_p[0];//from Lambda/AntiLambda
    lambda_pars[2] = Purities_p[1]*Fraction_p[1]*Purities_p[1]*Fraction_p[1];//both from LAMBDA
    lambda_pars[3] = Purities_p[0]*Fraction_p[0]*Purities_p[2]*Fraction_p[2]+Purities_p[2]*Fraction_p[2]*Purities_p[0]*Fraction_p[0];//from Sigma+/AntiSigma+
    lambda_pars[4] = Purities_p[2]*Fraction_p[2]*Purities_p[2]*Fraction_p[2];//both from Sigma+
    lambda_pars[5] = Purities_p[1]*Fraction_p[1]*Purities_p[2]*Fraction_p[2]+Purities_p[2]*Fraction_p[2]*Purities_p[1]*Fraction_p[1];//From Lambda and Sigma+
    lambda_pars[6] = Purities_p[0]*Fraction_p[0]*Purities_p[3]*Fraction_p[3]+Purities_p[3]*Fraction_p[3]*Purities_p[0]*Fraction_p[0];//one misID and one primary
    lambda_pars[7] = Purities_p[1]*Fraction_p[1]*Purities_p[3]*Fraction_p[3]+Purities_p[3]*Fraction_p[3]*Purities_p[1]*Fraction_p[1];//one misID and one from Lambda
    lambda_pars[8] = Purities_p[2]*Fraction_p[2]*Purities_p[3]*Fraction_p[3]+Purities_p[3]*Fraction_p[3]*Purities_p[2]*Fraction_p[2];//one misID and one from Sigma+
    lambda_pars[9] = Purities_p[3]*Fraction_p[3]*Purities_p[3]*Fraction_p[3];//both misID


    //double SUM=0;
    //for(unsigned uLam=0; uLam<4; uLam++){
    //    printf("λ(pp)_%u = %.1f\n",uLam,lambda_pars[uLam]*100.);
    //    SUM+=lambda_pars[uLam]*100.;
    //}
    //printf("SUM: %.1f\n------------\n",SUM);
}

//0 is primary
//1 is pSigma0->pL
//2 is pXim->pL
//3 is the flat feeddown
//4 is missid
void DLM_CommonAnaFunctions::SetUpLambdaParsContrib_pAL(const TString& DataSample, const int& Variation_p, const int& Variation_AL, double* lambda_pars){
  double Purities_p[4];
  double Fraction_p[4];
  double Purities_AL[5];
  double Fraction_AL[5];
  GetPurities_p(DataSample,Variation_p,Purities_p);
  GetFractions_p(DataSample,Variation_p,Fraction_p);
  GetPurities_AL(DataSample,Variation_AL,Purities_AL);
  GetFractions_AL(DataSample,Variation_AL,Fraction_AL);
  lambda_pars[0] =    Purities_p[0]*Fraction_p[0]*Purities_AL[0]*Fraction_AL[0];//primaries
  lambda_pars[1] =    Purities_p[0]*Fraction_p[0]*Purities_AL[2]*Fraction_AL[2];//p and Lambda from Xim == Xi0
  lambda_pars[2] =    Purities_p[0]*Fraction_p[0]*Purities_AL[2]*Fraction_AL[2];//p and Lambda from Xim == Xi0
  lambda_pars[3] =    Purities_p[0]*Fraction_p[0]*Purities_AL[1]*Fraction_AL[1];//p and Lambda from Sigma0
  lambda_pars[4] =    Purities_p[1]*Fraction_p[1]*Purities_AL[0]*Fraction_AL[0];//p from Lambda and prim Lambda
  lambda_pars[5] =    Purities_p[1]*Fraction_p[1]*Purities_AL[2]*Fraction_AL[2];//p from Lambda and Lambda from Xim
  lambda_pars[6] =    Purities_p[1]*Fraction_p[1]*Purities_AL[2]*Fraction_AL[2];//p from Lambda and Lambda from Xi0
  lambda_pars[7] =    Purities_p[1]*Fraction_p[1]*Purities_AL[1]*Fraction_AL[1];//p from Lambda and Lambda from Sigma0
  lambda_pars[8] =    Purities_p[2]*Fraction_p[2]*Purities_AL[0]*Fraction_AL[0];//p from Sigma+ and prim Lambda
  lambda_pars[9] =    Purities_p[2]*Fraction_p[2]*Purities_AL[2]*Fraction_AL[2];//p from Sigma+ and Lambda from Xim
  lambda_pars[10] =    Purities_p[2]*Fraction_p[2]*Purities_AL[2]*Fraction_AL[2];//p from Sigma+ and Lambda from Xi0
  lambda_pars[11] =    Purities_p[2]*Fraction_p[2]*Purities_AL[1]*Fraction_AL[1];//p from Sigma+ and Lambda from Sigma0

  lambda_pars[12] =    Purities_p[3]*Fraction_p[3]*Purities_AL[0]*Fraction_AL[0];//p misID and prim Lambda
  lambda_pars[13] =    Purities_p[3]*Fraction_p[3]*Purities_AL[2]*Fraction_AL[2];//p misID and Lambda from Xim
  lambda_pars[14] =    Purities_p[3]*Fraction_p[3]*Purities_AL[2]*Fraction_AL[2];//p misID and Lambda from Xi0
  lambda_pars[15] =    Purities_p[3]*Fraction_p[3]*Purities_AL[1]*Fraction_AL[1];//p misID and Lambda from Sigma0

  lambda_pars[16] =    Purities_p[0]*Fraction_p[0]*Purities_AL[4]*Fraction_AL[4];//p prim and misID Lambda
  lambda_pars[17] =    Purities_p[1]*Fraction_p[1]*Purities_AL[4]*Fraction_AL[4];//p from Lambda and misID Lambda
  lambda_pars[18] =    Purities_p[2]*Fraction_p[2]*Purities_AL[4]*Fraction_AL[4];//p from Sigma+ and misID Lambda
  lambda_pars[19] =    Purities_p[3]*Fraction_p[3]*Purities_AL[4]*Fraction_AL[4];//both misID

}

void DLM_CommonAnaFunctions::SetUpLambdaParsContrib_LAL(const TString& DataSample, const int& Variation_L, const int& Variation_AL, double* lambda_pars){
  double Purities_L[5];
  double Fraction_L[5];
  double Purities_AL[5];
  double Fraction_AL[5];
  GetPurities_L(DataSample,Variation_L,Purities_L);
  GetFractions_L(DataSample,Variation_L,Fraction_L);
  GetPurities_AL(DataSample,Variation_AL,Purities_AL);
  GetFractions_AL(DataSample,Variation_AL,Fraction_AL);
  lambda_pars[0] =    Purities_L[0]*Fraction_L[0]*Purities_AL[0]*Fraction_AL[0];//primaries
  lambda_pars[1] =    Purities_L[0]*Fraction_L[0]*Purities_AL[1]*Fraction_AL[1]+Purities_L[1]*Fraction_L[1]*Purities_AL[0]*Fraction_AL[0];//one prim. one from Sigma0
  lambda_pars[2] =    Purities_L[1]*Fraction_L[1]*Purities_AL[1]*Fraction_AL[1];//both from Sigma0
  lambda_pars[3] =    Purities_L[0]*Fraction_L[0]*Purities_AL[2]*Fraction_AL[2]+Purities_L[2]*Fraction_L[2]*Purities_AL[0]*Fraction_AL[0];//one prim. one from Xi0
  lambda_pars[4] =    Purities_L[2]*Fraction_L[2]*Purities_AL[2]*Fraction_AL[2];//both from Xi0
  lambda_pars[5] =    Purities_L[0]*Fraction_L[0]*Purities_AL[2]*Fraction_AL[2]+Purities_L[2]*Fraction_L[2]*Purities_AL[0]*Fraction_AL[0];//one prim. one from Xi-
  lambda_pars[6] =    Purities_L[2]*Fraction_L[2]*Purities_AL[2]*Fraction_AL[2];//both from Xi-

  lambda_pars[7] =    Purities_L[1]*Fraction_L[1]*Purities_AL[2]*Fraction_AL[2]+Purities_L[2]*Fraction_L[2]*Purities_AL[1]*Fraction_AL[1];//one Sigma0 one from Xi0
  lambda_pars[8] =    Purities_L[1]*Fraction_L[1]*Purities_AL[2]*Fraction_AL[2]+Purities_L[2]*Fraction_L[2]*Purities_AL[1]*Fraction_AL[1];//one Sigma0 one from Xi-

  lambda_pars[9] =    Purities_L[2]*Fraction_L[2]*Purities_AL[2]*Fraction_AL[2]*2.;//one from Xi0 one from Xi-

  lambda_pars[10] =    Purities_L[0]*Fraction_L[0]*Purities_AL[4]*Fraction_AL[4]+Purities_L[4]*Fraction_L[4]*Purities_AL[0]*Fraction_AL[0];//one prim one misID
  lambda_pars[11] =    Purities_L[4]*Fraction_L[4]*Purities_AL[1]*Fraction_AL[1]+Purities_L[1]*Fraction_L[1]*Purities_AL[4]*Fraction_AL[4];//one misID one Sigma0
  lambda_pars[12] =    Purities_L[4]*Fraction_L[4]*Purities_AL[2]*Fraction_AL[2]+Purities_L[2]*Fraction_L[2]*Purities_AL[4]*Fraction_AL[4];//one misID one Xi0
  lambda_pars[13] =    Purities_L[4]*Fraction_L[4]*Purities_AL[2]*Fraction_AL[2]+Purities_L[2]*Fraction_L[2]*Purities_AL[4]*Fraction_AL[4];//one misID one Xi-

  lambda_pars[14] =    Purities_L[4]*Fraction_L[4]*Purities_AL[4]*Fraction_AL[4];//both misID

}


//0 is primary
//1 is from Xim1530
//2 is from Xin1530
//3
//4 is missid
void DLM_CommonAnaFunctions::SetUpLambdaPars_pXim(const TString& DataSample, const int& Variation_p, const int& Variation_Xim, double* lambda_pars){
    double Purities_p[4];
    double Fraction_p[4];
    double Purities_Xim[5];
    double Fraction_Xim[5];
    GetPurities_p(DataSample,Variation_p,Purities_p);
    GetFractions_p(DataSample,Variation_p,Fraction_p);
    GetPurities_Xim(DataSample,Variation_Xim,Purities_Xim);
    GetFractions_Xim(DataSample,Variation_Xim,Fraction_Xim);
    lambda_pars[0] =    Purities_p[0]*Fraction_p[0]*Purities_Xim[0]*Fraction_Xim[0];
    lambda_pars[1] =    Purities_p[0]*Fraction_p[0]*Purities_Xim[1]*Fraction_Xim[1];
    lambda_pars[2] =    Purities_p[0]*Fraction_p[0]*Purities_Xim[2]*Fraction_Xim[2];
    lambda_pars[4] =    Purities_p[0]*Purities_Xim[4]+Purities_p[3]*Purities_Xim[0]+Purities_p[3]*Purities_Xim[4];
    lambda_pars[3] =    1.-lambda_pars[0]-lambda_pars[1]-lambda_pars[2]-lambda_pars[4];

    //double SUM=0;
    //for(unsigned uLam=0; uLam<5; uLam++){
    //    printf("λ(pΞ)_%u = %.1f\n",uLam,lambda_pars[uLam]*100.);
    //    SUM+=lambda_pars[uLam]*100.;
    //}
    //printf("SUM: %.1f\n------------\n",SUM);
}

//        DataPeriod=="pp13TeV"?  :
//                                ;

TH2F* DLM_CommonAnaFunctions::GetResolutionMatrix(const TString& DataSample,const TString&& System){
    TString FileName;
    TString HistoName;

    if(DataSample=="pp13TeV_MB_Run2paper"){
        FileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/ResolutionMatrices/Sample6_MeV_compact.root";
    }
    else if(DataSample=="pp13TeV_HM_March19"){
        printf("\033[1;33mWARNING:\033[0m pp13TeV_HM_March19 is not available yet!\n");
        FileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/ResolutionMatrices/Sample6_MeV_compact.root";
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        FileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/ResolutionMatrices/Sample3_MeV_compact.root";
    }
    else if(DataSample=="pp13TeV_MB_BBar"){
        FileName = "/Users/sartozza/Desktop/SystematicsAndCalib/ppRun2_MB/Sample6_MeV_compact.root";
    }
    else if(DataSample=="pp13TeV_HM_BBar"){
        FileName = "/Users/sartozza/Desktop/SystematicsAndCalib/ppRun2_MB/Sample6_MeV_compact.root";
    }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        FileName = "";
    }

    if(System=="pp"){
        HistoName = "hSigmaMeV_Proton_Proton";
    }
    else if(System=="pLambda"){
        HistoName = "hSigmaMeV_Proton_Lambda";
    }
    else if(System=="LambdaLambda"){
        HistoName = "hSigmaMeV_Lambda_Lambda";
    }
    else if(System=="pXim"){
        HistoName = "hSigmaMeV_Proton_Xim";
    }
    else if(System=="ppbar"){
        HistoName = "hSigmaMeV_Proton_AntiProton";
    }
    else if(System=="pLambdabar"){
        HistoName = "hSigmaMeV_Proton_AntiLambda";
    }
    else if(System=="LambdaLambdabar"){
        HistoName = "hSigmaMeV_Lambda_AntiLambda";
    }
    else if(System=="pSigmabar"){
        HistoName = "hSigmaMeV_Proton_AntiLambda";
    }
    else if(System=="pXibar"){
        HistoName = "hSigmaMeV_Proton_AntiXim";
    }
    else{
        printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
    }
    ///FUCKING ROOT SUCKS!!!!!! SUCK MY COCK!!!! SUCK IT YOU BITCH!!!!!!!!!!
    //so we need to copy our histogram, as else we lose it when we delete the file
    //and we need to change to the "central" root directory, as else histoCopy will also be lost
    //and we need to play with the name a little bit, else we are fucked!
    TFile* FileROOT = new TFile(FileName, "read");
    TH2F* histo = (TH2F*)FileROOT->Get(HistoName);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' does not exist\n",HistoName.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH2F *histoCopy = (TH2F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}
TH2F* DLM_CommonAnaFunctions::GetResidualMatrix(const TString&& FinalSystem, const TString& InitialSystem){
    TString FileName;
    TString HistoName;

//    FileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/DecayMatrices/run2_decay_matrices_old.root";
//    FileName = "/Users/Valentina/Desktop/SystematicsAndCalib/ppRun2_MB/run2_decay_matrices_old.root";
    FileName = "/Users/sartozza/Desktop/SystematicsAndCalib/ppRun2_MB/DM_170719.root";

    if(FinalSystem=="pp"&&InitialSystem=="pLambda"){
        HistoName = "hRes_pp_pL";
    }
    else if(FinalSystem=="pLambda"&&InitialSystem=="pSigma0"){
        HistoName = "hRes_pL_pSigma0";
    }
    else if(FinalSystem=="pLambda"&&InitialSystem=="pXim"){
        HistoName = "hRes_pL_pXim";
    }
    else if(FinalSystem=="pXim"&&InitialSystem=="pXim1530"){
        HistoName = "hRes_pXim_pXim1530";
    }
    else if(FinalSystem=="ppbar"&&InitialSystem=="pLambdabar"){
        HistoName = "hRes_pp_pL";
    }
    else if(FinalSystem=="ppbar"&&InitialSystem=="LambdaLambdabar"){
        HistoName = "hRes_pp_LL";
    }
    else if(FinalSystem=="pLambdabar"&&InitialSystem=="LambdaLambdabar"){
        HistoName = "hRes_pL_LL";
    }
    else if(FinalSystem=="pLambdabar"&&InitialSystem=="pSigmabar"){
        FileName = "/Users/sartozza/Desktop/SystematicsAndCalib/ppRun2_MB/run2_decay_matrices_old.root";
        HistoName = "hRes_pL_pSigma0";
    }
    else if(FinalSystem=="pLambdabar"&&InitialSystem=="pXibar"){
        HistoName = "hRes_pL_pXim";
    }
    else{
        printf("\033[1;31mERROR:\033[0m The decay '%s->%s' does not exist\n",InitialSystem.Data(),FinalSystem.Data());
    }
    TFile* FileROOT = new TFile(FileName, "read");
    TH2F* histo = (TH2F*)FileROOT->Get(HistoName);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' does not exist\n",HistoName.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH2F *histoCopy = (TH2F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}

//iReb = 0 is 4 MeV, 1 is 8, 2 is 12, 3 is 16, 4 is 20
TH1F* DLM_CommonAnaFunctions::GetAliceExpCorrFun(const TString& DataSample,
		const TString& System,const int& iReb,const int& iSph){
    TString FileName;
    TString HistoName;


  if(DataSample=="pp13TeV_MB_BBar"){
            if(System=="ppbar"){
                FileName = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/data/MB/CFOutput_pAp_App_full.root";
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else if(System=="pLambdabar"){
                FileName = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/data/MB/CFOutput_pAL_ApL_full.root";
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else if(System=="LambdaLambdabar"){
                FileName = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/data/MB/CFOutput_LAL_ALL_full.root";
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
                // std::cout<<"Reading LLbar root file, works fine?\n"<<"HistoName = "<<HistoName<<std::endl;
            }
            else{
                printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
            }
        }else if(DataSample=="pp13TeV_HM_BBar"){
                if(System=="ppbar"){
//                    FileName = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAp_App_full.root";
                    FileName = TString::Format("/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/CF/Raw_CF/Norm018_028/CFOutput_pAp_%i.root",
                    		iSph);
                    HistoName = "hCk_ReweightedMeV_0";
                }
                else if(System=="pLambdabar"){
//                    FileName = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAL_ApL_full.root";
                    FileName = TString::Format("/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/CF/Raw_CF/Norm018_028/CFOutput_pAL_%i.root",iSph);
                    HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
                }
                else if(System=="LambdaLambdabar"){
//                    FileName = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_LAL_ALL_full.root";
                    FileName = TString::Format("/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/CF/Raw_CF/Norm018_028/CFOutput_LAL_%i.root",iSph);
                    HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
                    // std::cout<<"Reading LLbar root file, works fine?\n"<<"HistoName = "<<HistoName<<std::endl;
                }
                else{
                    printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
                }
            }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        FileName = "";
    }

    TFile* FileROOT = new TFile(FileName, "read");
    TH1F* histo = (TH1F*)FileROOT->Get(HistoName);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' does not exist\n",HistoName.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH1F *histoCopy = (TH1F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}

TH1F* DLM_CommonAnaFunctions::GetAliceExpCorrFunVar(const TString& DataSample,const TString& System,const int& iVar){
    TString FileName;
    TString HistoName;

if(DataSample=="pp13TeV_HM_BBar"){
                if(System=="ppbar"){
                  if(iVar==0){
                    FileName = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/Systematics/HMNanoAOD/SyspAp/Systematics_pAp_reb1_def.root";
                    HistoName = TString::Format("histVar_%i",iVar);
                  }
                    FileName = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/Systematics/HMNanoAOD/SyspAp/Systematics_pAp_reb1_def.root";
                    HistoName = TString::Format("histVar_%i",iVar);
                }
                else if(System=="pLambdabar"){
                    FileName = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/Systematics/HMNanoAOD/SyspAL/Systematics_pAL_rebin4.root";
                    HistoName = TString::Format("histVar_%i",iVar);
                }
                else if(System=="LambdaLambdabar"){
                    FileName = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/Systematics/HMNanoAOD/SysLAL/Systematics_LAL_reb4_def.root";
                    HistoName = TString::Format("histVar_%i",iVar);
                    // std::cout<<"Reading LLbar root file, works fine?\n"<<"HistoName = "<<HistoName<<std::endl;
                }
                else{
                    printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
                }
            }

    TFile* FileROOT = new TFile(FileName, "read");
    TDirectory* dirROOT = FileROOT->GetDirectory("Raw");
    TH1F* histo = (TH1F*)dirROOT->Get(HistoName);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' does not exist\n",HistoName.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH1F *histoCopy = (TH1F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}

TH1F* DLM_CommonAnaFunctions::GetAliceExpCorrFunmT(const TString& DataSample,
		const TString& System,const int& imT){
    TString FileName;
    TString HistoName;

if(DataSample=="pp13TeV_HM_BBar"){
                if(System=="ppbar"){
                    FileName = TString::Format("/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/mT_Analysis/Data/pAp_pair/CFOutput_mT_pApDef_HMBBar_%i.root",imT);
                    HistoName = "hCorrected_pAp";
                }
                else if(System=="pLambdabar"){
                    FileName = TString::Format("/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/mT_Analysis/Data/pAL_pair/CFOutput_mT_pALDef_HMBBar_%i.root",imT);
                    HistoName = "hCorrected_pAL";
                }
                else if(System=="LambdaLambdabar"){
                    FileName = TString::Format("/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/mT_Analysis/Data/LAL_pair/CFOutput_mT_LALDef_HMBBar_%i.root",imT);
                    HistoName = "hCorrected_LAL";
                }
                else{
                    printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
                }
            }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        FileName = "";
    }

    TFile* FileROOT = new TFile(FileName, "read");
    TH1F* histo = (TH1F*)FileROOT->Get(HistoName);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' does not exist\n",HistoName.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH1F *histoCopy = (TH1F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}

//iReb = 0 is 4 MeV, 1 is 8, 2 is 12, 3 is 16, 4 is 20
TH1F* DLM_CommonAnaFunctions::GetAliceExpCorrFunTemplate(const TString& DataSample,
		const TString& System,const int& iReb){
    TString FileName;
    TString HistoName;

if(DataSample=="pp13TeV_HM_BBar"){
                if(System=="ppbar"){
                    FileName = TString::Format("/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/CF/TemplateCorrected_CF/fOutput.root");
                    HistoName = "hCorrected_pAp";
                }
                else if(System=="pLambdabar"){
                    FileName = TString::Format("/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/CF/TemplateCorrected_CF/fOutput.root");
                    HistoName = "hCorrected_pAL";
                }
                else if(System=="LambdaLambdabar"){
                    FileName = TString::Format("/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/CF/TemplateCorrected_CF/fOutput.root");
                    HistoName = "hCorrected_LAL";
                }
                else{
                    printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
                }
            }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        FileName = "";
    }

    TFile* FileROOT = new TFile(FileName, "read");
    TH1F* histo = (TH1F*)FileROOT->Get(HistoName);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' does not exist\n",HistoName.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH1F *histoCopy = (TH1F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}

TH1F* DLM_CommonAnaFunctions::GetAliceExpCorrFunTemplatemT(const TString& DataSample,
		const TString& System,const int& imT){
    TString FileName;
    TString HistoName;

if(DataSample=="pp13TeV_HM_BBar"){
                if(System=="ppbar"){
                    FileName = TString::Format("/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/CF/TemplateCorrected_CF/mT_Corrected/fOutputpAp_%i.root",imT);
                    HistoName = "hCorrected_pAp";
                }
                else if(System=="pLambdabar"){
                    FileName = TString::Format("/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/CF/TemplateCorrected_CF/mT_Corrected/fOutputpAL_%i.root",imT);
                    HistoName = "hCorrected_pAL";
                }
                else if(System=="LambdaLambdabar"){
                    FileName = TString::Format("/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/CF/TemplateCorrected_CF/mT_Corrected/fOutputLAL_%i.root",imT);
                    HistoName = "hCorrected_LAL";
                }
                else{
                    printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
                }
            }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        FileName = "";
    }

    TFile* FileROOT = new TFile(FileName, "read");
    TH1F* histo = (TH1F*)FileROOT->Get(HistoName);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' does not exist\n",HistoName.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH1F *histoCopy = (TH1F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}

TH1F* DLM_CommonAnaFunctions::GetAliceExpCorrFunTemplatemTLocal(const TString& DataSample,
		const TString& System,const int& imT, TString InputFolder){
    TString FileName;
    TString HistoName;

if(DataSample=="pp13TeV_HM_BBar"){
                if(System=="ppbar"){
                    FileName = TString::Format(InputFolder + "/fOutputpAp_%i.root",imT);
                    HistoName = "hCorrected_pAp";
                }
                else if(System=="pLambdabar"){
                    FileName = TString::Format(InputFolder +"fOutputpAL_%i.root",imT);
                    HistoName = "hCorrected_pAL";
                }
                else if(System=="LambdaLambdabar"){
                    FileName = TString::Format(InputFolder +"fOutputLAL_%i.root",imT);
                    HistoName = "hCorrected_LAL";
                }
                else{
                    printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
                }
            }
    else{
        printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        FileName = "";
    }

    TFile* FileROOT = new TFile(FileName, "read");
    TH1F* histo = (TH1F*)FileROOT->Get(HistoName);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' does not exist\n",HistoName.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH1F *histoCopy = (TH1F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}

void DLM_CommonAnaFunctions::SetCatsFilesFolder(const TString& folder){
    CatsFilesFolder[0] = folder;
}

void DLM_CommonAnaFunctions::Clean_CommonAnaFunctions(){
    //for(unsigned uLevy=0; uLevy<NumCleverLevyObjects; uLevy++){
    //    delete CleverLevy[uLevy];
    //    CleverLevy[uLevy] = NULL;
    //}
    //delete [] CleverLevy;
    //CleverLevy=NULL;
}

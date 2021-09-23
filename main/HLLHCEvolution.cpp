#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TBox.h"
#include "TLatex.h"

#define k 8.617E-05  // Boltzmann constant
#define timeStep 60. // in minutes
#define minInYear 525600. 



double funcAlpha(double* x, double* par)
{
  double xx = x[0];
  double val = 0.;

  double p0=39.17; double p1=-3.243; double p2=-0.1545; double p3=0.01374; 
  for(int i = 0; i < 20; ++i)
  {
    double tau = pow(10,i-5);
    double alpha0 = (p0+p1*(i-5)+p2*pow(i-5,2)+p3*pow(i-5,3))/20.;
    val += alpha0*exp(-xx/tau);
  }
  
  return val;
}



double timeScale(double T_a, double T_r)
{
  T_a = T_a + 273.15;
  T_r = T_r + 273.15;
  return exp( (-1./k) * ((1.1692-0.00049*T_a*T_a/(T_a+655.))/T_a-(1.1692-0.00049*T_r*T_r/(T_r+655.))/T_r));

}



double myfunc(double* x, double* par)
{
  double xx = x[0];
  
  double x0 = par[4];

  double A = par[0];
  double C = par[1];
  double D = par[2];
  double E = par[3];
  double B = D*E*exp(E*x0) - 2.*C*x0;
  double F = A + B*x0 + C*x0*x0 - D*exp(E*x0);

  if(xx < x0) return A + B*xx + C*xx*xx;
  else        return D*exp(E*xx)+F;
}






int main(int argc, char** argv)
{
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  std::string SiPMType = opts.GetOpt<std::string>("Input.SiPMType");
    
  double T_ref = 21.;
  double T_ann = opts.GetOpt<double>("Input.Tann");
  double T_CO2 = opts.GetOpt<double>("Input.TCO2");

  int useTECs = opts.GetOpt<int>("Input.useTECs");
  double TECsDeltaT = opts.GetOpt<double>("Input.TECsDeltaT");
  
  double T_op = T_CO2;
  if( useTECs ) T_op -= TECsDeltaT;

  int interfillAnnealing = opts.GetOpt<int>("Input.interfillAnnealing");
  double interfillTemp = opts.GetOpt<double>("Input.interfillTemp");
  
  double maxPower = opts.GetOpt<double>("Input.maxPower");
  
  int year, month, irr;
  double lum, lumi;
  double alpha = 0., tau, alpha0;
  double  norm = 0.;

  double vett[20];
  for(int j = 0; j < 20; ++j) vett[j] = 0;


  std::string HLLHCScheduleLabel = opts.GetOpt<std::string>("Input.HLLHCSchedule");
  std::ifstream HLLHCSchedule;
  HLLHCSchedule.open(Form("data/%s.txt",HLLHCScheduleLabel.c_str()), std::ios::in);
  TGraph g_instLumi_vs_time, g_alpha_vs_time, g_alphaNorm_vs_time, g_instLumidark_vs_time, g_dcr_vs_time;
  TGraph g_fluence_vs_time;
  TGraph g_intLumi_vs_time;
  TGraph g_temp_vs_time;
  
  
  
  
  //----------------------
  // parse HL-LHC schedule
  std::string line;
  while(true)
  {
    getline(HLLHCSchedule,line);
    if( !HLLHCSchedule.good() ) break;

    if( line[0] == '#' ) continue;

    std::stringstream ss(line);
    ss >> year >> month >> irr >> lum;
    
    
    double t = (year*12+month) * 43200; // in minutes
    double theta = 1.;
    double temp;
    
    // define luminosity for each period
    if( irr == 0 ) // technical stops - assume room temperature (15°C)
    {
      lumi = 0.;
      theta = timeScale(15., T_ref);
      temp = T_ann;
    }
    
    else if( irr == 4 ) // technical stops - assume annealing temperature T_ann
    {
      lumi = 0.;
      theta = timeScale(T_ann, T_ref);
      temp = T_ann;
    }
    
    else if( irr == 1 ) // beam commissioning - temperature is T_op
    {
      lumi = lum*0.1;
      theta = timeScale(T_op, T_ref);
      temp = T_op;
    }
    
    else if( irr == 3 ) // Pb-Pb collisions - temperature is T_op
    {
      lumi = lum*0.001;
      theta=timeScale(T_op, T_ref);
      temp = T_op;
    }
    
    else if( irr == 2 || irr == 9 || irr == 91 ) // pp collisions - temperature is T_op
    {
      lumi = lum;
      theta = timeScale(T_op, T_ref);
      temp = T_op;
    }
    
    

    //------------
    // alpha graph
    double p0 = 39.17;
    double p1 = -3.243;
    double p2 = -0.1545;
    double p3 = 0.01374; 
    
    for(int j = 0; j < 30*24; ++j) // loop over 30*24 hours in a month 
    { 
      float effLumi = lumi;
      float effTemp = temp;

      if( interfillAnnealing && irr == 2 )
      {
        if(
          (j >= 110 && j < 122) ||
          (j >= 232 && j < 244) ||
          (j >= 354 && j < 366) ||
          (j >= 476 && j < 488) ||
          (j >= 598 && j < 610) &&
          useTECs
          )
        {
          effLumi = 0;
          theta = timeScale(interfillTemp, T_ref);
          effTemp = interfillTemp;
        }
        else
        {
          effLumi = lumi;
          theta = timeScale(T_op, T_ref);
          effTemp = T_op;
        }
      }
      
      // // only 12 days in a month of collisions. 
      // if( lumi > 0 )
      // {
      //   if( j%24 < 12 ) effLumi = lumi * 0.919;
      //   else            effLumi = 0;
      // }

      if( irr == 5 || irr == 51 ) // assume 2 weeks of annealing at 0° C during the data taking
      {
        if( j < 288 ) effLumi = lumi * 0.919;
        else          effLumi = 0;
        
        if( j >=360 )
        {
          effLumi = 0.;
          theta = timeScale(0., T_ref);
          effTemp = 0.;
        }
      }
      if( irr == 6 || irr == 61 ) // assume 2 weeks of annealing at 15° C during the data taking
      {
        if( j < 288 ) effLumi = lumi * 0.919;
        else          effLumi = 0;
        
        if( j >= 360 )
        {
          effLumi = 0.;
          theta = timeScale(15., T_ref);
          effTemp = 15.;
        }
      }
      if( irr == 7 || irr == 71 ) // assume 2 weeks of annealing at 40° C during the data taking
      {
        if( j < 288 ) effLumi = lumi * 0.919;
        else          effLumi = 0;
        
        if( j >= 360 )
        {
          effLumi = 0.;
          theta = timeScale(40., T_ref);
          effTemp = 40.;
        }
      }
      if( irr == 8 || irr == 81 ) // assume 2 weeks of annealing at 100° C during the data taking
      {
        if( j < 288 ) effLumi = lumi * 0.919;
        else          effLumi = 0;
        
        if( j >= 360 )
        {
          effLumi = 0.;
          theta = timeScale(70., T_ref);
          effTemp = 70.;
        }
      }
      if( irr == 9 || irr == 91 ) // assume 4 days every ~6 weeks of annealing at 40° C (15° C without TECs) during the data taking
      {
        if( !interfillAnnealing && j < 624 )
        {
          effLumi = lumi;
          theta = timeScale(T_op, T_ref);
          effTemp = T_op;
        }
        
        if( interfillAnnealing && j < 624 )
        {
          if(
            (j >= 110 && j < 122) ||
            (j >= 232 && j < 244) ||
            (j >= 354 && j < 366) ||
            (j >= 476 && j < 488) ||
            (j >= 598 && j < 610) &&
            useTECs
            )
          {
            effLumi = 0;
            theta = timeScale(interfillTemp, T_ref);
            effTemp = interfillTemp;
          }
          else
          {
            effLumi = lumi;
            theta = timeScale(T_op, T_ref);
            effTemp = T_op;
          }          
        }
        
        if( j >= 624 )
        {
          effLumi = 0.;
          if( useTECs )
          {
            if( irr == 9 )
            {
              theta = timeScale(interfillTemp, T_ref);
              effTemp = interfillTemp;
            }
            if( irr == 91 )
            {
              theta = timeScale(40., T_ref);
              effTemp = 40.;
            }
          }
          else
          {
            if( irr == 9 )
            {
              theta = timeScale(T_op, T_ref);
              effTemp = T_op;
            }
            if( irr == 91 )
            {
              theta = timeScale(0., T_ref);
              effTemp = 0.;
            }
          }
        }
      }
      
      g_instLumi_vs_time.SetPoint(g_instLumi_vs_time.GetN(), (j*timeStep+t)/minInYear, effLumi);
      g_temp_vs_time.SetPoint(g_temp_vs_time.GetN(),(j*timeStep+t)/minInYear,effTemp);
      
      norm += (effLumi*timeStep);
      
      for(int i = 0; i < 20; ++i)
      {
        tau = pow(10,i-5);
        alpha0 = (p0+p1*(i-5)+p2*pow(i-5,2)+p3*pow(i-5,3)) / 20.;
        
        vett[i] = alpha0*tau/theta*effLumi*(1-exp(-theta*timeStep/tau)) + (vett[i])*exp(-theta*timeStep/tau);
        alpha+=vett[i];
      }
      
      g_alpha_vs_time.SetPoint(g_alpha_vs_time.GetN(), (j*timeStep+t)/minInYear, alpha);
      g_fluence_vs_time.SetPoint(g_fluence_vs_time.GetN(), (j*timeStep+t)/minInYear, norm);
      g_intLumi_vs_time.SetPoint(g_intLumi_vs_time.GetN(), (j*timeStep+t)/minInYear, norm);
      
      alpha = 0.; 
    }
  }
  
  
  // normalize to total fluence
  for(int point = 0; point < g_alpha_vs_time.GetN(); ++point)
  {
    double x,y;
    g_alpha_vs_time.GetPoint(point,x,y);
    g_alphaNorm_vs_time.SetPoint(point,x,y/norm);
    
    g_fluence_vs_time.GetPoint(point,x,y);
    g_fluence_vs_time.SetPoint(point,x,y/norm*1.9E14);

    g_intLumi_vs_time.GetPoint(point,x,y);
    g_intLumi_vs_time.SetPoint(point,x,y/norm*3000);
  }
  
  
  
  
  //------------------------------
  // compute optimum working point
  std::string TECsLabel = useTECs ? "_TECs" : "";
  std::string interfillLabel = interfillAnnealing ? "_interfillAnnealing" : "";
  TFile* outFile = TFile::Open(Form("plots/outFile_%s_Top%d_Tann%d%s%s_%s_maxPower%.0f.root",SiPMType.c_str(),int(T_op),int(T_ann),TECsLabel.c_str(),interfillLabel.c_str(),HLLHCScheduleLabel.c_str(),maxPower),"RECREATE");
  
  TGraph g_tResBest_vs_time;
  TGraph g_VbiasBest_vs_time;
  TGraph g_VovBest_vs_time;
  TGraph g_DCRBest_vs_time;
  TGraph g_PDEBest_vs_time;
  TGraph g_nPEBest_vs_time;
  TGraph g_gainBest_vs_time;
  TGraph g_SoNBest_vs_time;
  TGraph g_powerBest_vs_time;
  TGraph g_staticPowerBest_vs_time;
  TGraph g_dynamicPowerBest_vs_time;
  TGraph g_TECsPowerBest_vs_time;
  TGraph g_currentBest_vs_time;
  TGraph g_staticCurrentBest_vs_time;
  TGraph g_dynamicCurrentBest_vs_time;
  
  TGraph g_tResBest_vs_intLumi;
  TGraph g_VbiasBest_vs_intLumi;
  TGraph g_VovBest_vs_intLumi;
  TGraph g_DCRBest_vs_intLumi;
  TGraph g_PDEBest_vs_intLumi;
  TGraph g_nPEBest_vs_intLumi;
  TGraph g_gainBest_vs_intLumi;
  TGraph g_SoNBest_vs_intLumi;
  TGraph g_powerBest_vs_intLumi;
  TGraph g_staticPowerBest_vs_intLumi;
  TGraph g_dynamicPowerBest_vs_intLumi;
  TGraph g_TECsPowerBest_vs_intLumi;
  TGraph g_currentBest_vs_intLumi;
  TGraph g_staticCurrentBest_vs_intLumi;
  TGraph g_dynamicCurrentBest_vs_intLumi;
  
  TF1* f_PDE = new TF1("func","[0]*(1-exp(-1.*[1]*x))",0.,10.);
  if( SiPMType == "HPK_PDELoss07p5" || SiPMType == "HPK_PDELoss15p0" ) f_PDE -> SetParameters(0.394321,0.738063);
  if( SiPMType == "FBK" ) f_PDE -> SetParameters(0.454306,0.363804);
  TFile* sipmParams;
  if( SiPMType == "HPK_PDELoss07p5" || SiPMType == "HPK_PDELoss15p0" ) sipmParams = TFile::Open("data/sipm_spec_input_HDR2-015-v2-1e13.root","READ");
  if( SiPMType == "FBK" ) sipmParams = TFile::Open("data/sipm_spec_input_FBK-W7C-1e13.root","READ");
  TF1* f_gain = (TF1*)( sipmParams->Get("fGain_vs_OV") );
  TF1* f_ENF = (TF1*)( sipmParams->Get("fENF_vs_OV") );

  float LCE = 0.15;
  if( SiPMType == "FBK" ) LCE *= 1.17;
  
  TFile* inFile_DCRParams = TFile::Open("data/DCRParams_new.root","READ");
  TGraph* g_DCR_par0 = (TGraph*)( inFile_DCRParams->Get("g_par0"));
  TGraph* g_DCR_par1 = (TGraph*)( inFile_DCRParams->Get("g_par1"));
  TGraph* g_DCR_par2 = (TGraph*)( inFile_DCRParams->Get("g_par2"));
  TGraph* g_DCR_par3 = (TGraph*)( inFile_DCRParams->Get("g_par3"));
  TGraph* g_DCR_par4 = (TGraph*)( inFile_DCRParams->Get("g_par4"));

  TFile* inFile_TOFHIRPerformance = TFile::Open("data/TOFHIR2BPerformance.root");
  std::vector<int> DCRVals;
  DCRVals.push_back(10);
  DCRVals.push_back(25);
  DCRVals.push_back(40);
  DCRVals.push_back(60);
  DCRVals.push_back(80);
  std::map<int,TF1*> tRes_tot_vs_nPE;
  for(auto DCRVal : DCRVals)
  {
    tRes_tot_vs_nPE[DCRVal] = (TF1*)( inFile_TOFHIRPerformance->Get(Form("tRes_tot_vs_nPE_DCR%dGHz",DCRVal)) );
  }
  
  TFile* inFile_TOFHIRPerformance_gainCorrection = TFile::Open("data/TOFHIR2BPerformance_gainCorrection.root");
  TGraph* g_TOFHIRPerformance_gainCorrection = (TGraph*)( inFile_TOFHIRPerformance_gainCorrection->Get("gainCorrection"));

  TFile* inFile_TECsPower = TFile::Open("data/TECsPower.root");
  TF1* f_TECsPower_noSiPMLoad = (TF1*)( inFile_TECsPower->Get("f_noSiPMLoad") );
  
  // TF1* turnOn_PDE = new TF1("turnOn_PDE","0.15/(1+exp(-1.*(x-1E14)/1E13))",0.,2E14);
  // TF1* turnOn_gain = new TF1("turnOn_gain","0.30/(1+exp(-1.*(x-1E14)/1E13))",0.,2E14);
   // TF1* turnOn_PDE = new TF1("turnOn_PDE","1.-0.15/2E14*x",0.,2E14);
   // TF1* turnOn_gain = new TF1("turnOn_gain","1-0.30/2E14*x",0.,2E14);
  TF1* turnOn_PDE;
  TF1* turnOn_gain;
  if( SiPMType == "HPK_PDELoss07p5" )
  {
    turnOn_PDE = new TF1("turnOn_PDE","1.-0.075/2E14*x",0.,2E14);
    turnOn_gain = new TF1("turnOn_gain","1-0.30/2E14*x",0.,2E14);
  }
  if( SiPMType == "HPK_PDELoss15p0" )
  {
    turnOn_PDE = new TF1("turnOn_PDE","1.-0.150/2E14*x",0.,2E14);
    turnOn_gain = new TF1("turnOn_gain","1-0.30/2E14*x",0.,2E14);
  }
  if( SiPMType == "FBK" )
  {
    turnOn_PDE = new TF1("turnOn_PDE","1.-0./2E14*x",0.,2E14);
    turnOn_gain = new TF1("turnOn_gain","1-0.25/2E14*x",0.,2E14);
  }
  
  float tResAvg = 0.;
  int nTResAvg = 0;
  
  for(int point = 0; point < g_fluence_vs_time.GetN(); ++point)
  {
    double time, fluence, alpha, instLumi, intLumi;
    g_instLumi_vs_time.GetPoint(point,time,instLumi);
    g_fluence_vs_time.GetPoint(point,time,fluence);
    g_intLumi_vs_time.GetPoint(point,time,intLumi);
    g_alphaNorm_vs_time.GetPoint(point,time,alpha);
    
    if( instLumi == 0. ) continue;
    
    
    //------------------
    // evaluate SiPM Vbr
    float intercept = 0.;
    float slope = 0.;

    if( SiPMType == "HPK_PDELoss07p5" || SiPMType == "HPK_PDELoss15p0" )
    { 
      intercept = 37.95 + 1.19E-14*fluence;   // at 0° C, from Carlos data
      slope = 37.5 + 1.1E-14*fluence;   // mV/° C, from Carlos data
    }
    if( SiPMType == "FBK" )
    {
      intercept = 32.04 + 0.28E-14*fluence;   // at 0° C, from SiPM specs document
      slope = 31.2 + 1.1E-14*fluence;   // mV/° C, from SiPM specs document
    }
    
    float Vbr = intercept + 0.001*slope*(T_op-0.);
    
    float DCRRef_par0 = g_DCR_par0->Eval(fluence);
    float DCRRef_par1 = g_DCR_par1->Eval(fluence);
    float DCRRef_par2 = g_DCR_par2->Eval(fluence);
    float DCRRef_par3 = g_DCR_par3->Eval(fluence);
    float DCRRef_par4 = g_DCR_par4->Eval(fluence);
    // TF1* f_DCRRef_vs_Vov = new TF1("f_DCRRef_vs_Vov","pol2",0.,7.);
    // f_DCRRef_vs_Vov -> SetParameters(DCRRef_par0,DCRRef_par1,DCRRef_par2);
    TF1* f_DCRRef_vs_Vov = new TF1("f_DCRRef_vs_Vov",myfunc,0.,7.,5);
    f_DCRRef_vs_Vov -> SetParameters(DCRRef_par0,DCRRef_par1,DCRRef_par2,DCRRef_par3,DCRRef_par4);
      
    float DCRRef = f_DCRRef_vs_Vov -> Eval(1.); // DCR at 1 V and - 30° C as per Carlos data
    
    float tResBest = 999999.;
    float VbiasBest = -1.;
    float VovBest = -1.;
    float DCRBest = -1.;
    float PDEBest = 0.;
    float nPEBest = -1.;
    float gainBest = -1.;
    float SoNBest = -1.;
    float powerBest = 0.;
    float dynamicPowerBest = 0.;
    float staticPowerBest = 0.;
    float TECsPowerBest = 0.;
    float currentBest = 0.;
    float dynamicCurrentBest = 0.;
    float staticCurrentBest = 0.;
    
    for(float Vov = 0.2; Vov < 5.; Vov += 0.01)
    {
      //-----------------------------------------------
      // evaluate effective DCR for a given OV and T_op
      float DCR = ( 9280.19 * alpha * 1.E-17 * 1.9E14 ) * // from Sasha, at -35° C and 1 V
        (f_DCRRef_vs_Vov->Eval(Vov)/DCRRef)               // morphing vs. OV using Carlos data
        ;

      if( SiPMType == "FBK" ) DCR *= 1.10;
      
      float B = -0.00416498*alpha + 0.0798623;   // DCR scaling with temperature, including dependence of scaling factor on alpha
      
      DCR = DCR * exp(B*(T_op+35.));
      
      
      //-----------------------------------------
      // evaluate total power dissipated per SiPM
      float staticCurrent = DCR*1E09 * f_ENF->Eval(Vov) * f_gain->Eval(Vov)*(turnOn_gain->Eval(fluence)) * 1.602E-19;
      float staticPower = staticCurrent * (Vbr + Vov) * 1000.;   // in mW
      if( staticPower > maxPower ) continue;
      
      float dynamicCurrent = (2.3*1E06 * instLumi * 5/7.5) * (1.0 * 4.2 * 40000 * LCE * f_PDE->Eval(Vov)*(turnOn_PDE->Eval(fluence))) * f_ENF->Eval(Vov) * f_gain->Eval(Vov)*(turnOn_gain->Eval(fluence)) * 1.602E-19;   // 2.3 MHz equivalent MIP rate at 200 PU.
      float dynamicPower = dynamicCurrent * (Vbr + Vov) * 1000.;   // in mW
      if( (staticPower+dynamicPower) > maxPower ) continue;

      float TECsPower = 0.; // in mW per channel
      if( useTECs )
      {
        f_TECsPower_noSiPMLoad -> SetParameter(0.,(3./430.)*16.*(staticPower+dynamicPower));
        TECsPower = f_TECsPower_noSiPMLoad -> Eval(-1.*TECsDeltaT) / 16.;
      }
      if( (staticPower+dynamicPower+TECsPower) > maxPower ) continue;
      
      
      //-------------------------
      // evaluate time resolution
      float nPE = 1.0 * 4.2 * 40000 * LCE * f_PDE->Eval(Vov)*(turnOn_PDE->Eval(fluence));
      
      float tResCurr = sqrt( 
        pow(22.3*sqrt(9500./nPE),2) +
        pow(sqrt(pow(47.3,2)-pow(22.3*sqrt(9500./7000.),2)) * 7000/nPE*pow(DCR/45.,0.28),2)
        ) * 
        g_TOFHIRPerformance_gainCorrection->Eval(f_gain->Eval(Vov)*(turnOn_gain->Eval(fluence)));

      tResCurr = sqrt(tResCurr*tResCurr + 15*15); // add clock contribution in quadrature
      
      // TGraph* temp = new TGraph();
      // for(auto DCRVal : DCRVals)
      // {
      //   temp -> SetPoint(temp->GetN(),DCRVal,tRes_tot_vs_nPE[DCRVal]->Eval(nPE));
      // }
      // float tResCurr = temp->Eval(DCR);
      // delete temp;
      
      if( tResCurr < tResBest )
      {
        tResBest = tResCurr;
        VbiasBest = Vov + Vbr;
        VovBest = Vov;
        DCRBest = DCR;
        PDEBest = f_PDE->Eval(Vov)*(turnOn_PDE->Eval(fluence));
        nPEBest = nPE;
        gainBest = f_gain->Eval(Vov)*(turnOn_gain->Eval(fluence));
        SoNBest = f_PDE->Eval(Vov)*(turnOn_PDE->Eval(fluence))/sqrt(DCR);
        powerBest = staticPower+dynamicPower+TECsPower;
        dynamicPowerBest = dynamicPower;
        staticPowerBest = staticPower;
        TECsPowerBest = TECsPower;
        currentBest = (staticPower+dynamicPower) / (Vbr + Vov);
        dynamicCurrentBest = (dynamicPower) / (Vbr + Vov);
        staticCurrentBest = (staticPower) / (Vbr + Vov);
      }
    }

    std::cout << "Vov: " << VovBest << "   nPEBest: " << nPEBest << "   gainBest: " << gainBest << "   DCRBest: " << DCRBest << "   tResBest: " << tResBest << std::endl;
    g_tResBest_vs_time.SetPoint(g_tResBest_vs_time.GetN(),time,tResBest);
    g_VbiasBest_vs_time.SetPoint(g_VbiasBest_vs_time.GetN(),time,VbiasBest);
    g_VovBest_vs_time.SetPoint(g_VovBest_vs_time.GetN(),time,VovBest);
    g_DCRBest_vs_time.SetPoint(g_DCRBest_vs_time.GetN(),time,DCRBest);
    g_PDEBest_vs_time.SetPoint(g_PDEBest_vs_time.GetN(),time,PDEBest);
    g_nPEBest_vs_time.SetPoint(g_nPEBest_vs_time.GetN(),time,nPEBest);
    g_gainBest_vs_time.SetPoint(g_gainBest_vs_time.GetN(),time,gainBest);
    g_SoNBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,SoNBest);
    g_powerBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,powerBest);
    g_dynamicPowerBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,dynamicPowerBest);
    g_staticPowerBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,staticPowerBest);
    g_TECsPowerBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,TECsPowerBest);
    g_currentBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,currentBest);
    g_dynamicCurrentBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,dynamicCurrentBest);
    g_staticCurrentBest_vs_time.SetPoint(g_SoNBest_vs_time.GetN(),time,staticCurrentBest);
    
    g_tResBest_vs_intLumi.SetPoint(g_tResBest_vs_intLumi.GetN(),intLumi,tResBest);
    g_VbiasBest_vs_intLumi.SetPoint(g_VbiasBest_vs_intLumi.GetN(),intLumi,VbiasBest);
    g_VovBest_vs_intLumi.SetPoint(g_VovBest_vs_intLumi.GetN(),intLumi,VovBest);
    g_DCRBest_vs_intLumi.SetPoint(g_DCRBest_vs_intLumi.GetN(),intLumi,DCRBest);
    g_PDEBest_vs_intLumi.SetPoint(g_PDEBest_vs_intLumi.GetN(),intLumi,PDEBest);
    g_nPEBest_vs_intLumi.SetPoint(g_nPEBest_vs_intLumi.GetN(),intLumi,nPEBest);
    g_gainBest_vs_intLumi.SetPoint(g_gainBest_vs_intLumi.GetN(),intLumi,gainBest);
    g_SoNBest_vs_intLumi.SetPoint(g_SoNBest_vs_intLumi.GetN(),intLumi,SoNBest);
    g_powerBest_vs_intLumi.SetPoint(g_powerBest_vs_intLumi.GetN(),intLumi,powerBest);
    g_dynamicPowerBest_vs_intLumi.SetPoint(g_dynamicPowerBest_vs_intLumi.GetN(),intLumi,dynamicPowerBest);
    g_staticPowerBest_vs_intLumi.SetPoint(g_staticPowerBest_vs_intLumi.GetN(),intLumi,staticPowerBest);
    g_TECsPowerBest_vs_intLumi.SetPoint(g_TECsPowerBest_vs_intLumi.GetN(),intLumi,TECsPowerBest);
    g_currentBest_vs_intLumi.SetPoint(g_currentBest_vs_intLumi.GetN(),intLumi,currentBest);
    g_dynamicCurrentBest_vs_intLumi.SetPoint(g_dynamicCurrentBest_vs_intLumi.GetN(),intLumi,dynamicCurrentBest);
    g_staticCurrentBest_vs_intLumi.SetPoint(g_staticCurrentBest_vs_intLumi.GetN(),intLumi,staticCurrentBest);
    
    tResAvg += tResBest;
    ++nTResAvg;
  }
  
  tResAvg /= nTResAvg;
  
  


  //------------
  // save graphs
  outFile -> cd();
  
  g_alphaNorm_vs_time.GetXaxis()->SetTitle("years from 2027");
  g_alphaNorm_vs_time.GetYaxis()->SetTitle("#alpha [10^{-17} A/cm]");
  g_alphaNorm_vs_time.SetMarkerStyle(20);
  g_alphaNorm_vs_time.SetMarkerSize(0.5);
  g_alphaNorm_vs_time.Write("g_alphaNorm_vs_time");

  g_fluence_vs_time.GetXaxis()->SetTitle("years from 2027");
  g_fluence_vs_time.GetYaxis()->SetTitle("fluence [cm^{-2}]");
  g_fluence_vs_time.SetMarkerStyle(20);
  g_fluence_vs_time.SetMarkerSize(0.5);
  g_fluence_vs_time.Write("g_fluence_vs_time");
  
  g_intLumi_vs_time.GetXaxis()->SetTitle("years from 2027");
  g_intLumi_vs_time.GetYaxis()->SetTitle("int. luminosity [fb^{-1}]");
  g_intLumi_vs_time.SetMarkerStyle(20);
  g_intLumi_vs_time.SetMarkerSize(0.5);
  g_intLumi_vs_time.Write("g_intLumi_vs_time");
  
  g_instLumi_vs_time.GetXaxis()->SetTitle("years from 2027");
  g_instLumi_vs_time.GetYaxis()->SetTitle("inst. luminosity [a.u.]");
  g_instLumi_vs_time.SetMarkerStyle(20);
  g_instLumi_vs_time.SetMarkerSize(0.5);
  g_instLumi_vs_time.Write("g_instLumi_vs_time");
  
  g_temp_vs_time.GetXaxis()->SetTitle("years from 2027");
  g_temp_vs_time.GetYaxis()->SetTitle("temperature [#circ C]");
  g_temp_vs_time.SetMarkerStyle(20);
  g_temp_vs_time.SetMarkerSize(0.5);
  g_temp_vs_time.Write("g_temp_vs_time");
  
  g_tResBest_vs_time.SetTitle(";years from 2027;#sigma_{t}^{best} [ps]");
  g_tResBest_vs_time.Write("g_tResBest_vs_time");
  g_VbiasBest_vs_time.SetTitle(";years from 2027;V_{bias} at best #sigma_{t} [V]");
  g_VbiasBest_vs_time.Write("g_VbiasBest_vs_time");
  g_VovBest_vs_time.SetTitle(";years from 2027;V_{OV} at best #sigma_{t} [V]");
  g_VovBest_vs_time.Write("g_VovBest_vs_time");
  g_DCRBest_vs_time.SetTitle(";years from 2027;DCR at best #sigma_{t} [GHz]");
  g_DCRBest_vs_time.Write("g_DCRBest_vs_time");
  g_PDEBest_vs_time.SetTitle(";years from 2027;PDE at best #sigma_{t}");
  g_PDEBest_vs_time.Write("g_PDEBest_vs_time");
  g_nPEBest_vs_time.SetTitle(";years from 2027;N_{p.e.} at best #sigma_{t}");
  g_nPEBest_vs_time.Write("g_nPEBest_vs_time");
  g_gainBest_vs_time.SetTitle(";years from 2027;SiPM gain at best #sigma_{t}");
  g_gainBest_vs_time.Write("g_gainBest_vs_time");
  g_SoNBest_vs_time.SetTitle(";years from 2027;S/N at best #sigma_{t}");
  g_SoNBest_vs_time.Write("g_SoNBest_vs_time");
  g_powerBest_vs_time.SetTitle(";years from 2027;total power per ch. at best #sigma_{t} [mW]");
  g_powerBest_vs_time.Write("g_powerBest_vs_time");
  g_dynamicPowerBest_vs_time.SetTitle(";years from 2027;dynamic power per ch. at best #sigma_{t} [mW]");
  g_dynamicPowerBest_vs_time.Write("g_dynamicPowerBest_vs_time");
  g_staticPowerBest_vs_time.SetTitle(";years from 2027;static power per ch. at best #sigma_{t} [mW]");
  g_staticPowerBest_vs_time.Write("g_staticPowerBest_vs_time");
  g_TECsPowerBest_vs_time.SetTitle(";years from 2027;TECs power per ch. at best #sigma_{t} [mW]");
  g_TECsPowerBest_vs_time.Write("g_TECsPowerBest_vs_time");
  g_currentBest_vs_time.SetTitle(";years from 2027;total current per ch. at best #sigma_{t} [mA]");
  g_currentBest_vs_time.Write("g_currentBest_vs_time");
  g_dynamicCurrentBest_vs_time.SetTitle(";years from 2027;dynamic current per ch. at best #sigma_{t} [mA]");
  g_dynamicCurrentBest_vs_time.Write("g_dynamicCurrentBest_vs_time");
  g_staticCurrentBest_vs_time.SetTitle(";years from 2027;static current per ch. at best #sigma_{t} [mA]");
  g_staticCurrentBest_vs_time.Write("g_staticCurrentBest_vs_time");
  
  g_tResBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];#sigma_{t}^{best} [ps]");
  g_tResBest_vs_intLumi.Write("g_tResBest_vs_intLumi");
  g_VbiasBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];V_{OV} at best #sigma_{t} [V]");
  g_VbiasBest_vs_intLumi.Write("g_VbiasBest_vs_intLumi");
  g_VovBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];V_{OV} at best #sigma_{t} [V]");
  g_VovBest_vs_intLumi.Write("g_VovBest_vs_intLumi");
  g_DCRBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];DCR at best #sigma_{t} [GHz]");
  g_DCRBest_vs_intLumi.Write("g_DCRBest_vs_intLumi");
  g_PDEBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];PDE at best #sigma_{t}");
  g_PDEBest_vs_intLumi.Write("g_PDEBest_vs_intLumi");
  g_nPEBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];N_{p.e.} at best #sigma_{t}");
  g_nPEBest_vs_intLumi.Write("g_nPEBest_vs_intLumi");
  g_gainBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];SiPM gain at best #sigma_{t}");
  g_gainBest_vs_intLumi.Write("g_gainBest_vs_intLumi");
  g_SoNBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];S/N at best #sigma_{t}");
  g_SoNBest_vs_intLumi.Write("g_SoNBest_vs_intLumi");
  g_powerBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];total power per ch. at best #sigma_{t} [mW]");
  g_powerBest_vs_intLumi.Write("g_powerBest_vs_intLumi");
  g_dynamicPowerBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];dynamic power per ch. at best #sigma_{t} [mW]");
  g_dynamicPowerBest_vs_intLumi.Write("g_dynamicPowerBest_vs_intLumi");
  g_staticPowerBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];static power per ch. at best #sigma_{t} [mW]");
  g_staticPowerBest_vs_intLumi.Write("g_staticPowerBest_vs_intLumi");
  g_TECsPowerBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];TECs power per ch. at best #sigma_{t} [mW]");
  g_TECsPowerBest_vs_intLumi.Write("g_TECsPowerBest_vs_intLumi");
  g_currentBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];total current per ch. at best #sigma_{t} [mA]");
  g_currentBest_vs_intLumi.Write("g_currentBest_vs_intLumi");
  g_dynamicCurrentBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];dynamic current per ch. at best #sigma_{t} [mA]");
  g_dynamicCurrentBest_vs_intLumi.Write("g_dynamicCurrentBest_vs_intLumi");
  g_staticCurrentBest_vs_intLumi.SetTitle(";int. luminosity [fb^{-1}];static current per ch. at best #sigma_{t} [mA]");
  g_staticCurrentBest_vs_intLumi.Write("g_staticCurrentBest_vs_intLumi");
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
  return 0;
}

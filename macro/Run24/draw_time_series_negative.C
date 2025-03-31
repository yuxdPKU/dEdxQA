#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <sPhenixStyle.C>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

void draw_time_series_negative(std::string infile = "run_ratio_negative.txt") {

  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  std::ifstream inputFile("/phenix/u/jinhuang/links/sPHENIX_work/forXudong/ZDC_run_export_from_1711944000.0_to_1728964740.0.csv");
    
  if (!inputFile.is_open()) {
    std::cerr << "can not open file!" << std::endl;
  }

  std::map<int,float> runtimestamp;
  std::string line;

  std::getline(inputFile, line);

  while (std::getline(inputFile, line)) {
    std::stringstream ss(line);
    std::string token;
    int runnumber;
    float brtime;

    std::getline(ss, token, ',');  // skip
    std::getline(ss, token, ',');  // get runnumber
    runnumber = std::stoi(token);
    std::getline(ss, token, ',');  // get brtime
    brtime = std::stod(token);

    runtimestamp[runnumber] = brtime;
  }

  inputFile.close();

  TCanvas *c1 = new TCanvas("c1", "TGraph Example", 800, 600);
  
  std::ifstream inputFile2(infile);

  if (!inputFile2) {
      std::cerr << "Can not open file!" << std::endl;
      return;
  }

  std::vector<float> vec_runNo, vec_runtime, vec_ratio;
  float runNo, runtime, ratio;
  int runmin=1000000, runmax=0;

  while (inputFile2 >> runNo >> ratio)
  {
    if (std::isnan(runNo) || std::isnan(ratio))
    {
      continue;
    }
    //if (ratio>2) continue;
    auto it = runtimestamp.find(runNo);
    if (it != runtimestamp.end())
    {
      if (runNo < runmin) {runmin = runNo;}
      if (runNo > runmax) {runmax = runNo;}
      vec_runtime.push_back(it->second);
      vec_runNo.push_back(runNo);
      vec_ratio.push_back(ratio);
    }
  }

  std::cout<<"valid run number = "<<vec_runNo.size()<<std::endl;

  TGraph *graph = new TGraph(vec_runNo.size(), vec_runtime.data(), vec_ratio.data());
  
  graph->SetTitle("dE/dx time series in pp Run;Date;dE/dx ratio");
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.8);
  graph->SetMarkerColor(kBlue);
  graph->SetLineColor(kRed);
  graph->SetLineWidth(2);
  
  graph->Draw("ALP");

  //TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
  //legend->AddEntry(graph, "y = x^2", "lp");
  //legend->Draw();

  auto xminIt = std::min_element(vec_runtime.begin(), vec_runtime.end());
  auto xmaxIt = std::max_element(vec_runtime.begin(), vec_runtime.end());
  float xmin = *xminIt;
  float xmax = *xmaxIt;

  TLine *line1 = new TLine(xmin-200000, 1, xmax+200000, 1);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->SetLineColor(kRed);
  line1->Draw();

  TAxis *xAxis = graph->GetXaxis();
  xAxis->SetNdivisions(5);
  xAxis->SetTimeDisplay(1);
  xAxis->SetTimeFormat("%m-%d");
  xAxis->SetTimeOffset(0, "gmt");

  TPaveText *pt = new TPaveText(0.48, 0.55, 0.9, 0.9, "brNDC");
  pt->AddText("#it{#bf{sPHENIX}} Internal");
  pt->AddText("p+p #sqrt{s}=200 GeV");
  pt->AddText(Form("Run %d - %d",runmin,runmax));
  pt->AddText(Form("Laudau Fit to dE/dx"));
  pt->AddText(Form("Mean HighZ / LowZ ratio"));
  pt->AddText(Form("HighZ: [-100,-80] cm"));
  pt->AddText(Form("LowZ: [-20,0] cm"));
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.04);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->Draw();

  c1->Print("figure/pp_dedx_time_series_negative.pdf");
}

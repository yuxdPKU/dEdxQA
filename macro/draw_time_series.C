#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <sPhenixStyle.C>
#include <iostream>
#include <fstream>

void draw_time_series(std::string infile = "run_ratio.txt") {

  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1", "TGraph Example", 800, 600);
  
  std::ifstream inputFile(infile);

  if (!inputFile) {
      std::cerr << "Can not open file!" << std::endl;
      return;
  }

  std::vector<float> vec_runNo, vec_ratio;
  float runNo, ratio;

  while (inputFile >> runNo >> ratio) {
      if (std::isnan(runNo) || std::isnan(ratio)) {
        continue;
      }
      vec_runNo.push_back(runNo);
      vec_ratio.push_back(ratio);
  }

  std::cout<<"valid run number = "<<vec_runNo.size()<<std::endl;

  TGraph *graph = new TGraph(vec_runNo.size(), vec_runNo.data(), vec_ratio.data());
  
  graph->SetTitle("dE/dx time series in pp Run;Run Number;dE/dx ratio");
  graph->SetMarkerStyle(21);
  graph->SetMarkerColor(kBlue);
  graph->SetLineColor(kRed);
  graph->SetLineWidth(2);
  
  graph->Draw("ALP");

  //TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
  //legend->AddEntry(graph, "y = x^2", "lp");
  //legend->Draw();

  TLine *line = new TLine(51300, 1, 53200, 1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  line->Draw();

  TAxis *xAxis = graph->GetXaxis();
  xAxis->SetNdivisions(8);

  c1->Print("pp_dedx_time_series.pdf");
}

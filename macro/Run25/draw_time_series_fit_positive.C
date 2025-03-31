#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

void draw_time_series_fit_positive(std::string infile = "run_ratio_positive.txt") {

  gStyle->SetOptStat(0);

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

  std::vector<float> vec_runNo, vec_runtime, vec_ratio, vec_ratio_err, vec_ratio_high, vec_ratio_low;
  float runNo, runtime, ratio, ratio_err;
  int runmin=1000000, runmax=0;

  while (inputFile2 >> runNo >> ratio >> ratio_err)
  {
    if (std::isnan(runNo) || std::isnan(ratio))
    {
      continue;
    }
    if (ratio>2) continue;
//cout<<"runNo = "<<runNo<<endl;
    auto it = runtimestamp.find(runNo);
    if (it != runtimestamp.end())
    {
      if (runNo < runmin) {runmin = runNo;}
      if (runNo > runmax) {runmax = runNo;}
      vec_runtime.push_back(it->second);
      vec_runNo.push_back(runNo);
      vec_ratio.push_back(ratio);
      vec_ratio_err.push_back(ratio_err);
      vec_ratio_high.push_back(ratio+ratio_err);
      vec_ratio_low.push_back(ratio-ratio_err);
    }
    else
    {
      if (runNo < runmin) {runmin = runNo;}
      if (runNo > runmax) {runmax = runNo;}
      vec_runNo.push_back(runNo);
      vec_ratio.push_back(ratio);
      vec_ratio_err.push_back(ratio_err);
      vec_ratio_high.push_back(ratio+ratio_err);
      vec_ratio_low.push_back(ratio-ratio_err);
    }
  }

  int n = vec_runNo.size();
  std::cout<<"valid run number = "<<n<<std::endl;

  //TGraph *graph = new TGraph(vec_runNo.size(), vec_runtime.data(), vec_ratio.data());
  TGraph *graph = new TGraph(vec_runNo.size(), vec_runNo.data(), vec_ratio.data());
  TGraph* graph_high = new TGraph(vec_runNo.size(), vec_runNo.data(), vec_ratio_high.data());
  TGraph* graph_low= new TGraph(vec_runNo.size(), vec_runNo.data(), vec_ratio_low.data());

  //graph->SetTitle("dE/dx time series in pp Run;Date;dE/dx ratio");
  graph->SetTitle("dE/dx time series in pp Run;runNo;dE/dx ratio");
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.8);
  graph->SetMarkerColor(kBlue);
  graph->SetLineColor(kRed);
  graph->SetLineWidth(2);
  
  graph->Draw("ALP");

  TGraph* graph_band = new TGraph(2 * n);
  for (int i = 0; i < n; i++) {
    graph_band->SetPoint(i, vec_runNo[i], vec_ratio_high[i]);
    graph_band->SetPoint(n + i, vec_runNo[n - 1 - i], vec_ratio_low[n - 1 - i]);
  }
  graph_band->SetFillColorAlpha(kBlue, 0.3);
  graph_band->SetFillStyle(1001);
  graph_band->Draw("F");

  graph->Draw("LP SAME");

  //TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
  //legend->AddEntry(graph, "y = x^2", "lp");
  //legend->Draw();

  //auto xminIt = std::min_element(vec_runtime.begin(), vec_runtime.end());
  //auto xmaxIt = std::max_element(vec_runtime.begin(), vec_runtime.end());
  auto xminIt = std::min_element(vec_runNo.begin(), vec_runNo.end());
  auto xmaxIt = std::max_element(vec_runNo.begin(), vec_runNo.end());
  float xmin = *xminIt;
  float xmax = *xmaxIt;

  /*
  TF1 *expFitFunc = new TF1("expFitFunc", "[0]*exp([1]*(x-1723760957))", xmin, xmax);

  expFitFunc->SetParameters(1.5,-4.7000363e-05);
  expFitFunc->SetParLimits(0, -10, 10);
  expFitFunc->SetParLimits(1, -1e-02, 0);
  graph->Fit(expFitFunc);
  expFitFunc->Draw("same");

  float p0 = expFitFunc->GetParameter(0);
  float p1 = expFitFunc->GetParameter(1);
  */

  //TLine *line1 = new TLine(xmin-3000, 1, xmax+3000, 1);
  TLine *line1 = new TLine(xmin, 1, xmax, 1);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->SetLineColor(kRed);
  line1->Draw();

  TAxis *xAxis = graph->GetXaxis();
  xAxis->SetNdivisions(10);
  //xAxis->SetNdivisions(4);
  //xAxis->SetTimeDisplay(1);
  //xAxis->SetTimeFormat("%m/%d-%H");
  //xAxis->SetTimeOffset(0, "gmt");

  TPaveText *pt = new TPaveText(0.30, 0.65, 0.70, 0.9, "brNDC");
  pt->AddText("#it{#bf{sPHENIX}} Internal");
  pt->AddText("Cosmics");
  pt->AddText(Form("Run %d - %d",runmin,runmax));
  pt->AddText("Positive Z");
  //pt->AddText(Form("Exp Fit Function"));
  //pt->AddText(Form("%.2f#timesexp[-(x-1723760957)/%.0f]",p0,-1/p1));
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.04);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->Draw();

  c1->Print("figure/cosmics_dedx_time_series_zoomin_positive.pdf");
}

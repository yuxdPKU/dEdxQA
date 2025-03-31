#include <filesystem>
#include <sPhenixStyle.C>
#include <iostream>
#include <fstream>

std::vector<int> readNumberFromText(std::string infile);
void draw_dEdx_clusgz_1D_positive(std::string inputfilelist = "runlist")
{
  std::ofstream outputFile("run_ratio_positive.txt");

  if (!outputFile) {
      std::cerr << "Can not open file" << std::endl;
      return;
  }

  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  std::vector<int> runlist = readNumberFromText(inputfilelist);
  const int nruns = runlist.size();

  for (int irun = 0; irun < nruns; irun++)
  {
    int runnumber = runlist.at(irun);
    std::cout<<"Run "<<runnumber<<endl;

    TChain* chain = new TChain("tree");
    chain->Add(Form("root/clusters_seeds_%d-*_dedxqa.root",runnumber));

    int nevent  = chain->GetEntries();
    cout<<"total nevent = "<<nevent<<endl;

    if (nevent < 100)
    {
      cout<<"Statistics too low. Skip"<<endl;
      continue;
    }

    TH1F* h_track_dEdx_z9 = new TH1F("h_track_dEdx_z9","track dEdx",50,0,3000);
    TH1F* h_track_dEdx_z5 = new TH1F("h_track_dEdx_z5","track dEdx",50,0,3000);

    TCut cut = Form("charge*p<-0.2 && nintt>0");
    chain->Draw("ClusAdcPerLayerThickness_z9>>h_track_dEdx_z9",cut);
    chain->Draw("ClusAdcPerLayerThickness_z5>>h_track_dEdx_z5",cut);

    if (h_track_dEdx_z9->Integral() < 100 || h_track_dEdx_z5->Integral() < 100)
    {
      cout<<"h_track_dEdx_z9->Integral() = "<<h_track_dEdx_z9->Integral()<<" , h_track_dEdx_z5->Integral() = "<<h_track_dEdx_z5->Integral()<<endl;
      delete h_track_dEdx_z9;
      delete h_track_dEdx_z5;
      continue;
    }

    h_track_dEdx_z9->SetMinimum();
    h_track_dEdx_z5->SetMinimum();

    h_track_dEdx_z9->Scale(h_track_dEdx_z5->Integral() / h_track_dEdx_z9->Integral());

    float ymax = h_track_dEdx_z9->GetMaximum() > h_track_dEdx_z5->GetMaximum() ? h_track_dEdx_z9->GetMaximum() : h_track_dEdx_z5->GetMaximum();

    TCanvas* can = new TCanvas("can","",800,600);
    can->cd();

    h_track_dEdx_z9->SetMaximum(1.1*ymax);
    h_track_dEdx_z9->SetLineColor(kRed);
    h_track_dEdx_z9->Draw("hist");
    h_track_dEdx_z9->GetXaxis()->SetTitle("Mean cluster Adc corrected by path length");
    h_track_dEdx_z9->GetYaxis()->SetTitle("Events");
    h_track_dEdx_z5->SetLineColor(kBlue);
    h_track_dEdx_z5->Draw("hist,same");

    TF1 *landauFit_z9 = new TF1("landauFit_z9", "[2] * TMath::Landau(x,[0],[1])", 0, 3000);
    TF1 *landauFit_z5 = new TF1("landauFit_z5", "[2] * TMath::Landau(x,[0],[1])", 0, 3000);
    landauFit_z9->SetParameters(0.5*h_track_dEdx_z9->GetMean(), 0.5*h_track_dEdx_z9->GetRMS(), h_track_dEdx_z9->Integral());
    landauFit_z5->SetParameters(0.5*h_track_dEdx_z5->GetMean(), 0.5*h_track_dEdx_z5->GetRMS(), h_track_dEdx_z5->Integral());
    landauFit_z9->SetParLimits(0, 0, 3000);
    landauFit_z5->SetParLimits(0, 0, 3000);
    landauFit_z9->SetParLimits(1, 0.1*h_track_dEdx_z9->GetRMS(), 5*h_track_dEdx_z9->GetRMS());
    landauFit_z5->SetParLimits(1, 0.1*h_track_dEdx_z5->GetRMS(), 5*h_track_dEdx_z5->GetRMS());
    landauFit_z9->SetParLimits(2, 0, 2*h_track_dEdx_z9->Integral());
    landauFit_z5->SetParLimits(2, 0, 2*h_track_dEdx_z5->Integral());
    h_track_dEdx_z9->Fit("landauFit_z9");
    h_track_dEdx_z5->Fit("landauFit_z5");

    landauFit_z9->SetLineColor(kRed+1);
    landauFit_z9->SetLineWidth(2);
    landauFit_z5->SetLineColor(kBlue+1);
    landauFit_z5->SetLineWidth(2);

    landauFit_z9->Draw("same");
    landauFit_z5->Draw("same");

    //float xpos_z9 = h_track_dEdx_z9->GetBinCenter( h_track_dEdx_z9->GetMaximumBin() );
    //float xpos_z5 = h_track_dEdx_z5->GetBinCenter( h_track_dEdx_z5->GetMaximumBin() );
    float xpos_z9 = landauFit_z9->GetParameter(0);
    float xpos_z5 = landauFit_z5->GetParameter(0);
    float z9_to_z5_ratio = xpos_z9 / xpos_z5;

    if (isnan(z9_to_z5_ratio))
    {
      continue;
    }
    outputFile << runnumber << " " << z9_to_z5_ratio << std::endl;

    TPaveText *pt = new TPaveText(0.48, 0.45, 0.9, 0.9, "brNDC");
    pt->AddText("#it{#bf{sPHENIX}} Internal");
    pt->AddText("p+p #sqrt{s}=200 GeV");
    pt->AddText(Form("Run %d",runnumber));
    pt->AddText("Negative tracks p > 0.2 GeV");
    pt->AddText("Silicon matching (nintt > 0)");
    pt->AddText(Form("Laudau Fit"));
    pt->AddText(Form("Mean #color[2]{HighZ} / #color[4]{LowZ} = %.2f",z9_to_z5_ratio));
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.05);
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->Draw();

    TLegend *legend = new TLegend(0.48, 0.30, 0.9, 0.45);
    legend->AddEntry(h_track_dEdx_z9, "Cluster Z in [80,100] cm", "l");
    legend->AddEntry(h_track_dEdx_z5, "Cluster Z in [0.20] cm", "l");
    legend->SetTextSize(0.04);
    legend->Draw();

    TLine *line_z9 = new TLine(xpos_z9, 0, xpos_z9, 1.1*ymax);
    line_z9->SetLineStyle(2);
    line_z9->SetLineWidth(2);
    line_z9->SetLineColor(kRed);
    line_z9->Draw();

    TLine *line_z5 = new TLine(xpos_z5, 0, xpos_z5, 1.1*ymax);
    line_z5->SetLineStyle(2);
    line_z5->SetLineWidth(2);
    line_z5->SetLineColor(kBlue);
    line_z5->Draw();

    can->Update();
    can->Print(Form("figure/track_dEdx_run%d_positive.pdf",runnumber));
  
    delete can;
    delete h_track_dEdx_z9;
    delete h_track_dEdx_z5;
  }

  outputFile.close();

}

std::vector<int> readNumberFromText(std::string infile) {

    std::vector<int> vector;
    vector.clear();

    std::ifstream inputFile(infile);
    if (!inputFile.is_open()) {
        std::cerr << "Can not open file!" << std::endl;
        return vector;
    }

    int number;
    while (inputFile >> number) {
        vector.push_back(number);
    }

    inputFile.close();
    return vector;
}

#include <filesystem>
#include <sPhenixStyle.C>
#include <iostream>
#include <fstream>

std::vector<int> readNumberFromText(std::string infile);
void draw_dEdx_clusgz_1D_negative(std::string inputfilelist = "runlist")
{
  std::ofstream outputFile("run_ratio_negative.txt");

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
if (runnumber!=51245 && runnumber!=51249 && runnumber!=51255) continue;

    TChain* chain = new TChain("tree");
    chain->Add(Form("root/clusters_seeds_%d-*_dedxqa.root",runnumber));

    int nevent  = chain->GetEntries();
    cout<<"total nevent = "<<nevent<<endl;

    if (nevent < 100)
    {
      cout<<"Statistics too low. Skip"<<endl;
      continue;
    }

    TH1F* h_track_dEdx_z0 = new TH1F("h_track_dEdx_z0","track dEdx",50,0,3000);
    TH1F* h_track_dEdx_z4 = new TH1F("h_track_dEdx_z4","track dEdx",50,0,3000);

    TCut cut = Form("charge*p<-0.2 && nintt>0");
    chain->Draw("ClusAdcPerLayerThickness_z0>>h_track_dEdx_z0",cut);
    chain->Draw("ClusAdcPerLayerThickness_z4>>h_track_dEdx_z4",cut);

    if (h_track_dEdx_z0->Integral() < 100 || h_track_dEdx_z4->Integral() < 100)
    {
      cout<<"h_track_dEdx_z0->Integral() = "<<h_track_dEdx_z0->Integral()<<" , h_track_dEdx_z4->Integral() = "<<h_track_dEdx_z4->Integral()<<endl;
      delete h_track_dEdx_z0;
      delete h_track_dEdx_z4;
      continue;
    }

    h_track_dEdx_z0->SetMinimum();
    h_track_dEdx_z4->SetMinimum();

    h_track_dEdx_z0->Scale(h_track_dEdx_z4->Integral() / h_track_dEdx_z0->Integral());

    float ymax = h_track_dEdx_z0->GetMaximum() > h_track_dEdx_z4->GetMaximum() ? h_track_dEdx_z0->GetMaximum() : h_track_dEdx_z4->GetMaximum();

    TCanvas* can = new TCanvas("can","",800,600);
    can->cd();

    h_track_dEdx_z0->SetMaximum(1.1*ymax);
    h_track_dEdx_z0->SetLineColor(kRed);
    h_track_dEdx_z0->Draw("hist");
    h_track_dEdx_z0->GetXaxis()->SetTitle("Mean cluster Adc corrected by path length");
    h_track_dEdx_z0->GetYaxis()->SetTitle("Events");
    h_track_dEdx_z4->SetLineColor(kBlue);
    h_track_dEdx_z4->Draw("hist,same");

    TF1 *landauFit_z0 = new TF1("landauFit_z0", "[2] * TMath::Landau(x,[0],[1])", 0, 3000);
    TF1 *landauFit_z4 = new TF1("landauFit_z4", "[2] * TMath::Landau(x,[0],[1])", 0, 3000);
    landauFit_z0->SetParameters(0.8*h_track_dEdx_z0->GetMean(), 0.5*h_track_dEdx_z0->GetRMS(), h_track_dEdx_z0->Integral());
    landauFit_z4->SetParameters(0.8*h_track_dEdx_z4->GetMean(), 0.5*h_track_dEdx_z4->GetRMS(), h_track_dEdx_z4->Integral());
    landauFit_z0->SetParLimits(0, 0.05*h_track_dEdx_z0->GetMean(), 2*h_track_dEdx_z0->GetMean());
    landauFit_z4->SetParLimits(0, 0.05*h_track_dEdx_z4->GetMean(), 2*h_track_dEdx_z4->GetMean());
    landauFit_z0->SetParLimits(1, 0.1*h_track_dEdx_z0->GetRMS(), 5*h_track_dEdx_z0->GetRMS());
    landauFit_z4->SetParLimits(1, 0.1*h_track_dEdx_z4->GetRMS(), 5*h_track_dEdx_z4->GetRMS());
    landauFit_z0->SetParLimits(2, 0, 2*h_track_dEdx_z0->Integral());
    landauFit_z4->SetParLimits(2, 0, 2*h_track_dEdx_z4->Integral());
    h_track_dEdx_z0->Fit("landauFit_z0");
    h_track_dEdx_z4->Fit("landauFit_z4");

    landauFit_z0->SetLineColor(kRed+1);
    landauFit_z0->SetLineWidth(2);
    landauFit_z4->SetLineColor(kBlue+1);
    landauFit_z4->SetLineWidth(2);

    landauFit_z0->Draw("same");
    landauFit_z4->Draw("same");

    //float xpos_z0 = h_track_dEdx_z0->GetBinCenter( h_track_dEdx_z0->GetMaximumBin() );
    //float xpos_z4 = h_track_dEdx_z4->GetBinCenter( h_track_dEdx_z4->GetMaximumBin() );
    float xpos_z0 = landauFit_z0->GetParameter(0);
    float xpos_z4 = landauFit_z4->GetParameter(0);
    float z0_to_z4_ratio = xpos_z0 / xpos_z4;

    if (isnan(z0_to_z4_ratio))
    {
      continue;
    }
    outputFile << runnumber << " " << z0_to_z4_ratio << std::endl;

    TPaveText *pt = new TPaveText(0.48, 0.45, 0.9, 0.9, "brNDC");
    pt->AddText("#it{#bf{sPHENIX}} Internal");
    pt->AddText("p+p #sqrt{s}=200 GeV");
    pt->AddText(Form("Run %d",runnumber));
    pt->AddText("Negative tracks p > 0.2 GeV");
    pt->AddText("Silicon matching (nintt > 0)");
    pt->AddText(Form("Laudau Fit"));
    pt->AddText(Form("Mean #color[2]{HighZ} / #color[4]{LowZ} = %.2f",z0_to_z4_ratio));
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.05);
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->Draw();

    TLegend *legend = new TLegend(0.48, 0.30, 0.9, 0.45);
    legend->AddEntry(h_track_dEdx_z0, "Cluster Z in [-100,-80] cm", "l");
    legend->AddEntry(h_track_dEdx_z4, "Cluster Z in [-20,0] cm", "l");
    legend->SetTextSize(0.04);
    legend->Draw();

    TLine *line_z0 = new TLine(xpos_z0, 0, xpos_z0, 1.1*ymax);
    line_z0->SetLineStyle(2);
    line_z0->SetLineWidth(2);
    line_z0->SetLineColor(kRed);
    line_z0->Draw();

    TLine *line_z4 = new TLine(xpos_z4, 0, xpos_z4, 1.1*ymax);
    line_z4->SetLineStyle(2);
    line_z4->SetLineWidth(2);
    line_z4->SetLineColor(kBlue);
    line_z4->Draw();

    can->Update();
    can->Print(Form("figure/track_dEdx_run%d_negative.pdf",runnumber));
  
    delete can;
    delete h_track_dEdx_z0;
    delete h_track_dEdx_z4;
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

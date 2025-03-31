#include <filesystem>
#include <iostream>
#include <fstream>

std::vector<int> readNumberFromText(std::string infile);
void draw_dEdx_clusgz_1D(std::string inputfilelist = "runlist")
{
  gStyle->SetOptStat(0);

  std::ofstream outputFile_neg("run_ratio_negative.txt");
  std::ofstream outputFile_pos("run_ratio_positive.txt");

  if (!outputFile_neg || !outputFile_pos) {
      std::cerr << "Can not open file" << std::endl;
      return;
  }

  std::vector<int> runlist = readNumberFromText(inputfilelist);
  const int nruns = runlist.size();

  for (int irun = 0; irun < nruns; irun++)
  {
    int runnumber = runlist.at(irun);
    std::cout<<"Run "<<runnumber<<endl;

    std::stringstream nice_runnumber;
    nice_runnumber << std::setw(8) << std::setfill('0') << to_string(runnumber);

    int rounded_up = 100*(std::ceil((float) runnumber/100));
    std::stringstream nice_rounded_up;
    nice_rounded_up << std::setw(8) << std::setfill('0') << to_string(rounded_up);

    int rounded_down = 100*(std::floor((float) runnumber/100));
    std::stringstream nice_rounded_down;
    nice_rounded_down << std::setw(8) << std::setfill('0') << to_string(rounded_down);

    std::string filename = "/sphenix/data/data02/sphnxpro/QAhtml/aggregated/run3auau/cosmics/new_2024p017_v000/DST_TRKR_SEED/run_" + nice_rounded_down.str() + "_" + nice_rounded_up.str() + "/HIST_DST_TRKR_SEED_run3auau_new_2024p017_v000-" + nice_runnumber.str() + "-9999.root";
    TFile* file = new TFile(filename.c_str(),"");
    if (!file)
    {
      std::cout << "No root file for " << runnumber << std::endl;
      continue;
    }

    TH2 *h_dedx_pq[10];
    for (int i = 0; i < 10; i++)
    {
      h_dedx_pq[i] = (TH2*) file->Get(Form("h_TpcSeedsQA_dedx_pq_%d",i));
    }

    if (h_dedx_pq[0] && h_dedx_pq[4])
    {
      int xBinMin = h_dedx_pq[0]->GetXaxis()->FindBin(-3.);
      int xBinMax = h_dedx_pq[0]->GetXaxis()->FindBin(0.2);

      TH1D *h_dedx_z0 = h_dedx_pq[0]->ProjectionY("h_dedx_z0", xBinMin, xBinMax);
      TH1D *h_dedx_z4 = h_dedx_pq[4]->ProjectionY("h_dedx_z4", xBinMin, xBinMax);

      h_dedx_z0->GetXaxis()->SetRangeUser(0, 3000);
      h_dedx_z4->GetXaxis()->SetRangeUser(0, 3000);

      h_dedx_z0->Scale(h_dedx_z4->Integral() / h_dedx_z0->Integral());

      float ymax = h_dedx_z0->GetMaximum() > h_dedx_z4->GetMaximum() ? h_dedx_z0->GetMaximum() : h_dedx_z4->GetMaximum();

      TCanvas* can = new TCanvas("can","",800,600);
      can->cd();
      gPad->SetRightMargin(0.05);

      h_dedx_z0->SetMaximum(1.1*ymax);
      h_dedx_z0->SetLineColor(kRed);
      h_dedx_z0->Draw("hist");
      h_dedx_z0->GetXaxis()->SetTitle("Mean cluster Adc corrected by path length");
      h_dedx_z0->GetYaxis()->SetTitle("Events");
      h_dedx_z4->SetLineColor(kBlue);
      h_dedx_z4->Draw("hist,same");

      TF1 *landauFit_z0 = new TF1("landauFit_z0", "[2] * TMath::Landau(x,[0],[1])", 0, 3000);
      TF1 *landauFit_z4 = new TF1("landauFit_z4", "[2] * TMath::Landau(x,[0],[1])", 0, 3000);
      landauFit_z0->SetParameters(0.7*h_dedx_z0->GetMean(), 0.5*h_dedx_z0->GetRMS(), h_dedx_z0->Integral());
      landauFit_z4->SetParameters(0.7*h_dedx_z4->GetMean(), 0.5*h_dedx_z4->GetRMS(), h_dedx_z4->Integral());
      landauFit_z0->SetParLimits(0, 0.05*h_dedx_z0->GetMean(), 2*h_dedx_z0->GetMean());
      landauFit_z4->SetParLimits(0, 0.05*h_dedx_z4->GetMean(), 2*h_dedx_z4->GetMean());
      landauFit_z0->SetParLimits(1, 0.1*h_dedx_z0->GetRMS(), 5*h_dedx_z0->GetRMS());
      landauFit_z4->SetParLimits(1, 0.1*h_dedx_z4->GetRMS(), 5*h_dedx_z4->GetRMS());
      landauFit_z0->SetParLimits(2, 0, 2*h_dedx_z0->Integral());
      landauFit_z4->SetParLimits(2, 0, 2*h_dedx_z4->Integral());
      h_dedx_z0->Fit("landauFit_z0","Q");
      h_dedx_z4->Fit("landauFit_z4","Q");

      landauFit_z0->SetLineColor(kRed+1);
      landauFit_z0->SetLineWidth(2);
      landauFit_z4->SetLineColor(kBlue+1);
      landauFit_z4->SetLineWidth(2);

      landauFit_z0->Draw("same");
      landauFit_z4->Draw("same");

      float xpos_z0 = landauFit_z0->GetParameter(0);
      float xpos_z4 = landauFit_z4->GetParameter(0);
      float xpos_z0_err = landauFit_z0->GetParError(0);
      float xpos_z4_err = landauFit_z4->GetParError(0);
      float z0_to_z4_ratio = xpos_z0 / xpos_z4;
      float z0_to_z4_ratio_err = z0_to_z4_ratio * sqrt(pow(xpos_z0_err/xpos_z0,2) + pow(xpos_z4_err/xpos_z4,2));
      if (std::isnan(z0_to_z4_ratio))
      {
	continue;
      }
      outputFile_neg << runnumber << " " << z0_to_z4_ratio << z0_to_z4_ratio_err << std::endl;

      TPaveText *pt = new TPaveText(0.48, 0.6, 0.9, 0.9, "brNDC");
      pt->AddText("#it{#bf{sPHENIX}} Internal");
      pt->AddText("Cosmics");
      pt->AddText(Form("Run %d",runnumber));
      pt->AddText("Negative tracks p > 0.2 GeV");
      //pt->AddText(Form("Laudau Fit"));
      pt->AddText(Form("Mean #color[2]{HighZ} / #color[4]{LowZ} = %.2f",z0_to_z4_ratio));
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->SetTextSize(0.05);
      pt->SetFillStyle(0);
      pt->SetBorderSize(0);
      pt->Draw();

      TLegend *legend = new TLegend(0.48, 0.45, 0.9, 0.6);
      legend->AddEntry(h_dedx_z0, "Cluster Z in [-100,-80] cm", "l");
      legend->AddEntry(h_dedx_z4, "Cluster Z in [-20,0] cm", "l");
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
      delete h_dedx_z0;
      delete h_dedx_z4;
    }
    else
    {
        // histogram is missing
	std::cout << "No histogram" << std::endl;
        continue;
    }

    if (h_dedx_pq[9] && h_dedx_pq[5])
    {
      int xBinMin = h_dedx_pq[9]->GetXaxis()->FindBin(-3.);
      int xBinMax = h_dedx_pq[9]->GetXaxis()->FindBin(0.2);

      TH1D *h_dedx_z9 = h_dedx_pq[9]->ProjectionY("h_dedx_z9", xBinMin, xBinMax);
      TH1D *h_dedx_z5 = h_dedx_pq[5]->ProjectionY("h_dedx_z5", xBinMin, xBinMax);

      h_dedx_z9->GetXaxis()->SetRangeUser(0, 3000);
      h_dedx_z5->GetXaxis()->SetRangeUser(0, 3000);

      h_dedx_z9->Scale(h_dedx_z5->Integral() / h_dedx_z9->Integral());

      float ymax = h_dedx_z9->GetMaximum() > h_dedx_z5->GetMaximum() ? h_dedx_z9->GetMaximum() : h_dedx_z5->GetMaximum();

      TCanvas* can = new TCanvas("can","",800,600);
      can->cd();
      gPad->SetRightMargin(0.05);

      h_dedx_z9->SetMaximum(1.1*ymax);
      h_dedx_z9->SetLineColor(kRed);
      h_dedx_z9->Draw("hist");
      h_dedx_z9->GetXaxis()->SetTitle("Mean cluster Adc corrected by path length");
      h_dedx_z9->GetYaxis()->SetTitle("Events");
      h_dedx_z5->SetLineColor(kBlue);
      h_dedx_z5->Draw("hist,same");

      TF1 *landauFit_z9 = new TF1("landauFit_z9", "[2] * TMath::Landau(x,[0],[1])", 0, 3000);
      TF1 *landauFit_z5 = new TF1("landauFit_z5", "[2] * TMath::Landau(x,[0],[1])", 0, 3000);
      landauFit_z9->SetParameters(0.8*h_dedx_z9->GetMean(), 0.5*h_dedx_z9->GetRMS(), h_dedx_z9->Integral());
      landauFit_z5->SetParameters(0.8*h_dedx_z5->GetMean(), 0.5*h_dedx_z5->GetRMS(), h_dedx_z5->Integral());
      landauFit_z9->SetParLimits(0, 0.05*h_dedx_z9->GetMean(), 2*h_dedx_z9->GetMean());
      landauFit_z5->SetParLimits(0, 0.05*h_dedx_z5->GetMean(), 2*h_dedx_z5->GetMean());
      landauFit_z9->SetParLimits(1, 0.1*h_dedx_z9->GetRMS(), 5*h_dedx_z9->GetRMS());
      landauFit_z5->SetParLimits(1, 0.1*h_dedx_z5->GetRMS(), 5*h_dedx_z5->GetRMS());
      landauFit_z9->SetParLimits(2, 0, 2*h_dedx_z9->Integral());
      landauFit_z5->SetParLimits(2, 0, 2*h_dedx_z5->Integral());
      h_dedx_z9->Fit("landauFit_z9","Q");
      h_dedx_z5->Fit("landauFit_z5","Q");

      landauFit_z9->SetLineColor(kRed+1);
      landauFit_z9->SetLineWidth(2);
      landauFit_z5->SetLineColor(kBlue+1);
      landauFit_z5->SetLineWidth(2);

      landauFit_z9->Draw("same");
      landauFit_z5->Draw("same");

      float xpos_z9 = landauFit_z9->GetParameter(0);
      float xpos_z5 = landauFit_z5->GetParameter(0);
      float xpos_z9_err = landauFit_z9->GetParError(0);
      float xpos_z5_err = landauFit_z5->GetParError(0);
      float z9_to_z5_ratio = xpos_z9 / xpos_z5;
      float z9_to_z5_ratio_err = z9_to_z5_ratio * sqrt(pow(xpos_z9_err/xpos_z9,2) + pow(xpos_z5_err/xpos_z5,2));
      if (std::isnan(z9_to_z5_ratio))
      {
	continue;
      }
      outputFile_pos << runnumber << " " << z9_to_z5_ratio << z9_to_z5_ratio_err << std::endl;

      TPaveText *pt = new TPaveText(0.48, 0.6, 0.9, 0.9, "brNDC");
      pt->AddText("#it{#bf{sPHENIX}} Internal");
      pt->AddText("Cosmics");
      pt->AddText(Form("Run %d",runnumber));
      pt->AddText("Negative tracks p > 0.2 GeV");
      //pt->AddText(Form("Laudau Fit"));
      pt->AddText(Form("Mean #color[2]{HighZ} / #color[4]{LowZ} = %.2f",z9_to_z5_ratio));
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->SetTextSize(0.05);
      pt->SetFillStyle(0);
      pt->SetBorderSize(0);
      pt->Draw();

      TLegend *legend = new TLegend(0.48, 0.45, 0.9, 0.6);
      legend->AddEntry(h_dedx_z9, "Cluster Z in [80,100] cm", "l");
      legend->AddEntry(h_dedx_z5, "Cluster Z in [0,20] cm", "l");
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
      delete h_dedx_z9;
      delete h_dedx_z5;
    }
    else
    {
        // histogram is missing
	std::cout << "No histogram" << std::endl;
	continue;
    }
  }
  outputFile_neg.close();
  outputFile_pos.close();
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

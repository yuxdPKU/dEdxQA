#include <filesystem>
#include <sPhenixStyle.C>

std::vector<int> readNumberFromText(std::string infile);
void draw_dEdx_clusgz_1D(std::string inputfilelist = "runList/runlist.partial")
{
  int verbosity = 0;

  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  std::vector<int> runlist = readNumberFromText(inputfilelist);
  const int nruns = runlist.size();

  for (int irun = 0; irun < nruns; irun++)
  {
    int runnumber = runlist.at(irun);

    TChain* chain = new TChain("tree");
    chain->Add(Form("root/clusters_seeds_%d-*_dedxqa.root",runnumber));

    int nevent  = chain->GetEntries();
    cout<<"total nevent = "<<nevent<<endl;

    TH1F* h_track_dEdx_z0 = new TH1F("h_track_dEdx_z0","track dEdx",50,0,3000);
    TH1F* h_track_dEdx_z4 = new TH1F("h_track_dEdx_z4","track dEdx",50,0,3000);

    TCut cut = Form("charge*p<-0.2 && nintt>0");
    chain->Draw("ClusAdcPerLayerThickness_z0>>h_track_dEdx_z0",cut);
    chain->Draw("ClusAdcPerLayerThickness_z4>>h_track_dEdx_z4",cut);

    h_track_dEdx_z0->SetMinimum();
    h_track_dEdx_z4->SetMinimum();

    float xpos_z0 = h_track_dEdx_z0->GetBinCenter( h_track_dEdx_z0->GetMaximumBin() );
    float xpos_z4 = h_track_dEdx_z4->GetBinCenter( h_track_dEdx_z4->GetMaximumBin() );
    float z0_to_z4_ratio = xpos_z0 / xpos_z4;

    h_track_dEdx_z0->Scale(h_track_dEdx_z4->Integral() / h_track_dEdx_z0->Integral());

    float ymax = h_track_dEdx_z0->GetMaximum() > h_track_dEdx_z4->GetMaximum() ? h_track_dEdx_z0->GetMaximum() : h_track_dEdx_z4->GetMaximum();

    TCanvas* can = new TCanvas("can","",800,600);
    can->cd();
  
    h_track_dEdx_z0->SetMaximum(1.1*ymax);
    h_track_dEdx_z0->SetLineColor(kRed);
    h_track_dEdx_z0->Draw("hist");
    h_track_dEdx_z0->GetXaxis()->SetTitle("Events");
    h_track_dEdx_z0->GetYaxis()->SetTitle("Mean cluster Adc corrected by path length");
    h_track_dEdx_z4->SetLineColor(kBlue);
    h_track_dEdx_z4->Draw("hist,same");
  
    TPaveText *pt = new TPaveText(0.48, 0.5, 0.9, 0.9, "brNDC");
    pt->AddText("#it{#bf{sPHENIX}} Internal");
    pt->AddText("p+p #sqrt{s}=200 GeV");
    pt->AddText(Form("Run %d",runnumber));
    pt->AddText("Negative tracks p > 0.2 GeV");
    pt->AddText("Silicon matching (nintt > 0)");
    pt->AddText(Form("#color[2]{HighZ} / #color[4]{LowZ} = %.2f",z0_to_z4_ratio));
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.05);
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->Draw();

    TLegend *legend = new TLegend(0.48, 0.35, 0.9, 0.5);
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
    can->Print(Form("figure/track_dEdx_run%d.pdf",runnumber));
  
    delete can;
    delete h_track_dEdx_z0;
    delete h_track_dEdx_z4;
  }

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

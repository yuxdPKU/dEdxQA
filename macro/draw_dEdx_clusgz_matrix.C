#include <filesystem>
#include <sPhenixStyle.C>

float Cal_Asy(float NF, float NB);

void draw_dEdx_clusgz_matrix(int runnumber=51245)
{
  int verbosity = 0;

  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  TChain* chain = new TChain("tree");
  chain->Add(Form("root/clusters_seeds_%d-*_dedxqa.root",runnumber));

  int nevent  = chain->GetEntries();
  cout<<"total nevent = "<<nevent<<endl;

  for (int i=0; i<10; i++)
  {
    TH2F* h_track_dEdx_pq = new TH2F("h_track_dEdx_pq","track dEdx vs p*q",100,-3,3,100,0,10000);
    float gzmin = -100+i*20;
    float gzmax = -100+(i+1)*20;
    TCut cut = Form("nintt>0");
    chain->Draw(Form("ClusAdcPerLayerThickness_z%d:p*charge>>h_track_dEdx_pq",i),cut);

    TCanvas* can = new TCanvas("can","",800,600);
    can->cd();

    h_track_dEdx_pq->Draw("colz");
    h_track_dEdx_pq->GetXaxis()->SetTitle("p#timesq");
    h_track_dEdx_pq->GetYaxis()->SetTitle("Corrected Mean Adc / cluster");

    TPaveText *pt = new TPaveText(0.2, 0.7, 0.45, 0.9, "brNDC");
    pt->AddText("#it{#bf{sPHENIX}} Internal");
    pt->AddText("p+p #sqrt{s}=200 GeV");
    pt->AddText(Form("Run %d",runnumber));
    pt->AddText(Form("Cluster gz in [%d,%d] cm",(int)gzmin,(int)gzmax));
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.05);
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->Draw();

    can->Update();
    can->Print(Form("figure/track_dEdx_pq_run%d_%d.pdf",runnumber,i));

    delete can;
    delete h_track_dEdx_pq;
  }

}

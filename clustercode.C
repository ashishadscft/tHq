void totalcountcluster1(){
  
  TFile *f = new TFile("./TDRplotscluster.root","RECREATE");
  
  f->cd();
  gDirectory->pwd();
  f->ls();
  
  TString inpath = "/home/ashish/TEPX_rootfiles/samples_17Feb2020/";
  gROOT->ProcessLine(".x /home/ashish/rootlogon.C");
  
  //string containing names of input sample files
  std::vector<std::string> pulist = {//"2023D42noPU" , 
    "2023D42PU0p5", "2023D42PU1", "2023D42PU1p5", "2023D42PU2", "2023D42PU10", "2023D42PU30", "2023D42PU50", "2023D42PU100", "2023D42PU140", "2023D42PU200" };
  std::map<std::string, float> pumap = {//{ "2023D42noPU", 0 },
    { "2023D42PU0p5", 0.5 }, { "2023D42PU1", 1 }, { "2023D42PU1p5", 1.5 }, { "2023D42PU2", 2 }, { "2023D42PU10", 10 }, { "2023D42PU30", 30 }, { "2023D42PU50", 50 }, { "2023D42PU100", 100 }, { "2023D42PU140", 140 }, { "2023D42PU200", 200 } };
  
  std::vector<std::string> disklist = {"-4", "-3", "-2", "-1", "1", "2", "3", "4"};
  
  //create the profiles to be filled below.
  TGraphErrors* TEPXClustersPerEvent; //number of clusters vs pu
  TGraphErrors* TEPXClustersPerEvent1[8][5];
  TGraphErrors* TEPXClustersPerEvent_disk[8];
  TGraphErrors* TEPXClustersPerEvent_diskcombined[4];
  TGraphErrors* TEPXClustersPerEvent_diskcombinedring[8][5];
  TGraphErrors* TEPXClustersPerEvent_D4R1;
  
  TH2F* Histogram2D[8][5];
  
  //Non-linearity graphs
  TGraphErrors* NonLinearity_TEPXClustersPerEvent; //number of clusters vs pu 
  NonLinearity_TEPXClustersPerEvent = new TGraphErrors();
  
  TEPXClustersPerEvent = new TGraphErrors();
  TEPXClustersPerEvent_D4R1 = new TGraphErrors();
  
  for (int d = 0; d < disklist.size(); d++){
    
    TEPXClustersPerEvent_disk[d]= new TGraphErrors();
    
    if(disklist[d] > 0) {
      TEPXClustersPerEvent_diskcombined[d]= new TGraphErrors();
    }
    
    for (int r = 0; r < 5; r++){
      
      if(disklist[d] > 0) {
	TEPXClustersPerEvent_diskcombinedring[d][r]= new TGraphErrors();	
      }
      TEPXClustersPerEvent1[d][r] = new TGraphErrors();
      
    }  
  }
  
  
  //read the histograms
  TProfile* Prof_TEPXClustersPerEvent[20][20]; //array pu,disk
  for (int pu = 0; pu < pulist.size(); pu++) {
    
    TFile F(inpath + pulist[pu].c_str() + ".root", "read");
    gROOT->cd();
    
    float totalcount = 0;
    float totalcounterror = 0;
    
    float totalcount_D4R1 = 0;
    float totalcounterror_D4R1 = 0;
    
    float totalcount_disk[8] = {0,0,0,0,0,0,0,0};
    float totalcounterror_disk[8] = {0,0,0,0,0,0,0,0};
    
    float totalcount_diskcombined[4] = {0,0,0,0};
    float totalcounterror_diskcombined[4] = {0,0,0,0};
    
    float totalcount_diskcombinedring[4][5] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
    float totalcounterror_diskcombinedring[4][5] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
    
    for (int d = 0; d < disklist.size(); d++) {
      
      TString histoname = "BRIL_IT_Analysis/TEPX/Clusters/Number of clusters for Disk ";
      
      TH2F* H = (TH2F*)F.Get(histoname + disklist[d].c_str());
      Histogram2D[pu][d] = (TH2F*)H->Clone(TString(H->GetName()));
      Prof_TEPXClustersPerEvent[pu][d] = (TProfile*)H->ProfileX()->Clone(TString(H->GetName()) + "Profile");
      
      for (int r = 0; r < 5; r++) {
	
	float count = Prof_TEPXClustersPerEvent[pu][d]->GetBinContent(r + 1);
	float counterror = Prof_TEPXClustersPerEvent[pu][d]->GetBinError(r + 1);
	
	TEPXClustersPerEvent1[d][r]->SetPoint(pu, pumap[pulist[pu]], count);
	TEPXClustersPerEvent1[d][r]->SetPointError(pu, 0, counterror);
	
      }
      
      totalcount_disk[d]+= count;     //not correct. Think how to find/sum counts per disk ??
      totalcounterror_disk[d]+=counterror;   //not correct. does not make sense an array equated with a number.
      
    }
    
    
    totalcount+=count;
    totalcounterror+= counterror * counterror;
    
    
    for (int d = 0; d < disklist.size(); d++) {
      for (int r = 0; r < 5; r++) {
	
	TEPXClustersPerEvent_disk[d]->SetPoint(pu, pumap[pulist[pu]], totalcount_disk[d]);
	TEPXClustersPerEvent_disk[d]->SetPointError(pu, 0, sqrt(totalcounterror_disk[d]));
	
	if(disklist[d] > 0) {
	  
	  TEPXClustersPerEvent_diskcombined[d]->SetPoint(pu, pumap[pulist[pu]], totalcount_diskcombined[d]);
	  TEPXClustersPerEvent_diskcombined[d]->SetPointError(pu, 0, sqrt(totalcounterror_diskcombined[d]));
	  
	  TEPXClustersPerEvent_diskcombinedring[d][r]->SetPoint(pu, pumap[pulist[pu]],totalcount_diskcombinedring[d][r]);
	  TEPXClustersPerEvent_diskcombinedring[d][r]->SetPointError(pu, 0, sqrt(totalcounterror_diskcombinedring[d][r]));
	  
	}
	
	TEPXClustersPerEvent->SetPoint(pu, pumap[pulist[pu]], totalcount);
	TEPXClustersPerEvent->SetPointError(pu, 0, sqrt(totalcounterror));
	
	TEPXClustersPerEvent_D4R1->SetPoint(pu, pumap[pulist[pu]], totalcount_D4R1);
	TEPXClustersPerEvent_D4R1->SetPointError(pu, 0, sqrt(totalcounterror_D4R1));
	
      }
    }
  }  
  
  TEPXClustersPerEvent->SetMarkerSize(2);
  TEPXClustersPerEvent->SetMarkerColor(2);
  TEPXClustersPerEvent->SetName(TString("clusters_TEPX"));  
  f->WriteTObject(TEPXClustersPerEvent);
  
  TEPXClustersPerEvent_D4R1->SetMarkerSize(2);
  TEPXClustersPerEvent_D4R1->SetMarkerColor(2);
  TEPXClustersPerEvent_D4R1->SetName(TString("clusters_TEPXD4R1"));  
  f->WriteTObject(TEPXClustersPerEvent_D4R1);
  
  for (int d = 0; d < disklist.size(); d++) {
    
    TEPXClustersPerEvent_disk[d]->SetMarkerSize(2);
    TEPXClustersPerEvent_disk[d]->SetMarkerColor(2);
    TEPXClustersPerEvent_disk[d]->SetName(TString("clusters_D") + disklist[d].c_str());
    f->WriteTObject(TEPXClustersPerEvent_disk[d]);
    
    if(disklist[d] > 0) {
      
      TEPXClustersPerEvent_diskcombined[d]-> SetMarkerSize(2);
      TEPXClustersPerEvent_diskcombined[d]-> SetMarkerColor(2);
      TEPXClustersPerEvent_diskcombined[d]-> SetName(TString("clusters_combinedD")+ disklist[d].c_str());
      f->WriteTObject(TEPXClustersPerEvent_diskcombined[d]);
      
    }
    
    for (int r = 0; r < 5; r++){
      
      TEPXClustersPerEvent1[d][r]->SetMarkerSize(2);
      TEPXClustersPerEvent1[d][r]->SetMarkerColor(2);
      
      if (d <= 3) {
	TEPXClustersPerEvent1[d][r]->SetName(TString("clusters_") + "D" + (d-4) + "R" + (r+1));
      } else {
	
	TEPXClustersPerEvent1[d][r]->SetName(TString("clusters_") + "D" + (d-3) + "R" + (r+1));
      }
      
      f->WriteTObject(TEPXClustersPerEvent1[d][r]);
      
      if(disklist[d] > 0) {
	
	TEPXClustersPerEvent_diskcombinedring[d][r]->SetMarkerSize(2);
	TEPXClustersPerEvent_diskcombinedring[d][r]->SetMarkerColor(2);
	TEPXClustersPerEvent_diskcombinedring[d][r]->SetName(TString("clusters_combinedD")+ disklist[d].c_str()+"R"+(r+1));  
	f->WriteTObject(TEPXClustersPerEvent_diskcombinedring[d][r]);
	
      }
    }
  }
}







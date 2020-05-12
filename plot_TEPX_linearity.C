void plot_TEPX_linearity() {

    int option = 3;
    //TString outputpath = "/home/ashish/";
    TString outputpath = "/home/ashish/TEPX_plot1/";
    TString inpath = "/home/ashish/TEPX/samples_17Feb2020/";
    //gROOT->ProcessLine(".x BRIL-upgrade/rootlogon.C");
    
    //string containing names of input sample files
    std::vector<std::string> pulist = { "2023D42PU0p5", "2023D42PU1", "2023D42PU1p5", "2023D42PU2", "2023D42PU10", "2023D42PU30", "2023D42PU50", "2023D42PU100", "2023D42PU140", "2023D42PU200" };
    std::map<std::string, float> pumap = { { "2023D42PU0p5", 0.5 }, { "2023D42PU1", 1 }, { "2023D42PU1p5", 1.5 }, { "2023D42PU2", 2 }, { "2023D42PU10", 10 }, { "2023D42PU30", 30 }, { "2023D42PU50", 50 }, { "2023D42PU100", 100 }, { "2023D42PU140", 140 }, { "2023D42PU200", 200 } };
    std::vector<std::string> disklist = { "-4", "-3", "-2", "-1", "1", "2", "3", "4" };

    //create the profiles to be filled below.
    TGraphErrors* TEPXClustersPerEvent[20][4]; //number of clusters vs pu
    TH2F* Histogram2D[11][16];
    TH1D* ProjectionY[20][20];
    TH1D* ProjectionY1[20][20];
    TH1D* ProjectionY2[20][20];
    TH1D* ProjectionY3[20][20];
    TH1D* ProjectionY4[20][20];
    
    
  for (int pu = 0; pu < pulist.size(); pu++) {
    for (int d = 0; d < disklist.size(); d++){
      for (int r = 0; r < 6; r++){
            TEPXClustersPerEvent[d][r] = new TGraphErrors();
	    ProjectionY[pu][d] = new TH1D();
            ProjectionY1[pu][d] = new TH1D();
            ProjectionY2[pu][d] = new TH1D();
            ProjectionY3[pu][d] = new TH1D();
	    ProjectionY4[pu][d] = new TH1D();
      }
    }
}

    //Non-linearity graphs
    TGraphErrors* NonLinearity_TEPXClustersPerEvent[20][4]; //number of clusters vs pu
    for (int d = 0; d < disklist.size(); d++) {
        for (int r = 0; r < 5; r++) {
            NonLinearity_TEPXClustersPerEvent[d][r] = new TGraphErrors();
           
        }
    }
    
    

    //read the histograms
    TProfile* Prof_TEPXClustersPerEvent[20][20]; //array pu,disk
    for (int pu = 0; pu < pulist.size(); pu++) {
        TFile F(inpath + pulist[pu].c_str() + ".root", "read");
        gROOT->cd();

	
        for (int d = 0; d < disklist.size(); d++) {
            if (option == 1){

            //draw the 2D Cluster distributions for all PU
            TString histoname = "BRIL_IT_Analysis/TEPX/Clusters/Number of clusters for Disk ";

            TH2F* H = (TH2F*)F.Get(histoname + disklist[d].c_str());
	    Histogram2D[pu][d] = (TH2F*)H->Clone(TString(H->GetName()));

            ProjectionY[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 1, 1)->Clone(TString(H->GetName()) +             "Projection");
            ProjectionY1[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 2,2)->Clone(TString(H->GetName()) +             "Projection");
            ProjectionY2[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 3, 3)->Clone(TString(H->GetName()) +            "Projection");
	    ProjectionY3[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 4, 4)->Clone(TString(H->GetName()) +            "Projection");
	    ProjectionY4[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 5, 5)->Clone(TString(H->GetName()) +            "Projection");
            
            Prof_TEPXClustersPerEvent[pu][d] = (TProfile*)H->ProfileX()->Clone(TString(H->GetName()) + "Profile"); //number of cluster vs ring

		if (pu < 4){
	        TCanvas *c1 = new TCanvas();
            	c1->cd();
                Histogram2D[pu][d]->GetYaxis()->SetRangeUser(0, 350);
      		Histogram2D[pu][d]->Draw("colz");
      		char *hist2Dname = new char[10];
      		Int_t pud = pu * 8 + d + 1;
      		sprintf(hist2Dname,"histo2D_%d.gif",pud);
      		cout<<"=========================="<<hist2Dname<<endl;
      		c1->Print(outputpath + hist2Dname);
      		c1->Clear();
                
                TCanvas *c2 = new TCanvas();
            	c2->cd();
                ProjectionY[pu][d]->GetXaxis()->SetRangeUser(0, 60);
                ProjectionY[pu][d]->GetYaxis()->SetRangeUser(0, 70000);
                ProjectionY[pu][d]->SetLineColor(5);
      		ProjectionY[pu][d]->Draw("same");
                ProjectionY1[pu][d]->GetXaxis()->SetRangeUser(0, 60);
                ProjectionY1[pu][d]->GetYaxis()->SetRangeUser(0, 70000);
		ProjectionY1[pu][d]->SetLineColor(2);
      		ProjectionY1[pu][d]->Draw("same");
                ProjectionY2[pu][d]->GetXaxis()->SetRangeUser(0, 60);
                ProjectionY2[pu][d]->GetYaxis()->SetRangeUser(0, 70000);
		ProjectionY2[pu][d]->SetLineColor(3);
      		ProjectionY2[pu][d]->Draw("same");
		ProjectionY3[pu][d]->GetXaxis()->SetRangeUser(0, 60);
                ProjectionY3[pu][d]->GetYaxis()->SetRangeUser(0, 70000);
		ProjectionY3[pu][d]->SetLineColor(4);
      		ProjectionY3[pu][d]->Draw("same");
		ProjectionY4[pu][d]->GetXaxis()->SetRangeUser(0, 60);
                ProjectionY4[pu][d]->GetYaxis()->SetRangeUser(0, 70000);
		ProjectionY4[pu][d]->SetLineColor(6);
      		ProjectionY4[pu][d]->Draw("same");
                
               

                auto legend9 = new TLegend(0.7,0.7,0.88,0.88);
		legend9->AddEntry(ProjectionY[pu][d], "Ring 1" ,"l");
		legend9->AddEntry(ProjectionY1[pu][d], "Ring 2" ,"l");
		legend9->AddEntry(ProjectionY2[pu][d], "Ring 3" ,"l");
		legend9->AddEntry(ProjectionY3[pu][d], "Ring 4" ,"l");
                legend9->AddEntry(ProjectionY4[pu][d], "Ring 5" ,"l");

                legend9->SetFillColor(0);
  		legend9->SetLineColor(0);
  		legend9->SetFillColor(0);
		legend9->Draw();
      		char *Projection = new char[12];
      		Int_t pro = pu * 8 + d + 1;
      		sprintf(Projection,"Projection_%d.gif",pro);
      		cout<<"=========================="<<Projection<<endl;
      		c2->Print(outputpath + Projection);
      		c2->Clear();
		}

                if (pu >= 4){

	        TCanvas *c3 = new TCanvas();
            	c3->cd();
                Histogram2D[pu][d]->GetYaxis()->SetRangeUser(0, 3500);
      		Histogram2D[pu][d]->Draw("colz");
      		char *hist2Dname = new char[10];
      		Int_t pud1 = pu * 8 + d + 1;
      		sprintf(hist2Dname,"histo2D_%d.gif",pud1);
      		cout<<"=========================="<<hist2Dname<<endl;
      		c3->Print(outputpath + hist2Dname);
      		c3->Clear();

                
                TCanvas *c4 = new TCanvas();
            	c4->cd();
                ProjectionY[pu][d]->GetXaxis()->SetRangeUser(0, 1500);
                ProjectionY[pu][d]->GetYaxis()->SetRangeUser(0, 500);
                ProjectionY[pu][d]->SetLineColor(5);
      		ProjectionY[pu][d]->Draw("same");
                ProjectionY1[pu][d]->GetXaxis()->SetRangeUser(0, 1500);
                ProjectionY1[pu][d]->GetYaxis()->SetRangeUser(0, 500);
		ProjectionY1[pu][d]->SetLineColor(2);
      		ProjectionY1[pu][d]->Draw("same");
                ProjectionY2[pu][d]->GetXaxis()->SetRangeUser(0, 1500);
                ProjectionY2[pu][d]->GetYaxis()->SetRangeUser(0, 500);
		ProjectionY2[pu][d]->SetLineColor(3);
      		ProjectionY2[pu][d]->Draw("same");
		ProjectionY3[pu][d]->GetXaxis()->SetRangeUser(0, 1500);
                ProjectionY3[pu][d]->GetYaxis()->SetRangeUser(0, 500);
		ProjectionY3[pu][d]->SetLineColor(4);
      		ProjectionY3[pu][d]->Draw("same");
		ProjectionY4[pu][d]->GetXaxis()->SetRangeUser(0, 1500);
                ProjectionY4[pu][d]->GetYaxis()->SetRangeUser(0, 500);
		ProjectionY4[pu][d]->SetLineColor(6);
      		ProjectionY4[pu][d]->Draw("same");

                auto legend10 = new TLegend(0.7,0.7,0.88,0.88);
		legend10->AddEntry(ProjectionY[pu][d], "Ring 1" ,"l");
		legend10->AddEntry(ProjectionY1[pu][d], "Ring 2" ,"l");
		legend10->AddEntry(ProjectionY2[pu][d], "Ring 3" ,"l");
		legend10->AddEntry(ProjectionY3[pu][d], "Ring 4" ,"l");
                legend10->AddEntry(ProjectionY4[pu][d], "Ring 5" ,"l");

                legend10->SetFillColor(0);
  		legend10->SetLineColor(0);
  		legend10->SetFillColor(0);
		legend10->Draw();
      		char *Projection = new char[12];
      		Int_t pro = pu * 8 + d + 1;
      		sprintf(Projection,"Projection_%d.gif",pro);
      		cout<<"=========================="<<Projection<<endl;
      		c4->Print(outputpath + Projection);
      		c4->Clear();
	      
	}

                
                if(pu < 4){

                TCanvas *c6 = new TCanvas();
            	c6->cd();
                c6->SetLogy();

                ProjectionY[0][d]->SetMinimum(0.011);
                ProjectionY[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[0][d]->SetLineColor(1);
               	ProjectionY[0][d]->Draw("same");

		ProjectionY[1][d]->SetMinimum(0.001);
                ProjectionY[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[1][d]->SetLineColor(2);
               	ProjectionY[1][d]->Draw("same");
               
                ProjectionY[2][d]->SetMinimum(0.001);
                ProjectionY[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[2][d]->SetLineColor(3);
               	ProjectionY[2][d]->Draw("same");

                ProjectionY[3][d]->SetMinimum(0.01);
                ProjectionY[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[3][d]->SetLineColor(4);
               	ProjectionY[3][d]->Draw("same");

		

		auto legend1 = new TLegend(0.7,0.7,0.88,0.88);
		legend1->AddEntry(ProjectionY[0][d], "pu = 0.5" ,"l");
		legend1->AddEntry(ProjectionY[1][d], "pu = 1" ,"l");
		legend1->AddEntry(ProjectionY[2][d], "pu = 1.5" ,"l");
		legend1->AddEntry(ProjectionY[3][d], "pu = 2" ,"l");

                legend1->SetFillColor(0);
  		legend1->SetLineColor(0);
  		legend1->SetFillColor(0);
		legend1->Draw();

		char *Projection11 = new char[10];
      		Int_t pro12 = d + 1;
      		sprintf(Projection11,"ProjectionRing1_%d.gif",pro12);
      		cout<<"=========================="<<Projection11<<endl;
      		c6->Print(outputpath + Projection11);
      		c6->Clear();

                
		TCanvas *c7 = new TCanvas();
            	c7->cd();
                c7->SetLogy();

		ProjectionY1[0][d]->SetMinimum(0.011);
                ProjectionY1[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[0][d]->SetLineColor(1);
               	ProjectionY1[0][d]->Draw("same");

		ProjectionY1[1][d]->SetMinimum(0.001);
                ProjectionY1[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[1][d]->SetLineColor(2);
               	ProjectionY1[1][d]->Draw("same");
               
                ProjectionY1[2][d]->SetMinimum(0.001);
                ProjectionY1[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[2][d]->SetLineColor(3);
               	ProjectionY1[2][d]->Draw("same");

                ProjectionY1[3][d]->SetMinimum(0.01);
                ProjectionY1[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[3][d]->SetLineColor(4);
               	ProjectionY1[3][d]->Draw("same");

		auto legend2 = new TLegend(0.7,0.7,0.88,0.88);
		legend2->AddEntry(ProjectionY1[0][d], "pu = 0.5" ,"l");
		legend2->AddEntry(ProjectionY1[1][d], "pu = 1" ,"l");
		legend2->AddEntry(ProjectionY1[2][d], "pu = 1.5" ,"l");
		legend2->AddEntry(ProjectionY1[3][d], "pu = 2" ,"l");

                legend2->SetFillColor(0);
  		legend2->SetLineColor(0);
  		legend2->SetFillColor(0);
		legend2->Draw();

		char *Projection14 = new char[10];
      		Int_t pro14 = d + 1;
      		sprintf(Projection14,"ProjectionRing2_%d.gif",pro14);
      		cout<<"=========================="<<Projection14<<endl;
      		c7->Print(outputpath + Projection14);
      		c7->Clear();

		TCanvas *c8 = new TCanvas();
            	c8->cd();
                c8->SetLogy();

                ProjectionY2[0][d]->SetMinimum(0.011);
                ProjectionY2[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[0][d]->SetLineColor(1);
               	ProjectionY2[0][d]->Draw("same");

		ProjectionY2[1][d]->SetMinimum(0.001);
                ProjectionY2[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[1][d]->SetLineColor(2);
               	ProjectionY2[1][d]->Draw("same");
               
                ProjectionY2[2][d]->SetMinimum(0.001);
                ProjectionY2[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[2][d]->SetLineColor(3);
               	ProjectionY2[2][d]->Draw("same");

                ProjectionY2[3][d]->SetMinimum(0.01);
                ProjectionY2[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[3][d]->SetLineColor(4);
               	ProjectionY2[3][d]->Draw("same");

		auto legend3 = new TLegend(0.7,0.7,0.88,0.88);
		legend3->AddEntry(ProjectionY1[0][d], "pu = 0.5" ,"l");
		legend3->AddEntry(ProjectionY1[1][d], "pu = 1" ,"l");
		legend3->AddEntry(ProjectionY1[2][d], "pu = 1.5" ,"l");
		legend3->AddEntry(ProjectionY1[3][d], "pu = 2" ,"l");

                legend3->SetFillColor(0);
  		legend3->SetLineColor(0);
  		legend3->SetFillColor(0);
		legend3->Draw();

		char *Projection15 = new char[10];
      		Int_t pro15 = d + 1;
      		sprintf(Projection15,"ProjectionRing3_%d.gif",pro15);
      		cout<<"=========================="<<Projection15<<endl;
      		c8->Print(outputpath + Projection15);
      		c8->Clear();


		TCanvas *c9 = new TCanvas();
            	c9->cd();
                c9->SetLogy();

                ProjectionY3[0][d]->SetMinimum(0.011);
                ProjectionY3[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[0][d]->SetLineColor(1);
               	ProjectionY3[0][d]->Draw("same");

		ProjectionY3[1][d]->SetMinimum(0.001);
                ProjectionY3[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[1][d]->SetLineColor(2);
               	ProjectionY3[1][d]->Draw("same");
               
                ProjectionY3[2][d]->SetMinimum(0.001);
                ProjectionY3[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[2][d]->SetLineColor(3);
               	ProjectionY3[2][d]->Draw("same");

                ProjectionY3[3][d]->SetMinimum(0.01);
                ProjectionY3[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[3][d]->SetLineColor(4);
               	ProjectionY3[3][d]->Draw("same");

		auto legend4 = new TLegend(0.7,0.7,0.88,0.88);
		legend4->AddEntry(ProjectionY3[0][d], "pu = 0.5" ,"l");
		legend4->AddEntry(ProjectionY3[1][d], "pu = 1" ,"l");
		legend4->AddEntry(ProjectionY3[2][d], "pu = 1.5" ,"l");
		legend4->AddEntry(ProjectionY3[3][d], "pu = 2" ,"l");

                legend4->SetFillColor(0);
  		legend4->SetLineColor(0);
  		legend4->SetFillColor(0);
		legend4->Draw();

		char *Projection16 = new char[10];
      		Int_t pro16 = d + 1;
      		sprintf(Projection16,"ProjectionRing4_%d.gif",pro16);
      		cout<<"=========================="<<Projection16<<endl;
      		c8->Print(outputpath + Projection16);
      		c8->Clear();



		TCanvas *c10 = new TCanvas();
            	c10->cd();
                c10->SetLogy();

                ProjectionY4[0][d]->SetMinimum(0.011);
                ProjectionY4[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[0][d]->SetLineColor(1);
               	ProjectionY4[0][d]->Draw("same");

		ProjectionY4[1][d]->SetMinimum(0.001);
                ProjectionY4[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[1][d]->SetLineColor(2);
               	ProjectionY4[1][d]->Draw("same");
               
                ProjectionY4[2][d]->SetMinimum(0.001);
                ProjectionY4[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[2][d]->SetLineColor(3);
               	ProjectionY4[2][d]->Draw("same");

                ProjectionY4[3][d]->SetMinimum(0.01);
                ProjectionY4[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[3][d]->SetLineColor(4);
               	ProjectionY4[3][d]->Draw("same");

		auto legend5 = new TLegend(0.7,0.7,0.88,0.88);
		legend5->AddEntry(ProjectionY4[0][d], "pu = 0.5" ,"l");
		legend5->AddEntry(ProjectionY4[1][d], "pu = 1" ,"l");
		legend5->AddEntry(ProjectionY4[2][d], "pu = 1.5" ,"l");
		legend5->AddEntry(ProjectionY4[3][d], "pu = 2" ,"l");

                legend5->SetFillColor(0);
  		legend5->SetLineColor(0);
  		legend5->SetFillColor(0);
		legend5->Draw();

		char *Projection17 = new char[10];
      		Int_t pro17 = d + 1;
      		sprintf(Projection17,"ProjectionRing5_%d.gif",pro17);
      		cout<<"=========================="<<Projection17<<endl;
      		c10->Print(outputpath + Projection17);
      		c10->Clear();

		}

                
		TCanvas *c7 = new TCanvas();
            	c7->cd();
                c7->SetLogy();

                ProjectionY[0][d]->SetMinimum(0.011);
                ProjectionY[0][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[0][d]->SetLineColor(1);
               	ProjectionY[0][d]->Draw("same");

		ProjectionY[1][d]->SetMinimum(0.001);
                ProjectionY[1][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[1][d]->SetLineColor(2);
               	ProjectionY[1][d]->Draw("same");
               
                ProjectionY[2][d]->SetMinimum(0.001);
                ProjectionY[2][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[2][d]->SetLineColor(3);
               	ProjectionY[2][d]->Draw("same");

                ProjectionY[3][d]->SetMinimum(0.01);
                ProjectionY[3][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[3][d]->SetLineColor(4);
               	ProjectionY[3][d]->Draw("same");

                ProjectionY[4][d]->SetMinimum(0.01);
                ProjectionY[4][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[4][d]->SetLineColor(5);
               	ProjectionY[4][d]->Draw("same");

		ProjectionY[5][d]->SetMinimum(0.01);
                ProjectionY[5][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[5][d]->SetLineColor(6);
               	ProjectionY[5][d]->Draw("same");

		ProjectionY[6][d]->SetMinimum(0.01);
                ProjectionY[6][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[6][d]->SetLineColor(7);
               	ProjectionY[6][d]->Draw("same");

		ProjectionY[7][d]->SetMinimum(0.01);
                ProjectionY[7][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[7][d]->SetLineColor(8);
               	ProjectionY[7][d]->Draw("same");

		ProjectionY[8][d]->SetMinimum(0.01);
                ProjectionY[8][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[8][d]->SetLineColor(9);
               	ProjectionY[8][d]->Draw("same");

		ProjectionY[9][d]->SetMinimum(0.01);
                ProjectionY[9][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[9][d]->SetLineColor(44);
               	ProjectionY[9][d]->Draw("same");

 

		auto legend12 = new TLegend(0.7,0.7,0.88,0.88);
		legend12->AddEntry(ProjectionY[0][d], "pu = 0.5" ,"l");
		legend12->AddEntry(ProjectionY[1][d], "pu = 1" ,"l");
		legend12->AddEntry(ProjectionY[2][d], "pu = 1.5" ,"l");
		legend12->AddEntry(ProjectionY[3][d], "pu = 2" ,"l");
		legend12->AddEntry(ProjectionY[4][d], "pu = 10" ,"l");
		legend12->AddEntry(ProjectionY[5][d], "pu = 30" ,"l");
		legend12->AddEntry(ProjectionY[6][d], "pu = 50" ,"l");
		legend12->AddEntry(ProjectionY[7][d], "pu = 100" ,"l");
		legend12->AddEntry(ProjectionY[8][d], "pu = 140" ,"l");
		legend12->AddEntry(ProjectionY[9][d], "pu = 200" ,"l");
                legend12->SetFillColor(0);
  		legend12->SetLineColor(0);
  		legend12->SetFillColor(0);
		legend12->Draw();


		char *Projection12 = new char[10];
      		Int_t pro13 = d + 1;
      		sprintf(Projection12,"ProjectionsRing1_%d.gif",pro13);
      		cout<<"=========================="<<Projection12<<endl;
      		c7->Print(outputpath + Projection12);
      		c7->Clear();





		TCanvas *c8 = new TCanvas();
            	c8->cd();
                c8->SetLogy();

                ProjectionY1[0][d]->SetMinimum(0.011);
                ProjectionY1[0][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[0][d]->SetLineColor(1);
               	ProjectionY1[0][d]->Draw("same");

		ProjectionY1[1][d]->SetMinimum(0.001);
                ProjectionY1[1][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[1][d]->SetLineColor(2);
               	ProjectionY1[1][d]->Draw("same");
               
                ProjectionY1[2][d]->SetMinimum(0.001);
                ProjectionY1[2][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[2][d]->SetLineColor(3);
               	ProjectionY1[2][d]->Draw("same");

                ProjectionY1[3][d]->SetMinimum(0.01);
                ProjectionY1[3][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[3][d]->SetLineColor(4);
               	ProjectionY1[3][d]->Draw("same");

                ProjectionY1[4][d]->SetMinimum(0.01);
                ProjectionY1[4][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[4][d]->SetLineColor(5);
               	ProjectionY1[4][d]->Draw("same");

		ProjectionY1[5][d]->SetMinimum(0.01);
                ProjectionY1[5][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[5][d]->SetLineColor(6);
               	ProjectionY1[5][d]->Draw("same");

		ProjectionY1[6][d]->SetMinimum(0.01);
                ProjectionY1[6][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[6][d]->SetLineColor(7);
               	ProjectionY1[6][d]->Draw("same");

		ProjectionY1[7][d]->SetMinimum(0.01);
                ProjectionY1[7][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[7][d]->SetLineColor(8);
               	ProjectionY1[7][d]->Draw("same");

		ProjectionY1[8][d]->SetMinimum(0.01);
                ProjectionY1[8][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[8][d]->SetLineColor(9);
               	ProjectionY1[8][d]->Draw("same");

		ProjectionY1[9][d]->SetMinimum(0.01);
                ProjectionY1[9][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[9][d]->SetLineColor(44);
               	ProjectionY1[9][d]->Draw("same");

 

		auto legend13 = new TLegend(0.7,0.7,0.88,0.88);
		legend13->AddEntry(ProjectionY[0][d], "pu = 0.5" ,"l");
		legend13->AddEntry(ProjectionY[1][d], "pu = 1" ,"l");
		legend13->AddEntry(ProjectionY[2][d], "pu = 1.5" ,"l");
		legend13->AddEntry(ProjectionY[3][d], "pu = 2" ,"l");
		legend13->AddEntry(ProjectionY[4][d], "pu = 10" ,"l");
		legend13->AddEntry(ProjectionY[5][d], "pu = 30" ,"l");
		legend13->AddEntry(ProjectionY[6][d], "pu = 50" ,"l");
		legend13->AddEntry(ProjectionY[7][d], "pu = 100" ,"l");
		legend13->AddEntry(ProjectionY[8][d], "pu = 140" ,"l");
		legend13->AddEntry(ProjectionY[9][d], "pu = 200" ,"l");
                legend13->SetFillColor(0);
  		legend13->SetLineColor(0);
  		legend13->SetFillColor(0);
		legend13->Draw();


		char *Projection14 = new char[10];
      		Int_t pro14 = d + 1;
      		sprintf(Projection14,"ProjectionsRing2_%d.gif",pro14);
      		cout<<"=========================="<<Projection14<<endl;
      		c8->Print(outputpath + Projection14);
      		c8->Clear();



		TCanvas *c9 = new TCanvas();
            	c9->cd();
                c9->SetLogy();

                ProjectionY2[0][d]->SetMinimum(0.011);
                ProjectionY2[0][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[0][d]->SetLineColor(1);
               	ProjectionY2[0][d]->Draw("same");

		ProjectionY2[1][d]->SetMinimum(0.001);
                ProjectionY2[1][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[1][d]->SetLineColor(2);
               	ProjectionY2[1][d]->Draw("same");
               
                ProjectionY2[2][d]->SetMinimum(0.001);
                ProjectionY2[2][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[2][d]->SetLineColor(3);
               	ProjectionY2[2][d]->Draw("same");

                ProjectionY2[3][d]->SetMinimum(0.01);
                ProjectionY2[3][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[3][d]->SetLineColor(4);
               	ProjectionY2[3][d]->Draw("same");

                ProjectionY2[4][d]->SetMinimum(0.01);
                ProjectionY2[4][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[4][d]->SetLineColor(5);
               	ProjectionY2[4][d]->Draw("same");

		ProjectionY2[5][d]->SetMinimum(0.01);
                ProjectionY2[5][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[5][d]->SetLineColor(6);
               	ProjectionY2[5][d]->Draw("same");

		ProjectionY2[6][d]->SetMinimum(0.01);
                ProjectionY2[6][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[6][d]->SetLineColor(7);
               	ProjectionY2[6][d]->Draw("same");

		ProjectionY2[7][d]->SetMinimum(0.01);
                ProjectionY2[7][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[7][d]->SetLineColor(8);
               	ProjectionY2[7][d]->Draw("same");

		ProjectionY2[8][d]->SetMinimum(0.01);
                ProjectionY2[8][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[8][d]->SetLineColor(9);
               	ProjectionY2[8][d]->Draw("same");

		ProjectionY2[9][d]->SetMinimum(0.01);
                ProjectionY2[9][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[9][d]->SetLineColor(44);
               	ProjectionY2[9][d]->Draw("same");

 

		auto legend14 = new TLegend(0.7,0.7,0.88,0.88);
		legend14->AddEntry(ProjectionY[0][d], "pu = 0.5" ,"l");
		legend14->AddEntry(ProjectionY[1][d], "pu = 1" ,"l");
		legend14->AddEntry(ProjectionY[2][d], "pu = 1.5" ,"l");
		legend14->AddEntry(ProjectionY[3][d], "pu = 2" ,"l");
		legend14->AddEntry(ProjectionY[4][d], "pu = 10" ,"l");
		legend14->AddEntry(ProjectionY[5][d], "pu = 30" ,"l");
		legend14->AddEntry(ProjectionY[6][d], "pu = 50" ,"l");
		legend14->AddEntry(ProjectionY[7][d], "pu = 100" ,"l");
		legend14->AddEntry(ProjectionY[8][d], "pu = 140" ,"l");
		legend14->AddEntry(ProjectionY[9][d], "pu = 200" ,"l");
                legend14->SetFillColor(0);
  		legend14->SetLineColor(0);
  		legend14->SetFillColor(0);
		legend14->Draw();


		char *Projection15 = new char[10];
      		Int_t pro15 = d + 1;
      		sprintf(Projection15,"ProjectionsRing3_%d.gif",pro15);
      		cout<<"=========================="<<Projection15<<endl;
      		c9->Print(outputpath + Projection15);
      		c9->Clear();



		TCanvas *c10 = new TCanvas();
            	c10->cd();
                c10->SetLogy();

                ProjectionY3[0][d]->SetMinimum(0.011);
                ProjectionY3[0][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[0][d]->SetLineColor(1);
               	ProjectionY3[0][d]->Draw("same");

		ProjectionY3[1][d]->SetMinimum(0.001);
                ProjectionY3[1][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[1][d]->SetLineColor(2);
               	ProjectionY3[1][d]->Draw("same");
               
                ProjectionY3[2][d]->SetMinimum(0.001);
                ProjectionY3[2][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[2][d]->SetLineColor(3);
               	ProjectionY3[2][d]->Draw("same");

                ProjectionY3[3][d]->SetMinimum(0.01);
                ProjectionY3[3][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[3][d]->SetLineColor(4);
               	ProjectionY3[3][d]->Draw("same");

                ProjectionY3[4][d]->SetMinimum(0.01);
                ProjectionY3[4][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[4][d]->SetLineColor(5);
               	ProjectionY3[4][d]->Draw("same");

		ProjectionY3[5][d]->SetMinimum(0.01);
                ProjectionY3[5][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[5][d]->SetLineColor(6);
               	ProjectionY3[5][d]->Draw("same");

		ProjectionY3[6][d]->SetMinimum(0.01);
                ProjectionY3[6][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[6][d]->SetLineColor(7);
               	ProjectionY3[6][d]->Draw("same");

		ProjectionY3[7][d]->SetMinimum(0.01);
                ProjectionY3[7][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[7][d]->SetLineColor(8);
               	ProjectionY3[7][d]->Draw("same");

		ProjectionY3[8][d]->SetMinimum(0.01);
                ProjectionY3[8][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[8][d]->SetLineColor(9);
               	ProjectionY3[8][d]->Draw("same");

		ProjectionY3[9][d]->SetMinimum(0.01);
                ProjectionY3[9][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[9][d]->SetLineColor(44);
               	ProjectionY3[9][d]->Draw("same");

 

		auto legend15 = new TLegend(0.7,0.7,0.88,0.88);
		legend15->AddEntry(ProjectionY3[0][d], "pu = 0.5" ,"l");
		legend15->AddEntry(ProjectionY3[1][d], "pu = 1" ,"l");
		legend15->AddEntry(ProjectionY3[2][d], "pu = 1.5" ,"l");
		legend15->AddEntry(ProjectionY3[3][d], "pu = 2" ,"l");
		legend15->AddEntry(ProjectionY3[4][d], "pu = 10" ,"l");
		legend15->AddEntry(ProjectionY3[5][d], "pu = 30" ,"l");
		legend15->AddEntry(ProjectionY3[6][d], "pu = 50" ,"l");
		legend15->AddEntry(ProjectionY3[7][d], "pu = 100" ,"l");
		legend15->AddEntry(ProjectionY3[8][d], "pu = 140" ,"l");
		legend15->AddEntry(ProjectionY3[9][d], "pu = 200" ,"l");
                legend15->SetFillColor(0);
  		legend15->SetLineColor(0);
  		legend15->SetFillColor(0);
		legend15->Draw();

		char *Projection16 = new char[10];
      		Int_t pro16 = d + 1;
      		sprintf(Projection16,"ProjectionsRing4_%d.gif",pro16);
      		cout<<"=========================="<<Projection16<<endl;
      		c10->Print(outputpath + Projection16);
      		c10->Clear();



		TCanvas *c11 = new TCanvas();
            	c11->cd();
                c11->SetLogy();

                ProjectionY4[0][d]->SetMinimum(0.011);
                ProjectionY4[0][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[0][d]->SetLineColor(1);
               	ProjectionY4[0][d]->Draw("same");

		ProjectionY4[1][d]->SetMinimum(0.001);
                ProjectionY4[1][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[1][d]->SetLineColor(2);
               	ProjectionY4[1][d]->Draw("same");
               
                ProjectionY4[2][d]->SetMinimum(0.001);
                ProjectionY4[2][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[2][d]->SetLineColor(3);
               	ProjectionY4[2][d]->Draw("same");

                ProjectionY4[3][d]->SetMinimum(0.01);
                ProjectionY4[3][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[3][d]->SetLineColor(4);
               	ProjectionY4[3][d]->Draw("same");

                ProjectionY4[4][d]->SetMinimum(0.01);
                ProjectionY4[4][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[4][d]->SetLineColor(5);
               	ProjectionY4[4][d]->Draw("same");

		ProjectionY4[5][d]->SetMinimum(0.01);
                ProjectionY4[5][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[5][d]->SetLineColor(6);
               	ProjectionY4[5][d]->Draw("same");

		ProjectionY4[6][d]->SetMinimum(0.01);
                ProjectionY4[6][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[6][d]->SetLineColor(7);
               	ProjectionY4[6][d]->Draw("same");

		ProjectionY4[7][d]->SetMinimum(0.01);
                ProjectionY4[7][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[7][d]->SetLineColor(8);
               	ProjectionY4[7][d]->Draw("same");

		ProjectionY4[8][d]->SetMinimum(0.01);
                ProjectionY4[8][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[8][d]->SetLineColor(9);
               	ProjectionY4[8][d]->Draw("same");

		ProjectionY4[9][d]->SetMinimum(0.01);
                ProjectionY4[9][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[9][d]->SetLineColor(44);
               	ProjectionY4[9][d]->Draw("same");

 

		auto legend16 = new TLegend(0.7,0.7,0.88,0.88);
		legend16->AddEntry(ProjectionY4[0][d], "pu = 0.5" ,"l");
		legend16->AddEntry(ProjectionY4[1][d], "pu = 1" ,"l");
		legend16->AddEntry(ProjectionY4[2][d], "pu = 1.5" ,"l");
		legend16->AddEntry(ProjectionY4[3][d], "pu = 2" ,"l");
		legend16->AddEntry(ProjectionY4[4][d], "pu = 10" ,"l");
		legend16->AddEntry(ProjectionY4[5][d], "pu = 30" ,"l");
		legend16->AddEntry(ProjectionY4[6][d], "pu = 50" ,"l");
		legend16->AddEntry(ProjectionY4[7][d], "pu = 100" ,"l");
		legend16->AddEntry(ProjectionY4[8][d], "pu = 140" ,"l");
		legend16->AddEntry(ProjectionY4[9][d], "pu = 200" ,"l");
                legend16->SetFillColor(0);
  		legend16->SetLineColor(0);
  		legend16->SetFillColor(0);
		legend16->Draw();


		char *Projection17 = new char[10];
      		Int_t pro17 = d + 1;
      		sprintf(Projection17,"ProjectionsRing5_%d.gif",pro17);
      		cout<<"=========================="<<Projection17<<endl;
      		c11->Print(outputpath + Projection17);
      		c11->Clear();

               
	       
 }

	    
           
            if (option == 2){

             //draw the 2D 2x Coincidences distributions for all PU
            TString histoname = "BRIL_IT_Analysis/TEPX/2xCoincidences/Number of 2x Coincidences for Disk ";

            TH2F* H = (TH2F*)F.Get(histoname + disklist[d].c_str());
	    Histogram2D[pu][d] = (TH2F*)H->Clone(TString(H->GetName()));

            ProjectionY[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 1, 1)->Clone(TString(H->GetName()) +             "Projection");
            ProjectionY1[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 2,2)->Clone(TString(H->GetName()) +             "Projection");
            ProjectionY2[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 3, 3)->Clone(TString(H->GetName()) +            "Projection");
	    ProjectionY3[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 4, 4)->Clone(TString(H->GetName()) +            "Projection");
	    ProjectionY4[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 5, 5)->Clone(TString(H->GetName()) +            "Projection");
            
            Prof_TEPXClustersPerEvent[pu][d] = (TProfile*)H->ProfileX()->Clone(TString(H->GetName()) + "Profile"); //number of cluster vs ring

		if (pu < 4){
	        TCanvas *c1 = new TCanvas();
            	c1->cd();
                Histogram2D[pu][d]->GetYaxis()->SetRangeUser(0, 100);
      		Histogram2D[pu][d]->Draw("colz");
      		char *hist2Dname = new char[10];
      		Int_t pud = pu * 8 + d + 1;
      		sprintf(hist2Dname,"histos2D_%d.gif",pud);
      		cout<<"=========================="<<hist2Dname<<endl;
      		c1->Print(outputpath + hist2Dname);
      		c1->Clear();
                
                TCanvas *c2 = new TCanvas();
            	c2->cd();
                ProjectionY[pu][d]->GetXaxis()->SetRangeUser(0, 20);
                ProjectionY[pu][d]->GetYaxis()->SetRangeUser(0, 20000);
                ProjectionY[pu][d]->SetLineColor(5);
      		ProjectionY[pu][d]->Draw("same");
                ProjectionY1[pu][d]->GetXaxis()->SetRangeUser(0, 20);
                ProjectionY1[pu][d]->GetYaxis()->SetRangeUser(0, 20000);
		ProjectionY1[pu][d]->SetLineColor(2);
      		ProjectionY1[pu][d]->Draw("same");
                ProjectionY2[pu][d]->GetXaxis()->SetRangeUser(0, 20);
                ProjectionY2[pu][d]->GetYaxis()->SetRangeUser(0, 20000);
		ProjectionY2[pu][d]->SetLineColor(3);
      		ProjectionY2[pu][d]->Draw("same");
		ProjectionY3[pu][d]->GetXaxis()->SetRangeUser(0, 20);
                ProjectionY3[pu][d]->GetYaxis()->SetRangeUser(0, 20000);
		ProjectionY3[pu][d]->SetLineColor(4);
      		ProjectionY3[pu][d]->Draw("same");
		ProjectionY4[pu][d]->GetXaxis()->SetRangeUser(0, 20);
                ProjectionY4[pu][d]->GetYaxis()->SetRangeUser(0, 20000);
		ProjectionY4[pu][d]->SetLineColor(6);
      		ProjectionY4[pu][d]->Draw("same");
                
               

                auto legend10 = new TLegend(0.7,0.7,0.88,0.88);
		legend10->AddEntry(ProjectionY[pu][d], "Ring 1" ,"l");
		legend10->AddEntry(ProjectionY1[pu][d], "Ring 2" ,"l");
		legend10->AddEntry(ProjectionY2[pu][d], "Ring 3" ,"l");
		legend10->AddEntry(ProjectionY3[pu][d], "Ring 4" ,"l");
                legend10->AddEntry(ProjectionY4[pu][d], "Ring 5" ,"l");

                legend10->SetFillColor(0);
  		legend10->SetLineColor(0);
  		legend10->SetFillColor(0);
		legend10->Draw();
      		char *Projection = new char[12];
      		Int_t pro = pu * 8 + d + 1;
      		sprintf(Projection,"Projection2x_%d.gif",pro);
      		cout<<"=========================="<<Projection<<endl;
      		c2->Print(outputpath + Projection);
      		c2->Clear();

}
                if (pu >= 4){

	        TCanvas *c3 = new TCanvas();
            	c3->cd();
                Histogram2D[pu][d]->GetYaxis()->SetRangeUser(0, 1000);
      		Histogram2D[pu][d]->Draw("colz");
      		char *hist2Dname = new char[10];
      		Int_t pud1 = pu * 8 + d + 1;
      		sprintf(hist2Dname,"histos2D_%d.gif",pud1);
      		cout<<"=========================="<<hist2Dname<<endl;
      		c3->Print(outputpath + hist2Dname);
      		c3->Clear();

                
                TCanvas *c4 = new TCanvas();
            	c4->cd();
                ProjectionY[pu][d]->GetXaxis()->SetRangeUser(0, 500);
                ProjectionY[pu][d]->GetYaxis()->SetRangeUser(0, 2500);
                ProjectionY[pu][d]->SetLineColor(5);
      		ProjectionY[pu][d]->Draw("same");
                ProjectionY1[pu][d]->GetXaxis()->SetRangeUser(0, 500);
                ProjectionY1[pu][d]->GetYaxis()->SetRangeUser(0, 2500);
		ProjectionY1[pu][d]->SetLineColor(2);
      		ProjectionY1[pu][d]->Draw("same");
                ProjectionY2[pu][d]->GetXaxis()->SetRangeUser(0, 500);
                ProjectionY2[pu][d]->GetYaxis()->SetRangeUser(0, 2500);
		ProjectionY2[pu][d]->SetLineColor(3);
      		ProjectionY2[pu][d]->Draw("same");
		ProjectionY3[pu][d]->GetXaxis()->SetRangeUser(0, 500);
                ProjectionY3[pu][d]->GetYaxis()->SetRangeUser(0, 2500);
		ProjectionY3[pu][d]->SetLineColor(4);
      		ProjectionY3[pu][d]->Draw("same");
		ProjectionY4[pu][d]->GetXaxis()->SetRangeUser(0, 500);
                ProjectionY4[pu][d]->GetYaxis()->SetRangeUser(0, 2500);
		ProjectionY4[pu][d]->SetLineColor(6);
      		ProjectionY4[pu][d]->Draw("same");

                auto legend10 = new TLegend(0.7,0.7,0.88,0.88);
		legend10->AddEntry(ProjectionY[pu][d], "Ring 1" ,"l");
		legend10->AddEntry(ProjectionY1[pu][d], "Ring 2" ,"l");
		legend10->AddEntry(ProjectionY2[pu][d], "Ring 3" ,"l");
		legend10->AddEntry(ProjectionY3[pu][d], "Ring 4" ,"l");
                legend10->AddEntry(ProjectionY4[pu][d], "Ring 5" ,"l");

                legend10->SetFillColor(0);
  		legend10->SetLineColor(0);
  		legend10->SetFillColor(0);
		legend10->Draw();
      		char *Projection = new char[12];
      		Int_t pro = pu * 8 + d + 1;
      		sprintf(Projection,"Projection2x_%d.gif",pro);
      		cout<<"=========================="<<Projection<<endl;
      		c4->Print(outputpath + Projection);
      		c4->Clear();
	      
	}

                
                if(pu < 4){

                TCanvas *c6 = new TCanvas();
            	c6->cd();
                c6->SetLogy();

                ProjectionY[0][d]->SetMinimum(0.011);
                ProjectionY[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[0][d]->SetLineColor(1);
               	ProjectionY[0][d]->Draw("same");

		ProjectionY[1][d]->SetMinimum(0.001);
                ProjectionY[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[1][d]->SetLineColor(2);
               	ProjectionY[1][d]->Draw("same");
               
                ProjectionY[2][d]->SetMinimum(0.001);
                ProjectionY[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[2][d]->SetLineColor(3);
               	ProjectionY[2][d]->Draw("same");

                ProjectionY[3][d]->SetMinimum(0.01);
                ProjectionY[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[3][d]->SetLineColor(4);
               	ProjectionY[3][d]->Draw("same");

		

		auto legend1 = new TLegend(0.7,0.7,0.88,0.88);
		legend1->AddEntry(ProjectionY[0][d], "pu = 0.5" ,"l");
		legend1->AddEntry(ProjectionY[1][d], "pu = 1" ,"l");
		legend1->AddEntry(ProjectionY[2][d], "pu = 1.5" ,"l");
		legend1->AddEntry(ProjectionY[3][d], "pu = 2" ,"l");

                legend1->SetFillColor(0);
  		legend1->SetLineColor(0);
  		legend1->SetFillColor(0);
		legend1->Draw();

		char *Projection11 = new char[10];
      		Int_t pro12 = d + 1;
      		sprintf(Projection11,"ProjectionasRing1_%d.gif",pro12);
      		cout<<"=========================="<<Projection11<<endl;
      		c6->Print(outputpath + Projection11);
      		c6->Clear();

                
		TCanvas *c7 = new TCanvas();
            	c7->cd();
                c7->SetLogy();

		ProjectionY1[0][d]->SetMinimum(0.011);
                ProjectionY1[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[0][d]->SetLineColor(1);
               	ProjectionY1[0][d]->Draw("same");

		ProjectionY1[1][d]->SetMinimum(0.001);
                ProjectionY1[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[1][d]->SetLineColor(2);
               	ProjectionY1[1][d]->Draw("same");
               
                ProjectionY1[2][d]->SetMinimum(0.001);
                ProjectionY1[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[2][d]->SetLineColor(3);
               	ProjectionY1[2][d]->Draw("same");

                ProjectionY1[3][d]->SetMinimum(0.01);
                ProjectionY1[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[3][d]->SetLineColor(4);
               	ProjectionY1[3][d]->Draw("same");

		auto legend2 = new TLegend(0.7,0.7,0.88,0.88);
		legend2->AddEntry(ProjectionY1[0][d], "pu = 0.5" ,"l");
		legend2->AddEntry(ProjectionY1[1][d], "pu = 1" ,"l");
		legend2->AddEntry(ProjectionY1[2][d], "pu = 1.5" ,"l");
		legend2->AddEntry(ProjectionY1[3][d], "pu = 2" ,"l");

                legend2->SetFillColor(0);
  		legend2->SetLineColor(0);
  		legend2->SetFillColor(0);
		legend2->Draw();

		char *Projection14 = new char[10];
      		Int_t pro14 = d + 1;
      		sprintf(Projection14,"ProjectionasRing2_%d.gif",pro14);
      		cout<<"=========================="<<Projection14<<endl;
      		c7->Print(outputpath + Projection14);
      		c7->Clear();

		TCanvas *c8 = new TCanvas();
            	c8->cd();
                c8->SetLogy();

                ProjectionY2[0][d]->SetMinimum(0.011);
                ProjectionY2[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[0][d]->SetLineColor(1);
               	ProjectionY2[0][d]->Draw("same");

		ProjectionY2[1][d]->SetMinimum(0.001);
                ProjectionY2[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[1][d]->SetLineColor(2);
               	ProjectionY2[1][d]->Draw("same");
               
                ProjectionY2[2][d]->SetMinimum(0.001);
                ProjectionY2[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[2][d]->SetLineColor(3);
               	ProjectionY2[2][d]->Draw("same");

                ProjectionY2[3][d]->SetMinimum(0.01);
                ProjectionY2[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[3][d]->SetLineColor(4);
               	ProjectionY2[3][d]->Draw("same");

		auto legend3 = new TLegend(0.7,0.7,0.88,0.88);
		legend3->AddEntry(ProjectionY1[0][d], "pu = 0.5" ,"l");
		legend3->AddEntry(ProjectionY1[1][d], "pu = 1" ,"l");
		legend3->AddEntry(ProjectionY1[2][d], "pu = 1.5" ,"l");
		legend3->AddEntry(ProjectionY1[3][d], "pu = 2" ,"l");

                legend3->SetFillColor(0);
  		legend3->SetLineColor(0);
  		legend3->SetFillColor(0);
		legend3->Draw();

		char *Projection15 = new char[10];
      		Int_t pro15 = d + 1;
      		sprintf(Projection15,"ProjectionasRing3_%d.gif",pro15);
      		cout<<"=========================="<<Projection15<<endl;
      		c8->Print(outputpath + Projection15);
      		c8->Clear();


		TCanvas *c9 = new TCanvas();
            	c9->cd();
                c9->SetLogy();

                ProjectionY3[0][d]->SetMinimum(0.011);
                ProjectionY3[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[0][d]->SetLineColor(1);
               	ProjectionY3[0][d]->Draw("same");

		ProjectionY3[1][d]->SetMinimum(0.001);
                ProjectionY3[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[1][d]->SetLineColor(2);
               	ProjectionY3[1][d]->Draw("same");
               
                ProjectionY3[2][d]->SetMinimum(0.001);
                ProjectionY3[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[2][d]->SetLineColor(3);
               	ProjectionY3[2][d]->Draw("same");

                ProjectionY3[3][d]->SetMinimum(0.01);
                ProjectionY3[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[3][d]->SetLineColor(4);
               	ProjectionY3[3][d]->Draw("same");

		auto legend4 = new TLegend(0.7,0.7,0.88,0.88);
		legend4->AddEntry(ProjectionY3[0][d], "pu = 0.5" ,"l");
		legend4->AddEntry(ProjectionY3[1][d], "pu = 1" ,"l");
		legend4->AddEntry(ProjectionY3[2][d], "pu = 1.5" ,"l");
		legend4->AddEntry(ProjectionY3[3][d], "pu = 2" ,"l");

                legend4->SetFillColor(0);
  		legend4->SetLineColor(0);
  		legend4->SetFillColor(0);
		legend4->Draw();

		char *Projection16 = new char[10];
      		Int_t pro16 = d + 1;
      		sprintf(Projection16,"ProjectionasRing4_%d.gif",pro16);
      		cout<<"=========================="<<Projection16<<endl;
      		c8->Print(outputpath + Projection16);
      		c8->Clear();



		TCanvas *c10 = new TCanvas();
            	c10->cd();
                c10->SetLogy();

                ProjectionY4[0][d]->SetMinimum(0.011);
                ProjectionY4[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[0][d]->SetLineColor(1);
               	ProjectionY4[0][d]->Draw("same");

		ProjectionY4[1][d]->SetMinimum(0.001);
                ProjectionY4[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[1][d]->SetLineColor(2);
               	ProjectionY4[1][d]->Draw("same");
               
                ProjectionY4[2][d]->SetMinimum(0.001);
                ProjectionY4[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[2][d]->SetLineColor(3);
               	ProjectionY4[2][d]->Draw("same");

                ProjectionY4[3][d]->SetMinimum(0.01);
                ProjectionY4[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[3][d]->SetLineColor(4);
               	ProjectionY4[3][d]->Draw("same");

		auto legend5 = new TLegend(0.7,0.7,0.88,0.88);
		legend5->AddEntry(ProjectionY4[0][d], "pu = 0.5" ,"l");
		legend5->AddEntry(ProjectionY4[1][d], "pu = 1" ,"l");
		legend5->AddEntry(ProjectionY4[2][d], "pu = 1.5" ,"l");
		legend5->AddEntry(ProjectionY4[3][d], "pu = 2" ,"l");

                legend5->SetFillColor(0);
  		legend5->SetLineColor(0);
  		legend5->SetFillColor(0);
		legend5->Draw();

		char *Projection17 = new char[10];
      		Int_t pro17 = d + 1;
      		sprintf(Projection17,"ProjectionasRing5_%d.gif",pro17);
      		cout<<"=========================="<<Projection17<<endl;
      		c10->Print(outputpath + Projection17);
      		c10->Clear();

		}

                
		TCanvas *c7 = new TCanvas();
            	c7->cd();
                c7->SetLogy();

                ProjectionY[0][d]->SetMinimum(0.011);
                ProjectionY[0][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[0][d]->SetLineColor(1);
               	ProjectionY[0][d]->Draw("same");

		ProjectionY[1][d]->SetMinimum(0.001);
                ProjectionY[1][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[1][d]->SetLineColor(2);
               	ProjectionY[1][d]->Draw("same");
               
                ProjectionY[2][d]->SetMinimum(0.001);
                ProjectionY[2][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[2][d]->SetLineColor(3);
               	ProjectionY[2][d]->Draw("same");

                ProjectionY[3][d]->SetMinimum(0.01);
                ProjectionY[3][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[3][d]->SetLineColor(4);
               	ProjectionY[3][d]->Draw("same");

                ProjectionY[4][d]->SetMinimum(0.01);
                ProjectionY[4][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[4][d]->SetLineColor(5);
               	ProjectionY[4][d]->Draw("same");

		ProjectionY[5][d]->SetMinimum(0.01);
                ProjectionY[5][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[5][d]->SetLineColor(6);
               	ProjectionY[5][d]->Draw("same");

		ProjectionY[6][d]->SetMinimum(0.01);
                ProjectionY[6][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[6][d]->SetLineColor(7);
               	ProjectionY[6][d]->Draw("same");

		ProjectionY[7][d]->SetMinimum(0.01);
                ProjectionY[7][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[7][d]->SetLineColor(8);
               	ProjectionY[7][d]->Draw("same");

		ProjectionY[8][d]->SetMinimum(0.01);
                ProjectionY[8][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[8][d]->SetLineColor(9);
               	ProjectionY[8][d]->Draw("same");

		ProjectionY[9][d]->SetMinimum(0.01);
                ProjectionY[9][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY[9][d]->SetLineColor(44);
               	ProjectionY[9][d]->Draw("same");

 

		auto legend12 = new TLegend(0.7,0.7,0.88,0.88);
		legend12->AddEntry(ProjectionY[0][d], "pu = 0.5" ,"l");
		legend12->AddEntry(ProjectionY[1][d], "pu = 1" ,"l");
		legend12->AddEntry(ProjectionY[2][d], "pu = 1.5" ,"l");
		legend12->AddEntry(ProjectionY[3][d], "pu = 2" ,"l");
		legend12->AddEntry(ProjectionY[4][d], "pu = 10" ,"l");
		legend12->AddEntry(ProjectionY[5][d], "pu = 30" ,"l");
		legend12->AddEntry(ProjectionY[6][d], "pu = 50" ,"l");
		legend12->AddEntry(ProjectionY[7][d], "pu = 100" ,"l");
		legend12->AddEntry(ProjectionY[8][d], "pu = 140" ,"l");
		legend12->AddEntry(ProjectionY[9][d], "pu = 200" ,"l");
                legend12->SetFillColor(0);
  		legend12->SetLineColor(0);
  		legend12->SetFillColor(0);
		legend12->Draw();


		char *Projection12 = new char[10];
      		Int_t pro13 = d + 1;
      		sprintf(Projection12,"ProjectionbsRing1_%d.gif",pro13);
      		cout<<"=========================="<<Projection12<<endl;
      		c7->Print(outputpath + Projection12);
      		c7->Clear();





		TCanvas *c8 = new TCanvas();
            	c8->cd();
                c8->SetLogy();

                ProjectionY1[0][d]->SetMinimum(0.011);
                ProjectionY1[0][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[0][d]->SetLineColor(1);
               	ProjectionY1[0][d]->Draw("same");

		ProjectionY1[1][d]->SetMinimum(0.001);
                ProjectionY1[1][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[1][d]->SetLineColor(2);
               	ProjectionY1[1][d]->Draw("same");
               
                ProjectionY1[2][d]->SetMinimum(0.001);
                ProjectionY1[2][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[2][d]->SetLineColor(3);
               	ProjectionY1[2][d]->Draw("same");

                ProjectionY1[3][d]->SetMinimum(0.01);
                ProjectionY1[3][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[3][d]->SetLineColor(4);
               	ProjectionY1[3][d]->Draw("same");

                ProjectionY1[4][d]->SetMinimum(0.01);
                ProjectionY1[4][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[4][d]->SetLineColor(5);
               	ProjectionY1[4][d]->Draw("same");

		ProjectionY1[5][d]->SetMinimum(0.01);
                ProjectionY1[5][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[5][d]->SetLineColor(6);
               	ProjectionY1[5][d]->Draw("same");

		ProjectionY1[6][d]->SetMinimum(0.01);
                ProjectionY1[6][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[6][d]->SetLineColor(7);
               	ProjectionY1[6][d]->Draw("same");

		ProjectionY1[7][d]->SetMinimum(0.01);
                ProjectionY1[7][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[7][d]->SetLineColor(8);
               	ProjectionY1[7][d]->Draw("same");

		ProjectionY1[8][d]->SetMinimum(0.01);
                ProjectionY1[8][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[8][d]->SetLineColor(9);
               	ProjectionY1[8][d]->Draw("same");

		ProjectionY1[9][d]->SetMinimum(0.01);
                ProjectionY1[9][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY1[9][d]->SetLineColor(44);
               	ProjectionY1[9][d]->Draw("same");

 

		auto legend13 = new TLegend(0.7,0.7,0.88,0.88);
		legend13->AddEntry(ProjectionY[0][d], "pu = 0.5" ,"l");
		legend13->AddEntry(ProjectionY[1][d], "pu = 1" ,"l");
		legend13->AddEntry(ProjectionY[2][d], "pu = 1.5" ,"l");
		legend13->AddEntry(ProjectionY[3][d], "pu = 2" ,"l");
		legend13->AddEntry(ProjectionY[4][d], "pu = 10" ,"l");
		legend13->AddEntry(ProjectionY[5][d], "pu = 30" ,"l");
		legend13->AddEntry(ProjectionY[6][d], "pu = 50" ,"l");
		legend13->AddEntry(ProjectionY[7][d], "pu = 100" ,"l");
		legend13->AddEntry(ProjectionY[8][d], "pu = 140" ,"l");
		legend13->AddEntry(ProjectionY[9][d], "pu = 200" ,"l");
                legend13->SetFillColor(0);
  		legend13->SetLineColor(0);
  		legend13->SetFillColor(0);
		legend13->Draw();


		char *Projection14 = new char[10];
      		Int_t pro14 = d + 1;
      		sprintf(Projection14,"ProjectionbsRing2_%d.gif",pro14);
      		cout<<"=========================="<<Projection14<<endl;
      		c8->Print(outputpath + Projection14);
      		c8->Clear();



		TCanvas *c9 = new TCanvas();
            	c9->cd();
                c9->SetLogy();

                ProjectionY2[0][d]->SetMinimum(0.011);
                ProjectionY2[0][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[0][d]->SetLineColor(1);
               	ProjectionY2[0][d]->Draw("same");

		ProjectionY2[1][d]->SetMinimum(0.001);
                ProjectionY2[1][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[1][d]->SetLineColor(2);
               	ProjectionY2[1][d]->Draw("same");
               
                ProjectionY2[2][d]->SetMinimum(0.001);
                ProjectionY2[2][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[2][d]->SetLineColor(3);
               	ProjectionY2[2][d]->Draw("same");

                ProjectionY2[3][d]->SetMinimum(0.01);
                ProjectionY2[3][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[3][d]->SetLineColor(4);
               	ProjectionY2[3][d]->Draw("same");

                ProjectionY2[4][d]->SetMinimum(0.01);
                ProjectionY2[4][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[4][d]->SetLineColor(5);
               	ProjectionY2[4][d]->Draw("same");

		ProjectionY2[5][d]->SetMinimum(0.01);
                ProjectionY2[5][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[5][d]->SetLineColor(6);
               	ProjectionY2[5][d]->Draw("same");

		ProjectionY2[6][d]->SetMinimum(0.01);
                ProjectionY2[6][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[6][d]->SetLineColor(7);
               	ProjectionY2[6][d]->Draw("same");

		ProjectionY2[7][d]->SetMinimum(0.01);
                ProjectionY2[7][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[7][d]->SetLineColor(8);
               	ProjectionY2[7][d]->Draw("same");

		ProjectionY2[8][d]->SetMinimum(0.01);
                ProjectionY2[8][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[8][d]->SetLineColor(9);
               	ProjectionY2[8][d]->Draw("same");

		ProjectionY2[9][d]->SetMinimum(0.01);
                ProjectionY2[9][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY2[9][d]->SetLineColor(44);
               	ProjectionY2[9][d]->Draw("same");

 

		auto legend14 = new TLegend(0.7,0.7,0.88,0.88);
		legend14->AddEntry(ProjectionY[0][d], "pu = 0.5" ,"l");
		legend14->AddEntry(ProjectionY[1][d], "pu = 1" ,"l");
		legend14->AddEntry(ProjectionY[2][d], "pu = 1.5" ,"l");
		legend14->AddEntry(ProjectionY[3][d], "pu = 2" ,"l");
		legend14->AddEntry(ProjectionY[4][d], "pu = 10" ,"l");
		legend14->AddEntry(ProjectionY[5][d], "pu = 30" ,"l");
		legend14->AddEntry(ProjectionY[6][d], "pu = 50" ,"l");
		legend14->AddEntry(ProjectionY[7][d], "pu = 100" ,"l");
		legend14->AddEntry(ProjectionY[8][d], "pu = 140" ,"l");
		legend14->AddEntry(ProjectionY[9][d], "pu = 200" ,"l");
                legend14->SetFillColor(0);
  		legend14->SetLineColor(0);
  		legend14->SetFillColor(0);
		legend14->Draw();


		char *Projection15 = new char[10];
      		Int_t pro15 = d + 1;
      		sprintf(Projection15,"ProjectionbsRing3_%d.gif",pro15);
      		cout<<"=========================="<<Projection15<<endl;
      		c9->Print(outputpath + Projection15);
      		c9->Clear();



		TCanvas *c10 = new TCanvas();
            	c10->cd();
                c10->SetLogy();

                ProjectionY3[0][d]->SetMinimum(0.011);
                ProjectionY3[0][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[0][d]->SetLineColor(1);
               	ProjectionY3[0][d]->Draw("same");

		ProjectionY3[1][d]->SetMinimum(0.001);
                ProjectionY3[1][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[1][d]->SetLineColor(2);
               	ProjectionY3[1][d]->Draw("same");
               
                ProjectionY3[2][d]->SetMinimum(0.001);
                ProjectionY3[2][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[2][d]->SetLineColor(3);
               	ProjectionY3[2][d]->Draw("same");

                ProjectionY3[3][d]->SetMinimum(0.01);
                ProjectionY3[3][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[3][d]->SetLineColor(4);
               	ProjectionY3[3][d]->Draw("same");

                ProjectionY3[4][d]->SetMinimum(0.01);
                ProjectionY3[4][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[4][d]->SetLineColor(5);
               	ProjectionY3[4][d]->Draw("same");

		ProjectionY3[5][d]->SetMinimum(0.01);
                ProjectionY3[5][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[5][d]->SetLineColor(6);
               	ProjectionY3[5][d]->Draw("same");

		ProjectionY3[6][d]->SetMinimum(0.01);
                ProjectionY3[6][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[6][d]->SetLineColor(7);
               	ProjectionY3[6][d]->Draw("same");

		ProjectionY3[7][d]->SetMinimum(0.01);
                ProjectionY3[7][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[7][d]->SetLineColor(8);
               	ProjectionY3[7][d]->Draw("same");

		ProjectionY3[8][d]->SetMinimum(0.01);
                ProjectionY3[8][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[8][d]->SetLineColor(9);
               	ProjectionY3[8][d]->Draw("same");

		ProjectionY3[9][d]->SetMinimum(0.01);
                ProjectionY3[9][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY3[9][d]->SetLineColor(44);
               	ProjectionY3[9][d]->Draw("same");

 

		auto legend15 = new TLegend(0.7,0.7,0.88,0.88);
		legend15->AddEntry(ProjectionY3[0][d], "pu = 0.5" ,"l");
		legend15->AddEntry(ProjectionY3[1][d], "pu = 1" ,"l");
		legend15->AddEntry(ProjectionY3[2][d], "pu = 1.5" ,"l");
		legend15->AddEntry(ProjectionY3[3][d], "pu = 2" ,"l");
		legend15->AddEntry(ProjectionY3[4][d], "pu = 10" ,"l");
		legend15->AddEntry(ProjectionY3[5][d], "pu = 30" ,"l");
		legend15->AddEntry(ProjectionY3[6][d], "pu = 50" ,"l");
		legend15->AddEntry(ProjectionY3[7][d], "pu = 100" ,"l");
		legend15->AddEntry(ProjectionY3[8][d], "pu = 140" ,"l");
		legend15->AddEntry(ProjectionY3[9][d], "pu = 200" ,"l");
                legend15->SetFillColor(0);
  		legend15->SetLineColor(0);
  		legend15->SetFillColor(0);
		legend15->Draw();

		char *Projection16 = new char[10];
      		Int_t pro16 = d + 1;
      		sprintf(Projection16,"ProjectionbsRing4_%d.gif",pro16);
      		cout<<"=========================="<<Projection16<<endl;
      		c10->Print(outputpath + Projection16);
      		c10->Clear();



		TCanvas *c11 = new TCanvas();
            	c11->cd();
                c11->SetLogy();

                ProjectionY4[0][d]->SetMinimum(0.011);
                ProjectionY4[0][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[0][d]->SetLineColor(1);
               	ProjectionY4[0][d]->Draw("same");

		ProjectionY4[1][d]->SetMinimum(0.001);
                ProjectionY4[1][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[1][d]->SetLineColor(2);
               	ProjectionY4[1][d]->Draw("same");
               
                ProjectionY4[2][d]->SetMinimum(0.001);
                ProjectionY4[2][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[2][d]->SetLineColor(3);
               	ProjectionY4[2][d]->Draw("same");

                ProjectionY4[3][d]->SetMinimum(0.01);
                ProjectionY4[3][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[3][d]->SetLineColor(4);
               	ProjectionY4[3][d]->Draw("same");

                ProjectionY4[4][d]->SetMinimum(0.01);
                ProjectionY4[4][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[4][d]->SetLineColor(5);
               	ProjectionY4[4][d]->Draw("same");

		ProjectionY4[5][d]->SetMinimum(0.01);
                ProjectionY4[5][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[5][d]->SetLineColor(6);
               	ProjectionY4[5][d]->Draw("same");

		ProjectionY4[6][d]->SetMinimum(0.01);
                ProjectionY4[6][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[6][d]->SetLineColor(7);
               	ProjectionY4[6][d]->Draw("same");

		ProjectionY4[7][d]->SetMinimum(0.01);
                ProjectionY4[7][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[7][d]->SetLineColor(8);
               	ProjectionY4[7][d]->Draw("same");

		ProjectionY4[8][d]->SetMinimum(0.01);
                ProjectionY4[8][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[8][d]->SetLineColor(9);
               	ProjectionY4[8][d]->Draw("same");

		ProjectionY4[9][d]->SetMinimum(0.01);
                ProjectionY4[9][d]->GetXaxis()->SetRangeUser(0, 2500);
               	ProjectionY4[9][d]->SetLineColor(44);
               	ProjectionY4[9][d]->Draw("same");

 

		auto legend16 = new TLegend(0.7,0.7,0.88,0.88);
		legend16->AddEntry(ProjectionY4[0][d], "pu = 0.5" ,"l");
		legend16->AddEntry(ProjectionY4[1][d], "pu = 1" ,"l");
		legend16->AddEntry(ProjectionY4[2][d], "pu = 1.5" ,"l");
		legend16->AddEntry(ProjectionY4[3][d], "pu = 2" ,"l");
		legend16->AddEntry(ProjectionY4[4][d], "pu = 10" ,"l");
		legend16->AddEntry(ProjectionY4[5][d], "pu = 30" ,"l");
		legend16->AddEntry(ProjectionY4[6][d], "pu = 50" ,"l");
		legend16->AddEntry(ProjectionY4[7][d], "pu = 100" ,"l");
		legend16->AddEntry(ProjectionY4[8][d], "pu = 140" ,"l");
		legend16->AddEntry(ProjectionY4[9][d], "pu = 200" ,"l");
                legend16->SetFillColor(0);
  		legend16->SetLineColor(0);
  		legend16->SetFillColor(0);
		legend16->Draw();


		char *Projection17 = new char[10];
      		Int_t pro17 = d + 1;
      		sprintf(Projection17,"ProjectionbsRing5_%d.gif",pro17);
      		cout<<"=========================="<<Projection17<<endl;
      		c11->Print(outputpath + Projection17);
      		c11->Clear();
               
	

	    
}
	    
            
            if (option == 3){
            
                 //draw the 2D 2x Coincidences distributions for all PU
            TString histoname = "BRIL_IT_Analysis/TEPX/3xCoincidences/Number of 3x Coincidences for Disk ";

            TH2F* H = (TH2F*)F.Get(histoname + disklist[d].c_str());
	    Histogram2D[pu][d] = (TH2F*)H->Clone(TString(H->GetName()));

            ProjectionY[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 1, 1)->Clone(TString(H->GetName()) +             "Projection");
            ProjectionY1[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 2,2)->Clone(TString(H->GetName()) +             "Projection");
            ProjectionY2[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 3, 3)->Clone(TString(H->GetName()) +            "Projection");
	    ProjectionY3[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 4, 4)->Clone(TString(H->GetName()) +            "Projection");
	    ProjectionY4[pu][d] = (TH1D*)H->ProjectionY(TString("ProjectionY_") + pu + "_" + d, 5, 5)->Clone(TString(H->GetName()) +            "Projection");
            
            Prof_TEPXClustersPerEvent[pu][d] = (TProfile*)H->ProfileX()->Clone(TString(H->GetName()) + "Profile"); //number of cluster vs ring

		if (pu < 4){
	        TCanvas *c1 = new TCanvas();
            	c1->cd();
                Histogram2D[pu][d]->GetYaxis()->SetRangeUser(0, 50);
      		Histogram2D[pu][d]->Draw("colz");
      		char *hist2Dname = new char[10];
      		Int_t pud = pu * 8 + d + 1;
      		sprintf(hist2Dname,"histoss2D_%d.gif",pud);
      		cout<<"=========================="<<hist2Dname<<endl;
      		c1->Print(outputpath + hist2Dname);
      		c1->Clear();
                
                TCanvas *c2 = new TCanvas();
            	c2->cd();
                ProjectionY[pu][d]->GetXaxis()->SetRangeUser(0, 60);
                ProjectionY[pu][d]->GetYaxis()->SetRangeUser(0, 70000);
                ProjectionY[pu][d]->SetLineColor(5);
      		ProjectionY[pu][d]->Draw("same");
                ProjectionY1[pu][d]->GetXaxis()->SetRangeUser(0, 60);
                ProjectionY1[pu][d]->GetYaxis()->SetRangeUser(0, 70000);
		ProjectionY1[pu][d]->SetLineColor(2);
      		ProjectionY1[pu][d]->Draw("same");
                ProjectionY2[pu][d]->GetXaxis()->SetRangeUser(0, 60);
                ProjectionY2[pu][d]->GetYaxis()->SetRangeUser(0, 70000);
		ProjectionY2[pu][d]->SetLineColor(3);
      		ProjectionY2[pu][d]->Draw("same");
		ProjectionY3[pu][d]->GetXaxis()->SetRangeUser(0, 60);
                ProjectionY3[pu][d]->GetYaxis()->SetRangeUser(0, 70000);
		ProjectionY3[pu][d]->SetLineColor(4);
      		ProjectionY3[pu][d]->Draw("same");
		ProjectionY4[pu][d]->GetXaxis()->SetRangeUser(0, 60);
                ProjectionY4[pu][d]->GetYaxis()->SetRangeUser(0, 70000);
		ProjectionY4[pu][d]->SetLineColor(6);
      		ProjectionY4[pu][d]->Draw("same");
                
               

                auto legend10 = new TLegend(0.7,0.7,0.88,0.88);
		legend10->AddEntry(ProjectionY[pu][d], "Ring 1" ,"l");
		legend10->AddEntry(ProjectionY1[pu][d], "Ring 2" ,"l");
		legend10->AddEntry(ProjectionY2[pu][d], "Ring 3" ,"l");
		legend10->AddEntry(ProjectionY3[pu][d], "Ring 4" ,"l");
                legend10->AddEntry(ProjectionY4[pu][d], "Ring 5" ,"l");

                legend10->SetFillColor(0);
  		legend10->SetLineColor(0);
  		legend10->SetFillColor(0);
		legend10->Draw();
      		char *Projection = new char[12];
      		Int_t pro = pu * 8 + d + 1;
      		sprintf(Projection,"Projection3x_%d.gif",pro);
      		cout<<"=========================="<<Projection<<endl;
      		c2->Print(outputpath + Projection);
      		c2->Clear();

	}
                if (pu >= 4){

	        TCanvas *c3 = new TCanvas();
            	c3->cd();
                Histogram2D[pu][d]->GetYaxis()->SetRangeUser(0, 100);
      		Histogram2D[pu][d]->Draw("colz");
      		char *hist2Dname = new char[10];
      		Int_t pud1 = pu * 8 + d + 1;
      		sprintf(hist2Dname,"histoss2D_%d.gif",pud1);
      		cout<<"=========================="<<hist2Dname<<endl;
      		c3->Print(outputpath + hist2Dname);
      		c3->Clear();

                
                TCanvas *c4= new TCanvas();
            	c4->cd();
                ProjectionY[pu][d]->GetXaxis()->SetRangeUser(0, 1500);
                ProjectionY[pu][d]->GetYaxis()->SetRangeUser(0, 500);
                ProjectionY[pu][d]->SetLineColor(5);
      		ProjectionY[pu][d]->Draw("same");
                ProjectionY1[pu][d]->GetXaxis()->SetRangeUser(0, 1500);
                ProjectionY1[pu][d]->GetYaxis()->SetRangeUser(0, 500);
		ProjectionY1[pu][d]->SetLineColor(2);
      		ProjectionY1[pu][d]->Draw("same");
                ProjectionY2[pu][d]->GetXaxis()->SetRangeUser(0, 1500);
                ProjectionY2[pu][d]->GetYaxis()->SetRangeUser(0, 500);
		ProjectionY2[pu][d]->SetLineColor(3);
      		ProjectionY2[pu][d]->Draw("same");
		ProjectionY3[pu][d]->GetXaxis()->SetRangeUser(0, 1500);
                ProjectionY3[pu][d]->GetYaxis()->SetRangeUser(0, 500);
		ProjectionY3[pu][d]->SetLineColor(4);
      		ProjectionY3[pu][d]->Draw("same");
		ProjectionY4[pu][d]->GetXaxis()->SetRangeUser(0, 1500);
                ProjectionY4[pu][d]->GetYaxis()->SetRangeUser(0, 500);
		ProjectionY4[pu][d]->SetLineColor(6);
      		ProjectionY4[pu][d]->Draw("same");

                auto legend10 = new TLegend(0.7,0.7,0.88,0.88);
		legend10->AddEntry(ProjectionY[pu][d], "Ring 1" ,"l");
		legend10->AddEntry(ProjectionY1[pu][d], "Ring 2" ,"l");
		legend10->AddEntry(ProjectionY2[pu][d], "Ring 3" ,"l");
		legend10->AddEntry(ProjectionY3[pu][d], "Ring 4" ,"l");
                legend10->AddEntry(ProjectionY4[pu][d], "Ring 5" ,"l");

                legend10->SetFillColor(0);
  		legend10->SetLineColor(0);
  		legend10->SetFillColor(0);
		legend10->Draw();
      		char *Projection = new char[12];
      		Int_t pro = pu * 8 + d + 1;
      		sprintf(Projection,"Projection3x_%d.gif",pro);
      		cout<<"=========================="<<Projection<<endl;
      		c4->Print(outputpath + Projection);
      		c4->Clear();
	      
	}

                
                if(pu < 4){

                TCanvas *c6 = new TCanvas();
            	c6->cd();
                c6->SetLogy();

                ProjectionY[0][d]->SetMinimum(0.011);
                ProjectionY[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[0][d]->SetLineColor(1);
               	ProjectionY[0][d]->Draw("same");

		ProjectionY[1][d]->SetMinimum(0.001);
                ProjectionY[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[1][d]->SetLineColor(2);
               	ProjectionY[1][d]->Draw("same");
               
                ProjectionY[2][d]->SetMinimum(0.001);
                ProjectionY[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[2][d]->SetLineColor(3);
               	ProjectionY[2][d]->Draw("same");

                ProjectionY[3][d]->SetMinimum(0.01);
                ProjectionY[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[3][d]->SetLineColor(4);
               	ProjectionY[3][d]->Draw("same");

		

		auto legend1 = new TLegend(0.7,0.7,0.88,0.88);
		legend1->AddEntry(ProjectionY[0][d], "pu = 0.5" ,"l");
		legend1->AddEntry(ProjectionY[1][d], "pu = 1" ,"l");
		legend1->AddEntry(ProjectionY[2][d], "pu = 1.5" ,"l");
		legend1->AddEntry(ProjectionY[3][d], "pu = 2" ,"l");

                legend1->SetFillColor(0);
  		legend1->SetLineColor(0);
  		legend1->SetFillColor(0);
		legend1->Draw();

		char *Projection11 = new char[10];
      		Int_t pro12 = d + 1;
      		sprintf(Projection11,"ProjectioncsRing1_%d.gif",pro12);
      		cout<<"=========================="<<Projection11<<endl;
      		c6->Print(outputpath + Projection11);
      		c6->Clear();

                
		TCanvas *c7 = new TCanvas();
            	c7->cd();
                c7->SetLogy();

		ProjectionY1[0][d]->SetMinimum(0.011);
                ProjectionY1[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[0][d]->SetLineColor(1);
               	ProjectionY1[0][d]->Draw("same");

		ProjectionY1[1][d]->SetMinimum(0.001);
                ProjectionY1[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[1][d]->SetLineColor(2);
               	ProjectionY1[1][d]->Draw("same");
               
                ProjectionY1[2][d]->SetMinimum(0.001);
                ProjectionY1[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[2][d]->SetLineColor(3);
               	ProjectionY1[2][d]->Draw("same");

                ProjectionY1[3][d]->SetMinimum(0.01);
                ProjectionY1[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[3][d]->SetLineColor(4);
               	ProjectionY1[3][d]->Draw("same");

		auto legend2 = new TLegend(0.7,0.7,0.88,0.88);
		legend2->AddEntry(ProjectionY1[0][d], "pu = 0.5" ,"l");
		legend2->AddEntry(ProjectionY1[1][d], "pu = 1" ,"l");
		legend2->AddEntry(ProjectionY1[2][d], "pu = 1.5" ,"l");
		legend2->AddEntry(ProjectionY1[3][d], "pu = 2" ,"l");

                legend2->SetFillColor(0);
  		legend2->SetLineColor(0);
  		legend2->SetFillColor(0);
		legend2->Draw();

		char *Projection14 = new char[10];
      		Int_t pro14 = d + 1;
      		sprintf(Projection14,"ProjectioncsRing2_%d.gif",pro14);
      		cout<<"=========================="<<Projection14<<endl;
      		c7->Print(outputpath + Projection14);
      		c7->Clear();

		TCanvas *c8 = new TCanvas();
            	c8->cd();
                c8->SetLogy();

                ProjectionY2[0][d]->SetMinimum(0.011);
                ProjectionY2[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[0][d]->SetLineColor(1);
               	ProjectionY2[0][d]->Draw("same");

		ProjectionY2[1][d]->SetMinimum(0.001);
                ProjectionY2[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[1][d]->SetLineColor(2);
               	ProjectionY2[1][d]->Draw("same");
               
                ProjectionY2[2][d]->SetMinimum(0.001);
                ProjectionY2[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[2][d]->SetLineColor(3);
               	ProjectionY2[2][d]->Draw("same");

                ProjectionY2[3][d]->SetMinimum(0.01);
                ProjectionY2[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[3][d]->SetLineColor(4);
               	ProjectionY2[3][d]->Draw("same");

		auto legend3 = new TLegend(0.7,0.7,0.88,0.88);
		legend3->AddEntry(ProjectionY1[0][d], "pu = 0.5" ,"l");
		legend3->AddEntry(ProjectionY1[1][d], "pu = 1" ,"l");
		legend3->AddEntry(ProjectionY1[2][d], "pu = 1.5" ,"l");
		legend3->AddEntry(ProjectionY1[3][d], "pu = 2" ,"l");

                legend3->SetFillColor(0);
  		legend3->SetLineColor(0);
  		legend3->SetFillColor(0);
		legend3->Draw();

		char *Projection15 = new char[10];
      		Int_t pro15 = d + 1;
      		sprintf(Projection15,"ProjectioncsRing3_%d.gif",pro15);
      		cout<<"=========================="<<Projection15<<endl;
      		c8->Print(outputpath + Projection15);
      		c8->Clear();


		TCanvas *c9 = new TCanvas();
            	c9->cd();
                c9->SetLogy();

                ProjectionY3[0][d]->SetMinimum(0.011);
                ProjectionY3[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[0][d]->SetLineColor(1);
               	ProjectionY3[0][d]->Draw("same");

		ProjectionY3[1][d]->SetMinimum(0.001);
                ProjectionY3[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[1][d]->SetLineColor(2);
               	ProjectionY3[1][d]->Draw("same");
               
                ProjectionY3[2][d]->SetMinimum(0.001);
                ProjectionY3[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[2][d]->SetLineColor(3);
               	ProjectionY3[2][d]->Draw("same");

                ProjectionY3[3][d]->SetMinimum(0.01);
                ProjectionY3[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[3][d]->SetLineColor(4);
               	ProjectionY3[3][d]->Draw("same");

		auto legend4 = new TLegend(0.7,0.7,0.88,0.88);
		legend4->AddEntry(ProjectionY3[0][d], "pu = 0.5" ,"l");
		legend4->AddEntry(ProjectionY3[1][d], "pu = 1" ,"l");
		legend4->AddEntry(ProjectionY3[2][d], "pu = 1.5" ,"l");
		legend4->AddEntry(ProjectionY3[3][d], "pu = 2" ,"l");

                legend4->SetFillColor(0);
  		legend4->SetLineColor(0);
  		legend4->SetFillColor(0);
		legend4->Draw();

		char *Projection16 = new char[10];
      		Int_t pro16 = d + 1;
      		sprintf(Projection16,"ProjectioncsRing4_%d.gif",pro16);
      		cout<<"=========================="<<Projection16<<endl;
      		c8->Print(outputpath + Projection16);
      		c8->Clear();



		TCanvas *c10 = new TCanvas();
            	c10->cd();
                c10->SetLogy();

                ProjectionY4[0][d]->SetMinimum(0.011);
                ProjectionY4[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[0][d]->SetLineColor(1);
               	ProjectionY4[0][d]->Draw("same");

		ProjectionY4[1][d]->SetMinimum(0.001);
                ProjectionY4[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[1][d]->SetLineColor(2);
               	ProjectionY4[1][d]->Draw("same");
               
                ProjectionY4[2][d]->SetMinimum(0.001);
                ProjectionY4[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[2][d]->SetLineColor(3);
               	ProjectionY4[2][d]->Draw("same");

                ProjectionY4[3][d]->SetMinimum(0.01);
                ProjectionY4[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[3][d]->SetLineColor(4);
               	ProjectionY4[3][d]->Draw("same");

		auto legend5 = new TLegend(0.7,0.7,0.88,0.88);
		legend5->AddEntry(ProjectionY4[0][d], "pu = 0.5" ,"l");
		legend5->AddEntry(ProjectionY4[1][d], "pu = 1" ,"l");
		legend5->AddEntry(ProjectionY4[2][d], "pu = 1.5" ,"l");
		legend5->AddEntry(ProjectionY4[3][d], "pu = 2" ,"l");

                legend5->SetFillColor(0);
  		legend5->SetLineColor(0);
  		legend5->SetFillColor(0);
		legend5->Draw();

		char *Projection17 = new char[10];
      		Int_t pro17 = d + 1;
      		sprintf(Projection17,"ProjectioncsRing5_%d.gif",pro17);
      		cout<<"=========================="<<Projection17<<endl;
      		c10->Print(outputpath + Projection17);
      		c10->Clear();

		}

                
		TCanvas *c7 = new TCanvas();
            	c7->cd();
                c7->SetLogy();

                ProjectionY[0][d]->SetMinimum(0.011);
                ProjectionY[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[0][d]->SetLineColor(1);
               	ProjectionY[0][d]->Draw("same");

		ProjectionY[1][d]->SetMinimum(0.001);
                ProjectionY[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[1][d]->SetLineColor(2);
               	ProjectionY[1][d]->Draw("same");
               
                ProjectionY[2][d]->SetMinimum(0.001);
                ProjectionY[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[2][d]->SetLineColor(3);
               	ProjectionY[2][d]->Draw("same");

                ProjectionY[3][d]->SetMinimum(0.01);
                ProjectionY[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[3][d]->SetLineColor(4);
               	ProjectionY[3][d]->Draw("same");

                ProjectionY[4][d]->SetMinimum(0.01);
                ProjectionY[4][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[4][d]->SetLineColor(5);
               	ProjectionY[4][d]->Draw("same");

		ProjectionY[5][d]->SetMinimum(0.01);
                ProjectionY[5][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[5][d]->SetLineColor(6);
               	ProjectionY[5][d]->Draw("same");

		ProjectionY[6][d]->SetMinimum(0.01);
                ProjectionY[6][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[6][d]->SetLineColor(7);
               	ProjectionY[6][d]->Draw("same");

		ProjectionY[7][d]->SetMinimum(0.01);
                ProjectionY[7][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[7][d]->SetLineColor(8);
               	ProjectionY[7][d]->Draw("same");

		ProjectionY[8][d]->SetMinimum(0.01);
                ProjectionY[8][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[8][d]->SetLineColor(9);
               	ProjectionY[8][d]->Draw("same");

		ProjectionY[9][d]->SetMinimum(0.01);
                ProjectionY[9][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY[9][d]->SetLineColor(44);
               	ProjectionY[9][d]->Draw("same");

 

		auto legend12 = new TLegend(0.7,0.7,0.88,0.88);
		legend12->AddEntry(ProjectionY[0][d], "pu = 0.5" ,"l");
		legend12->AddEntry(ProjectionY[1][d], "pu = 1" ,"l");
		legend12->AddEntry(ProjectionY[2][d], "pu = 1.5" ,"l");
		legend12->AddEntry(ProjectionY[3][d], "pu = 2" ,"l");
		legend12->AddEntry(ProjectionY[4][d], "pu = 10" ,"l");
		legend12->AddEntry(ProjectionY[5][d], "pu = 30" ,"l");
		legend12->AddEntry(ProjectionY[6][d], "pu = 50" ,"l");
		legend12->AddEntry(ProjectionY[7][d], "pu = 100" ,"l");
		legend12->AddEntry(ProjectionY[8][d], "pu = 140" ,"l");
		legend12->AddEntry(ProjectionY[9][d], "pu = 200" ,"l");
                legend12->SetFillColor(0);
  		legend12->SetLineColor(0);
  		legend12->SetFillColor(0);
		legend12->Draw();


		char *Projection12 = new char[10];
      		Int_t pro13 = d + 1;
      		sprintf(Projection12,"ProjectiondsRing1_%d.gif",pro13);
      		cout<<"=========================="<<Projection12<<endl;
      		c7->Print(outputpath + Projection12);
      		c7->Clear();





		TCanvas *c8 = new TCanvas();
            	c8->cd();
                c8->SetLogy();

                ProjectionY1[0][d]->SetMinimum(0.011);
                ProjectionY1[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[0][d]->SetLineColor(1);
               	ProjectionY1[0][d]->Draw("same");

		ProjectionY1[1][d]->SetMinimum(0.001);
                ProjectionY1[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[1][d]->SetLineColor(2);
               	ProjectionY1[1][d]->Draw("same");
               
                ProjectionY1[2][d]->SetMinimum(0.001);
                ProjectionY1[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[2][d]->SetLineColor(3);
               	ProjectionY1[2][d]->Draw("same");

                ProjectionY1[3][d]->SetMinimum(0.01);
                ProjectionY1[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[3][d]->SetLineColor(4);
               	ProjectionY1[3][d]->Draw("same");

                ProjectionY1[4][d]->SetMinimum(0.01);
                ProjectionY1[4][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[4][d]->SetLineColor(5);
               	ProjectionY1[4][d]->Draw("same");

		ProjectionY1[5][d]->SetMinimum(0.01);
                ProjectionY1[5][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[5][d]->SetLineColor(6);
               	ProjectionY1[5][d]->Draw("same");

		ProjectionY1[6][d]->SetMinimum(0.01);
                ProjectionY1[6][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[6][d]->SetLineColor(7);
               	ProjectionY1[6][d]->Draw("same");

		ProjectionY1[7][d]->SetMinimum(0.01);
                ProjectionY1[7][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[7][d]->SetLineColor(8);
               	ProjectionY1[7][d]->Draw("same");

		ProjectionY1[8][d]->SetMinimum(0.01);
                ProjectionY1[8][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[8][d]->SetLineColor(9);
               	ProjectionY1[8][d]->Draw("same");

		ProjectionY1[9][d]->SetMinimum(0.01);
                ProjectionY1[9][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY1[9][d]->SetLineColor(44);
               	ProjectionY1[9][d]->Draw("same");

 

		auto legend13 = new TLegend(0.7,0.7,0.88,0.88);
		legend13->AddEntry(ProjectionY[0][d], "pu = 0.5" ,"l");
		legend13->AddEntry(ProjectionY[1][d], "pu = 1" ,"l");
		legend13->AddEntry(ProjectionY[2][d], "pu = 1.5" ,"l");
		legend13->AddEntry(ProjectionY[3][d], "pu = 2" ,"l");
		legend13->AddEntry(ProjectionY[4][d], "pu = 10" ,"l");
		legend13->AddEntry(ProjectionY[5][d], "pu = 30" ,"l");
		legend13->AddEntry(ProjectionY[6][d], "pu = 50" ,"l");
		legend13->AddEntry(ProjectionY[7][d], "pu = 100" ,"l");
		legend13->AddEntry(ProjectionY[8][d], "pu = 140" ,"l");
		legend13->AddEntry(ProjectionY[9][d], "pu = 200" ,"l");
                legend13->SetFillColor(0);
  		legend13->SetLineColor(0);
  		legend13->SetFillColor(0);
		legend13->Draw();


		char *Projection14 = new char[10];
      		Int_t pro14 = d + 1;
      		sprintf(Projection14,"ProjectiondsRing2_%d.gif",pro14);
      		cout<<"=========================="<<Projection14<<endl;
      		c8->Print(outputpath + Projection14);
      		c8->Clear();



		TCanvas *c9 = new TCanvas();
            	c9->cd();
                c9->SetLogy();

                ProjectionY2[0][d]->SetMinimum(0.011);
                ProjectionY2[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[0][d]->SetLineColor(1);
               	ProjectionY2[0][d]->Draw("same");

		ProjectionY2[1][d]->SetMinimum(0.001);
                ProjectionY2[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[1][d]->SetLineColor(2);
               	ProjectionY2[1][d]->Draw("same");
               
                ProjectionY2[2][d]->SetMinimum(0.001);
                ProjectionY2[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[2][d]->SetLineColor(3);
               	ProjectionY2[2][d]->Draw("same");

                ProjectionY2[3][d]->SetMinimum(0.01);
                ProjectionY2[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[3][d]->SetLineColor(4);
               	ProjectionY2[3][d]->Draw("same");

                ProjectionY2[4][d]->SetMinimum(0.01);
                ProjectionY2[4][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[4][d]->SetLineColor(5);
               	ProjectionY2[4][d]->Draw("same");

		ProjectionY2[5][d]->SetMinimum(0.01);
                ProjectionY2[5][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[5][d]->SetLineColor(6);
               	ProjectionY2[5][d]->Draw("same");

		ProjectionY2[6][d]->SetMinimum(0.01);
                ProjectionY2[6][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[6][d]->SetLineColor(7);
               	ProjectionY2[6][d]->Draw("same");

		ProjectionY2[7][d]->SetMinimum(0.01);
                ProjectionY2[7][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[7][d]->SetLineColor(8);
               	ProjectionY2[7][d]->Draw("same");

		ProjectionY2[8][d]->SetMinimum(0.01);
                ProjectionY2[8][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[8][d]->SetLineColor(9);
               	ProjectionY2[8][d]->Draw("same");

		ProjectionY2[9][d]->SetMinimum(0.01);
                ProjectionY2[9][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY2[9][d]->SetLineColor(44);
               	ProjectionY2[9][d]->Draw("same");

 

		auto legend14 = new TLegend(0.7,0.7,0.88,0.88);
		legend14->AddEntry(ProjectionY[0][d], "pu = 0.5" ,"l");
		legend14->AddEntry(ProjectionY[1][d], "pu = 1" ,"l");
		legend14->AddEntry(ProjectionY[2][d], "pu = 1.5" ,"l");
		legend14->AddEntry(ProjectionY[3][d], "pu = 2" ,"l");
		legend14->AddEntry(ProjectionY[4][d], "pu = 10" ,"l");
		legend14->AddEntry(ProjectionY[5][d], "pu = 30" ,"l");
		legend14->AddEntry(ProjectionY[6][d], "pu = 50" ,"l");
		legend14->AddEntry(ProjectionY[7][d], "pu = 100" ,"l");
		legend14->AddEntry(ProjectionY[8][d], "pu = 140" ,"l");
		legend14->AddEntry(ProjectionY[9][d], "pu = 200" ,"l");
                legend14->SetFillColor(0);
  		legend14->SetLineColor(0);
  		legend14->SetFillColor(0);
		legend14->Draw();


		char *Projection15 = new char[10];
      		Int_t pro15 = d + 1;
      		sprintf(Projection15,"ProjectiondsRing3_%d.gif",pro15);
      		cout<<"=========================="<<Projection15<<endl;
      		c9->Print(outputpath + Projection15);
      		c9->Clear();



		TCanvas *c10 = new TCanvas();
            	c10->cd();
                c10->SetLogy();

                ProjectionY3[0][d]->SetMinimum(0.011);
                ProjectionY3[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[0][d]->SetLineColor(1);
               	ProjectionY3[0][d]->Draw("same");

		ProjectionY3[1][d]->SetMinimum(0.001);
                ProjectionY3[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[1][d]->SetLineColor(2);
               	ProjectionY3[1][d]->Draw("same");
               
                ProjectionY3[2][d]->SetMinimum(0.001);
                ProjectionY3[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[2][d]->SetLineColor(3);
               	ProjectionY3[2][d]->Draw("same");

                ProjectionY3[3][d]->SetMinimum(0.01);
                ProjectionY3[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[3][d]->SetLineColor(4);
               	ProjectionY3[3][d]->Draw("same");

                ProjectionY3[4][d]->SetMinimum(0.01);
                ProjectionY3[4][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[4][d]->SetLineColor(5);
               	ProjectionY3[4][d]->Draw("same");

		ProjectionY3[5][d]->SetMinimum(0.01);
                ProjectionY3[5][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[5][d]->SetLineColor(6);
               	ProjectionY3[5][d]->Draw("same");

		ProjectionY3[6][d]->SetMinimum(0.01);
                ProjectionY3[6][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[6][d]->SetLineColor(7);
               	ProjectionY3[6][d]->Draw("same");

		ProjectionY3[7][d]->SetMinimum(0.01);
                ProjectionY3[7][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[7][d]->SetLineColor(8);
               	ProjectionY3[7][d]->Draw("same");

		ProjectionY3[8][d]->SetMinimum(0.01);
                ProjectionY3[8][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[8][d]->SetLineColor(9);
               	ProjectionY3[8][d]->Draw("same");

		ProjectionY3[9][d]->SetMinimum(0.01);
                ProjectionY3[9][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY3[9][d]->SetLineColor(44);
               	ProjectionY3[9][d]->Draw("same");

 

		auto legend15 = new TLegend(0.7,0.7,0.88,0.88);
		legend15->AddEntry(ProjectionY3[0][d], "pu = 0.5" ,"l");
		legend15->AddEntry(ProjectionY3[1][d], "pu = 1" ,"l");
		legend15->AddEntry(ProjectionY3[2][d], "pu = 1.5" ,"l");
		legend15->AddEntry(ProjectionY3[3][d], "pu = 2" ,"l");
		legend15->AddEntry(ProjectionY3[4][d], "pu = 10" ,"l");
		legend15->AddEntry(ProjectionY3[5][d], "pu = 30" ,"l");
		legend15->AddEntry(ProjectionY3[6][d], "pu = 50" ,"l");
		legend15->AddEntry(ProjectionY3[7][d], "pu = 100" ,"l");
		legend15->AddEntry(ProjectionY3[8][d], "pu = 140" ,"l");
		legend15->AddEntry(ProjectionY3[9][d], "pu = 200" ,"l");
                legend15->SetFillColor(0);
  		legend15->SetLineColor(0);
  		legend15->SetFillColor(0);
		legend15->Draw();

		char *Projection16 = new char[10];
      		Int_t pro16 = d + 1;
      		sprintf(Projection16,"ProjectiondsRing4_%d.gif",pro16);
      		cout<<"=========================="<<Projection16<<endl;
      		c10->Print(outputpath + Projection16);
      		c10->Clear();



		TCanvas *c11 = new TCanvas();
            	c11->cd();
                c11->SetLogy();

                ProjectionY4[0][d]->SetMinimum(0.011);
                ProjectionY4[0][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[0][d]->SetLineColor(1);
               	ProjectionY4[0][d]->Draw("same");

		ProjectionY4[1][d]->SetMinimum(0.001);
                ProjectionY4[1][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[1][d]->SetLineColor(2);
               	ProjectionY4[1][d]->Draw("same");
               
                ProjectionY4[2][d]->SetMinimum(0.001);
                ProjectionY4[2][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[2][d]->SetLineColor(3);
               	ProjectionY4[2][d]->Draw("same");

                ProjectionY4[3][d]->SetMinimum(0.01);
                ProjectionY4[3][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[3][d]->SetLineColor(4);
               	ProjectionY4[3][d]->Draw("same");

                ProjectionY4[4][d]->SetMinimum(0.01);
                ProjectionY4[4][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[4][d]->SetLineColor(5);
               	ProjectionY4[4][d]->Draw("same");

		ProjectionY4[5][d]->SetMinimum(0.01);
                ProjectionY4[5][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[5][d]->SetLineColor(6);
               	ProjectionY4[5][d]->Draw("same");

		ProjectionY4[6][d]->SetMinimum(0.01);
                ProjectionY4[6][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[6][d]->SetLineColor(7);
               	ProjectionY4[6][d]->Draw("same");

		ProjectionY4[7][d]->SetMinimum(0.01);
                ProjectionY4[7][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[7][d]->SetLineColor(8);
               	ProjectionY4[7][d]->Draw("same");

		ProjectionY4[8][d]->SetMinimum(0.01);
                ProjectionY4[8][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[8][d]->SetLineColor(9);
               	ProjectionY4[8][d]->Draw("same");

		ProjectionY4[9][d]->SetMinimum(0.01);
                ProjectionY4[9][d]->GetXaxis()->SetRangeUser(0, 100);
               	ProjectionY4[9][d]->SetLineColor(44);
               	ProjectionY4[9][d]->Draw("same");

 

		auto legend16 = new TLegend(0.7,0.7,0.88,0.88);
		legend16->AddEntry(ProjectionY4[0][d], "pu = 0.5" ,"l");
		legend16->AddEntry(ProjectionY4[1][d], "pu = 1" ,"l");
		legend16->AddEntry(ProjectionY4[2][d], "pu = 1.5" ,"l");
		legend16->AddEntry(ProjectionY4[3][d], "pu = 2" ,"l");
		legend16->AddEntry(ProjectionY4[4][d], "pu = 10" ,"l");
		legend16->AddEntry(ProjectionY4[5][d], "pu = 30" ,"l");
		legend16->AddEntry(ProjectionY4[6][d], "pu = 50" ,"l");
		legend16->AddEntry(ProjectionY4[7][d], "pu = 100" ,"l");
		legend16->AddEntry(ProjectionY4[8][d], "pu = 140" ,"l");
		legend16->AddEntry(ProjectionY4[9][d], "pu = 200" ,"l");
                legend16->SetFillColor(0);
  		legend16->SetLineColor(0);
  		legend16->SetFillColor(0);
		legend16->Draw();


		char *Projection17 = new char[10];
      		Int_t pro17 = d + 1;
      		sprintf(Projection17,"ProjectiondsRing5_%d.gif",pro17);
      		cout<<"=========================="<<Projection17<<endl;
      		c11->Print(outputpath + Projection17);
      		c11->Clear();
	

	    }

	 
     
	    
            //Bin content and Bin error of the  Profile X of Histogram H
            for (int r = 0; r < 5; r++) {
	      
                TEPXClustersPerEvent[d][r]->SetPoint(pu, pumap[pulist[pu]], Prof_TEPXClustersPerEvent[pu][d]->GetBinContent(r + 1));
                TEPXClustersPerEvent[d][r]->SetPointError(pu, 0, Prof_TEPXClustersPerEvent[pu][d]->GetBinError(r + 1));
              
            }
        }
    }


    
    

    //Make the fit graphs for all Rings and Disks
    TLatex label;
    label.SetTextSize(0.2);
    TF1* FitTEPXClustersPerEvent[20][4];

    for (long d = 0; d < disklist.size(); d++)
        for (long r = 0; r < 5; r++) {

            //fix the pads
            Int_t l = d * 5 + r + 1;
            TCanvas C("C");
            C.cd();

            //fit
            FitTEPXClustersPerEvent[d][r] = new TF1(TString("Fit_") + d + "_" + r, "[0]+[1]*x", 0.5, 2);
            FitTEPXClustersPerEvent[d][r]->SetLineColor(4);
            TEPXClustersPerEvent[d][r]->Fit(FitTEPXClustersPerEvent[d][r], "", "", 0.5, 2);
            
            
             
            //draw the fit graphs for all disks and rings
            TEPXClustersPerEvent[d][r]->GetYaxis()->SetNdivisions(12);
            TEPXClustersPerEvent[d][r]->GetYaxis()->SetLabelSize(0.03);
            TEPXClustersPerEvent[d][r]->GetXaxis()->SetLabelSize(0.03);
            TEPXClustersPerEvent[d][r]->GetXaxis()->SetTitle("Pileup");
            TEPXClustersPerEvent[d][r]->SetMarkerStyle(21);
            TEPXClustersPerEvent[d][r]->SetMarkerSize(1);
            TEPXClustersPerEvent[d][r]->SetMarkerColor(2);
            gStyle->SetOptStat(1111);
	    gStyle->SetOptFit(1111);
            gPad->SetGrid(1, 1);

            if (option == 1){

            TEPXClustersPerEvent[d][r]->GetYaxis()->SetTitle("Mean Number of Clusters");
            TEPXClustersPerEvent[d][r]->GetYaxis()->SetRangeUser(0, 30);
            TEPXClustersPerEvent[d][r]->GetXaxis()->SetNdivisions(12);
            TEPXClustersPerEvent[d][r]->GetXaxis()->SetRangeUser(0, 2);
            TEPXClustersPerEvent[d][r]->Draw("ape");
            FitTEPXClustersPerEvent[d][r]->Draw("lsame");
            label.SetTextSize(0.1);
            label.DrawLatexNDC(0.15, 0.8, TString("D=") + (long)(d <= 3 ? d - 4 : d + 1 - 4) + ",  R=" + (long)(r+1));
            char* histname = new char[10];
            sprintf(histname, "histo%d_linearity1.gif", l);
            cout << "==========================" << histname << endl;
            C.Print(outputpath + histname);
            C.Clear();
            
        }
	

            if (option == 2){
            
            TEPXClustersPerEvent[d][r]->GetYaxis()->SetTitle("Mean Number of 2x Coincidences");
            TEPXClustersPerEvent[d][r]->GetYaxis()->SetRangeUser(0, 6);
            TEPXClustersPerEvent[d][r]->GetXaxis()->SetRangeUser(0, 2);
            TEPXClustersPerEvent[d][r]->Draw("ape");
            FitTEPXClustersPerEvent[d][r]->Draw("lsame");
            label.SetTextSize(0.1);
            label.DrawLatexNDC(0.15, 0.8, TString("D=") + (long)(d <= 3 ? d - 4 : d + 1 - 4) + ",  R=" + (long)(r + 1));
            char* histname = new char[10];
            sprintf(histname, "histo%d_linearity2.gif", l);
            cout << "==========================" << histname << endl;
            C.Print(outputpath + histname);
            C.Clear();
}
            
            if (option == 3){
            
            TEPXClustersPerEvent[d][r]->GetYaxis()->SetTitle("Mean Number of 3x Coincidences");
            TEPXClustersPerEvent[d][r]->GetYaxis()->SetRangeUser(0, 6);
            TEPXClustersPerEvent[d][r]->GetXaxis()->SetRangeUser(0, 2);
            TEPXClustersPerEvent[d][r]->Draw("ape");
            FitTEPXClustersPerEvent[d][r]->Draw("lsame");
            label.SetTextSize(0.1);
            label.DrawLatexNDC(0.15, 0.8, TString("D=") + (long)(d <= 3 ? d - 4 : d + 1 - 4) + ",  R=" + (long)(r + 1));
            char* histname = new char[10];
            sprintf(histname, "histo%d_linearity2.gif", l);
            cout << "==========================" << histname << endl;
            C.Print(outputpath + histname);
            C.Clear();
            
            
 
     }      
}
    
    
    //Set residuals points and residuals point error
    for (int d = 0; d < disklist.size(); d++) {
        for (int r = 0; r < 5; r++) {
            for (int pu = 0; pu < pulist.size(); pu++) {


	        float x = pumap[pulist[pu]];
                float y = (Prof_TEPXClustersPerEvent[pu][d]->GetBinContent(r + 1) - FitTEPXClustersPerEvent[d][r]->Eval(pumap[pulist[pu]])) / FitTEPXClustersPerEvent[d][r]->Eval(pumap[pulist[pu]]);
                float yerr = (Prof_TEPXClustersPerEvent[pu][d]->GetBinError(r + 1)) / FitTEPXClustersPerEvent[d][r]->Eval(pumap[pulist[pu]]);

             
                if (pu > 2) {


                    NonLinearity_TEPXClustersPerEvent[d][r]->SetPoint(pu, x -(r-2), y);
                    NonLinearity_TEPXClustersPerEvent[d][r]->SetPointError(pu, 0, yerr);
                } 
                else 
                {


		    NonLinearity_TEPXClustersPerEvent[d][r]->SetPoint(pu, x, y);
                    NonLinearity_TEPXClustersPerEvent[d][r]->SetPointError(pu, 0, yerr);
                }

                
            }
        }
    }


    TCanvas C3("C3");
    C3.cd();

    //draw residuals vs pileup graph for Disk 4 Ring 1
    if (option == 1){
    NonLinearity_TEPXClustersPerEvent[7][0]->GetYaxis()->SetNdivisions(10);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetYaxis()->SetLabelSize(0.03);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetYaxis()->SetTitle("Residual   (Data-Fit)/Fit");
    NonLinearity_TEPXClustersPerEvent[7][0]->GetYaxis()->SetMaxDigits(4);
    NonLinearity_TEPXClustersPerEvent[7][0]->SetLineColor(6);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetYaxis()->SetRangeUser(-0.015, 0.015);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetXaxis()->SetNdivisions(12);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetXaxis()->SetLabelSize(0.04);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetXaxis()->SetRangeUser(0, 210);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetXaxis()->SetTitle("Pileup");
    NonLinearity_TEPXClustersPerEvent[7][0]->SetMarkerStyle(21);
    NonLinearity_TEPXClustersPerEvent[7][0]->SetMarkerColor(2);
    NonLinearity_TEPXClustersPerEvent[7][0]->SetMarkerSize(0.001);
    NonLinearity_TEPXClustersPerEvent[7][0]->SetTitle("Disk 4 Ring 1");
    TLine* line1 = new TLine(0, 0, 210, 0);
    line1->SetLineColor(kBlack);

    gPad->SetGrid(1, 1);
    NonLinearity_TEPXClustersPerEvent[7][0]->Draw("ape");
    line1->Draw("same");

    C3.Print(outputpath + "NS1Brane.gif");
    C3.Clear();
    }

    if (option == 2){
    NonLinearity_TEPXClustersPerEvent[7][0]->GetYaxis()->SetNdivisions(10);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetYaxis()->SetLabelSize(0.03);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetYaxis()->SetTitle("Residual   (Data-Fit)/Fit");
    NonLinearity_TEPXClustersPerEvent[7][0]->GetYaxis()->SetMaxDigits(4);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetYaxis()->SetRangeUser(-0.15, 0.15);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetXaxis()->SetNdivisions(12);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetXaxis()->SetLabelSize(0.04);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetXaxis()->SetRangeUser(0, 210);
    NonLinearity_TEPXClustersPerEvent[7][0]->GetXaxis()->SetTitle("Pileup");
    NonLinearity_TEPXClustersPerEvent[7][0]->SetMarkerStyle(21);
    NonLinearity_TEPXClustersPerEvent[7][0]->SetMarkerColor(2);
    NonLinearity_TEPXClustersPerEvent[7][0]->SetMarkerSize(0.6);
    NonLinearity_TEPXClustersPerEvent[7][0]->SetTitle("Disk 4 Ring 1");
    TLine* line1 = new TLine(0, 0, 210, 0);
    line1->SetLineColor(kBlack);

    gPad->SetGrid(1, 1);
    NonLinearity_TEPXClustersPerEvent[7][0]->Draw("ape");
    line1->Draw("same");

    C3.Print(outputpath + "NS2Brane.gif");
    C3.Clear();
    }
    if (option == 3){
    NonLinearity_TEPXClustersPerEvent[7][1]->GetYaxis()->SetNdivisions(10);
    NonLinearity_TEPXClustersPerEvent[7][1]->GetYaxis()->SetLabelSize(0.03);
    NonLinearity_TEPXClustersPerEvent[7][1]->GetYaxis()->SetTitle("Residual   (Data-Fit)/Fit");
    NonLinearity_TEPXClustersPerEvent[7][1]->GetYaxis()->SetMaxDigits(4);
    NonLinearity_TEPXClustersPerEvent[7][1]->GetYaxis()->SetRangeUser(-5, 5);
    NonLinearity_TEPXClustersPerEvent[7][1]->GetXaxis()->SetNdivisions(12);
    NonLinearity_TEPXClustersPerEvent[7][1]->GetXaxis()->SetLabelSize(0.04);
    NonLinearity_TEPXClustersPerEvent[7][1]->GetXaxis()->SetRangeUser(0, 210);
    NonLinearity_TEPXClustersPerEvent[7][1]->GetXaxis()->SetTitle("Pileup");
    NonLinearity_TEPXClustersPerEvent[7][1]->SetMarkerStyle(21);
    NonLinearity_TEPXClustersPerEvent[7][1]->SetMarkerColor(2);
    NonLinearity_TEPXClustersPerEvent[7][1]->SetMarkerSize(0.6);
    NonLinearity_TEPXClustersPerEvent[7][1]->SetTitle("Disk 4 Ring 2");
    TLine* line1 = new TLine(0, 0, 210, 0);
    line1->SetLineColor(kBlack);

    gPad->SetGrid(1, 1);
    NonLinearity_TEPXClustersPerEvent[7][1]->Draw("ape");
    line1->Draw("same");

    C3.Print(outputpath + "NS3Brane.gif");
    C3.Clear();
    }
    
    TCanvas C2("C2");
    C2.cd();
    

    //draw residuals vs pileup graph
    if(option == 1){
    for (int d = 0; d < disklist.size(); d++) {

        NonLinearity_TEPXClustersPerEvent[d][0]->GetYaxis()->SetNdivisions(10);
        NonLinearity_TEPXClustersPerEvent[d][0]->GetYaxis()->SetLabelSize(0.02);
        NonLinearity_TEPXClustersPerEvent[d][0]->GetYaxis()->SetTitle("Residual   (Data-Fit)/Fit");
        NonLinearity_TEPXClustersPerEvent[d][0]->GetYaxis()->SetRangeUser(-0.015, 0.015);

        NonLinearity_TEPXClustersPerEvent[d][0]->GetXaxis()->SetNdivisions(12);
        NonLinearity_TEPXClustersPerEvent[d][0]->GetXaxis()->SetLabelSize(0.04);
        NonLinearity_TEPXClustersPerEvent[d][0]->GetXaxis()->SetRangeUser(0, 210);
        NonLinearity_TEPXClustersPerEvent[d][0]->GetXaxis()->SetTitle("Pileup");

        NonLinearity_TEPXClustersPerEvent[d][0]->SetMarkerStyle(20);
        NonLinearity_TEPXClustersPerEvent[d][0]->SetLineColor(1);
        NonLinearity_TEPXClustersPerEvent[d][0]->SetMarkerSize(0);
        if (d <= 3) {
            NonLinearity_TEPXClustersPerEvent[d][0]->SetTitle(TString("Disk ") + (d - 4));
        } else {
            NonLinearity_TEPXClustersPerEvent[d][0]->SetTitle(TString("Disk ") + (d - 3));
        }

        auto legend = new TLegend(0.7, 0.2, 0.8, 0.4);
        NonLinearity_TEPXClustersPerEvent[d][0]->Draw("ape");
        legend->AddEntry(NonLinearity_TEPXClustersPerEvent[d][0], TString(" Ring ") + 1, "l");

        for (int r = 1; r < 5; r++) {

            NonLinearity_TEPXClustersPerEvent[d][r]->SetMarkerStyle(20);
            NonLinearity_TEPXClustersPerEvent[d][r]->SetLineColor(r + 1);
            NonLinearity_TEPXClustersPerEvent[d][r]->SetMarkerSize(0);
            NonLinearity_TEPXClustersPerEvent[d][r]->Draw("pesame");
            legend->AddEntry(NonLinearity_TEPXClustersPerEvent[d][r], TString(" Ring ") + (r + 1), "l");
        }

        TLine* line = new TLine(0, 0, 220, 0);
        line->SetLineColor(kBlack);
        line->Draw("same");

        legend->SetFillColor(0);
        legend->SetLineColor(0);
        legend->SetFillColor(0);
        legend->Draw("same");

        C2.Print(outputpath + TString("Clusterdisk_") + d + ".gif");
    }
}
    
   if(option == 2){
    	for (int d = 0; d < disklist.size(); d++) {

        NonLinearity_TEPXClustersPerEvent[d][0]->GetYaxis()->SetNdivisions(10);
        NonLinearity_TEPXClustersPerEvent[d][0]->GetYaxis()->SetLabelSize(0.02);
        NonLinearity_TEPXClustersPerEvent[d][0]->GetYaxis()->SetTitle("Residual   (Data-Fit)/Fit");
        NonLinearity_TEPXClustersPerEvent[d][0]->GetYaxis()->SetRangeUser(-0.15, 0.15);

        NonLinearity_TEPXClustersPerEvent[d][0]->GetXaxis()->SetNdivisions(12);
        NonLinearity_TEPXClustersPerEvent[d][0]->GetXaxis()->SetLabelSize(0.04);
        NonLinearity_TEPXClustersPerEvent[d][0]->GetXaxis()->SetRangeUser(0, 210);
        NonLinearity_TEPXClustersPerEvent[d][0]->GetXaxis()->SetTitle("Pileup");

        NonLinearity_TEPXClustersPerEvent[d][0]->SetMarkerStyle(20);
        NonLinearity_TEPXClustersPerEvent[d][0]->SetLineColor(1);
        NonLinearity_TEPXClustersPerEvent[d][0]->SetMarkerSize(0);

        if (d <= 3) {
            NonLinearity_TEPXClustersPerEvent[d][0]->SetTitle(TString("Disk ") + (d - 4));
        } else {
            NonLinearity_TEPXClustersPerEvent[d][0]->SetTitle(TString("Disk ") + (d - 3));
        }

        auto legend = new TLegend(0.7, 0.2, 0.8, 0.4);
        NonLinearity_TEPXClustersPerEvent[d][0]->Draw("ape");
        legend->AddEntry(NonLinearity_TEPXClustersPerEvent[d][0], TString(" Ring ") + 1, "l");

        for (int r = 1; r < 5; r++) {

            NonLinearity_TEPXClustersPerEvent[d][r]->SetMarkerStyle(20);
            NonLinearity_TEPXClustersPerEvent[d][r]->SetLineColor(r + 1);
            NonLinearity_TEPXClustersPerEvent[d][r]->SetMarkerSize(0);
            NonLinearity_TEPXClustersPerEvent[d][r]->Draw("pesame");
            legend->AddEntry(NonLinearity_TEPXClustersPerEvent[d][r], TString(" Ring ") + (r + 1), "l");
        }

        TLine* line = new TLine(0, 0, 220, 0);
        line->SetLineColor(kBlack);
        line->Draw("same");

        legend->SetFillColor(0);
        legend->SetLineColor(0);
        legend->SetFillColor(0);
        legend->Draw("same");

        C2.Print(outputpath + TString("2xCoincidencesdisk_") + d + ".gif");
    }
} 

  
    if(option == 3){
    	for (int d = 0; d < disklist.size(); d++) {

        NonLinearity_TEPXClustersPerEvent[d][1]->GetYaxis()->SetNdivisions(10);
        NonLinearity_TEPXClustersPerEvent[d][1]->GetYaxis()->SetLabelSize(0.02);
        NonLinearity_TEPXClustersPerEvent[d][1]->GetYaxis()->SetTitle("Residual   (Data-Fit)/Fit");
        NonLinearity_TEPXClustersPerEvent[d][1]->GetYaxis()->SetRangeUser(-50, 50);

        NonLinearity_TEPXClustersPerEvent[d][1]->GetXaxis()->SetNdivisions(12);
        NonLinearity_TEPXClustersPerEvent[d][1]->GetXaxis()->SetLabelSize(0.04);
        NonLinearity_TEPXClustersPerEvent[d][1]->GetXaxis()->SetRangeUser(0, 210);
        NonLinearity_TEPXClustersPerEvent[d][1]->GetXaxis()->SetTitle("Pileup");

        NonLinearity_TEPXClustersPerEvent[d][1]->SetMarkerStyle(20);
        NonLinearity_TEPXClustersPerEvent[d][1]->SetLineColor(1);
        NonLinearity_TEPXClustersPerEvent[d][1]->SetMarkerSize(0);
        if (d <= 3) {
            NonLinearity_TEPXClustersPerEvent[d][1]->SetTitle(TString("Disk ") + (d - 4));
        } else {
            NonLinearity_TEPXClustersPerEvent[d][1]->SetTitle(TString("Disk ") + (d - 3));
        }

        auto legend = new TLegend(0.7, 0.2, 0.8, 0.4);
        NonLinearity_TEPXClustersPerEvent[d][1]->Draw("ape");
        legend->AddEntry(NonLinearity_TEPXClustersPerEvent[d][1], TString(" Ring ") + 2, "l");

        for (int r = 2; r < 5; r++) {

            NonLinearity_TEPXClustersPerEvent[d][r]->SetMarkerStyle(20);
            NonLinearity_TEPXClustersPerEvent[d][r]->SetLineColor(r + 1);
            NonLinearity_TEPXClustersPerEvent[d][r]->SetMarkerSize(0);
            NonLinearity_TEPXClustersPerEvent[d][r]->Draw("pesame");
            legend->AddEntry(NonLinearity_TEPXClustersPerEvent[d][r], TString(" Ring ") + (r + 1) , "l");
        }

        TLine* line = new TLine(0, 0, 220, 0);
        line->SetLineColor(kBlack);
        line->Draw("same");

        legend->SetFillColor(0);
        legend->SetLineColor(0);
        legend->SetFillColor(0);
        legend->Draw("same");

        C2.Print(outputpath + TString("3xCoincidencesdisk_") + d + ".gif");
    }
  }

   
}

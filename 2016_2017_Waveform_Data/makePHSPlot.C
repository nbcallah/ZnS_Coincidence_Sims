void doit() {
	TTree* tFG = new TTree("nameFG", "title");
	tFG->ReadFile("./PHSTextHist_tree.csv", "num/D:num1:num2");
	
	TTree* tBKG = new TTree("nameBKG", "title");
	tBKG->ReadFile("./PHSTextHist_tree_bkg.csv", "num/D:num1:num2");
	
	int numFG = tFG->Draw("num1:num2:num", "");
	int numBKG = tBKG->Draw("num1:num2:num", "");
	
	TH2D* phsFG = new TH2D("phsFG", "Foreground Pulse Height Spectrum", 100, 0, 100, 100, 0, 100);
	TH2D* phsBKG = new TH2D("phsFG", "Background Pulse Height Spectrum", 100, 0, 100, 100, 0, 100);
	
	for(int i = 0; i < numFG; i++) {
		phsFG->Fill(tFG->GetV1()[i], tFG->GetV2()[i], tFG->GetV3()[i]);
	}
	
	for(int i = 0; i < numBKG; i++) {
		phsBKG->Fill(tBKG->GetV1()[i], tBKG->GetV2()[i], tBKG->GetV3()[i]);
	}
	
	phsFG->GetXaxis()->SetLabelSize(0.05);
	phsFG->GetXaxis()->SetTitleSize(0.05);
	phsFG->GetXaxis()->SetTitle("PMT 1 Sum");
	phsFG->GetYaxis()->SetLabelSize(0.05);
	phsFG->GetYaxis()->SetTitleSize(0.05);
	phsFG->GetYaxis()->SetTitle("PMT 2 Sum");
	phsFG->SetStats(0);
	phsFG->GetZaxis()->SetLabelSize(0.05);
	
	phsBKG->GetXaxis()->SetLabelSize(0.05);
	phsBKG->GetXaxis()->SetTitleSize(0.05);
	phsBKG->GetXaxis()->SetTitle("PMT 1 Sum");
	phsBKG->GetYaxis()->SetLabelSize(0.05);
	phsBKG->GetYaxis()->SetTitleSize(0.05);
	phsBKG->GetYaxis()->SetTitle("PMT 2 Sum");
	phsBKG->SetStats(0);
	phsBKG->GetZaxis()->SetLabelSize(0.05);
	
	TCanvas* can = new TCanvas("can", "can2");
	//printf("%f\n", can->GetRightMargin());
	can->SetRightMargin(0.15);
	
	phsFG->Draw("COLZ");
	can->SaveAs("fgPHS.pdf");
	
	can->SetLogz();
	
	phsBKG->Draw("COLZ");
	can->SaveAs("bkgPHS.pdf");
}
//
// arXiv:1509.06176v1, Una parametrización del flujo de muones de rayos cósmicos a nivel del mar
// La ecuación (3) es para el cálculo del flujo de muones de 0,5 a 6000 GeV, para un ángulo cenital de 0 y 75 grados.
// Se traza la unidad de flujo sin escalar (m^{-2} s^{-1} sr^{-1} GeV^{-1}).

// cálculo del flujo Gaisser con correción
double getFlux(double muonEnergy, double theta)
{
	const double P1 = 0.102573;
	const double P2 = -0.068287;
	const double P3 = 0.958633;
	const double P4 = 0.0407253;
	const double P5 = 0.817285;

	double cosTheta = cos(theta); 
	double cosThetaStar2 = ( cosTheta * cosTheta + P1 * P1 + P2 * pow(cosTheta, P3) + P4 * pow(cosTheta, P5) ) /
		(1. + P1*P1 + P2 + P4);
	double cosThetaStar = sqrt(cosThetaStar2);

	double Emu = muonEnergy;       // GeV

	double term1 = 1./(1 + 1.1 * Emu * cosThetaStar/115);
	double term2 = 0.054/(1 + 1.1 * Emu * cosThetaStar/850);

	double flux = 0.14 * pow(Emu * (1. + 3.64/(Emu * pow(cosThetaStar, 1.29))), -2.7) * (term1 +  term2);

	return flux * 10000; // 1/cm^2 -> 1/m^2
}

double funFlux(double *x, double *par) {
	double muonEnergy = x[0];
	double theta = par[0]/180. * TMath::Pi();

	return getFlux(muonEnergy, theta);
}

// función main ROOT
void muonFlux2() {

   TCanvas* c1 = new TCanvas("c1", "  "); // llamado del Canvas
   c1->SetLogx();                         // escala LOG x
   c1->SetLogy();                         // escala LOG y
   c1->SetGrid();                         // rejilla

   TF1 *f1 = new TF1("f0", funFlux, 0.5, 6000, 1); // función tipo F1
   f1->SetParameter(0, 0.);                        // params.
   f1->GetYaxis()->SetRangeUser(0.00000001, 300);
   f1->SetLineColor(kRed);
   f1->SetLineWidth(3);
   f1->GetXaxis()->SetTitle("E_{#mu}(GeV)");
   f1->GetYaxis()->SetTitle("I_{#mu} (m^{-2} s^{-1} sr^{-1} GeV^{-1})");
   f1->SetTitle("");
   f1->Draw();

   TF1 *f2 = new TF1("f75", funFlux, 0.5, 6000, 1);
   f2->SetParameter(0, 75.);
   f2->SetLineColor(kBlue);
   f2->SetLineWidth(3);
   f2->Draw("same");

   TLegend* legend = new TLegend(0.75, 0.78, 0.9, 0.9);
   legend->AddEntry("f0","#theta = 0#circ", "l");
   legend->AddEntry("f75","#theta = 75#circ", "l");
   legend->Draw();

   c1->Print("muonFlux2.pdf");
}
//
// arXiv:1509.06176v1, Una parametrización del flujo de muones de rayos cósmicos al nivel del mar, ecuación (3).
// Se representa gráficamente la intensidad integral de los muones por encima de 1 GeV/c en función del ángulo cenital.
//
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
	double muonEnerg = x[0];
	double theta = par[0]/180. * TMath::Pi();

	return getFlux(muonEnergy, theta);
}

// código main de ROOT
// el proceso general calcula el flujo con la forma de Gaisser, con correción esférica
// y después integra 
void muonFlux3() {
   TCanvas* c1 = new TCanvas("c1", "  ");
   c1->SetGrid();

   const int nAngle = 90;
   double angles[nAngle];
   double fluxes[nAngle];

   for(int i = 0; i < nAngle; ++i) {
       TF1 f1("f0", funFlux, 0.2, 6000, 1);  // crea el tipo de función F1
       f1.SetParameter(0, i);                // ajusta el param. 0 como el i del FOR
       ROOT::Math::WrappedTF1 wf1(f1);       // crea la función de peso WF1
       ROOT::Math::GaussIntegrator ig;       // llama un integrador Gauss (mtd de integración)
       ig.SetFunction(wf1);                  // ajusta WF1 al integrador 
       ig.SetRelTolerance(0.001);            // tolerancia 

	   angles[i] = i;                         // ángulo θ

	   // rango de energía para la integración del flujo de muones
       const double El = 1.;           // lim. inferior Ene
       const double Eh = 5000.;        // lim. superior Ene
	   fluxes[i] = ig.Integral(El, Eh); // flujo integrado 
   }

   // gráfica del flujo inttegrado 𝑓(Ene)
   TGraph* gr = new TGraph(nAngle, angles, fluxes);
   gr->SetTitle("E_{#mu} > 1 GeV");
   gr->GetXaxis()->SetRangeUser(0, 89);
   gr->GetXaxis()->SetTitle("Zenith angle (#circ)");
   gr->GetYaxis()->SetTitle("Muon flux (m^{-2} s^{-1} sr^{-1})");
   gr->SetLineColor(kRed);
   gr->SetLineWidth(3);
   gr->SetMarkerColor(kRed);
   gr->Draw("ACP");

   c1->Print("muonFlux3.pdf");
}
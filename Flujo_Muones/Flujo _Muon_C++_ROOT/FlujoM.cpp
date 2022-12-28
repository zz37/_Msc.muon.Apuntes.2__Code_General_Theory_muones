#include <cmath>
#include <iostream>
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"


// prototipos
double muonFlujo(double muonEnergia, double theta);
double flujoTot(double *x, double *par);

int main() {
    TCanvas* c1 = new TCanvas("c1", "  ");
    c1->SetLogx(1);
    c1->SetLogy(1);
    c1->SetGrid();
    TF1 *f1 = new TF1("Fórmula de Gaisser", flujoTot, 1, 1000, 1);
	f1->GetXaxis()->SetTitle("E_{#mu}(GeV)");
   	f1->GetYaxis()->SetTitle("#Phi(E_{#mu},#theta) [cm^{-2} s^{-1} sr^{-1} (GeV)^{-1}]");
    f1->Draw();
	c1->SaveAs("muonFlujo.png");

    return 0;
}

// Formula de Gaisser
// REF :rpp2020-rev-cosmic-rays, 30.3.1 Muons
// Eq.(30.4), válido para Eµ > 100 GeV y ángulo cenital, θ > 70°
double muonFlujo(double muonEnergia, double theta) {
	double Emu = muonEnergia;      // en GeV
	double costh = cos(theta);     // rads
 	double pion = 1./(1 + 1.1 * Emu * costh / 115);  // parte del pion
	double kaon = 0.054/(1 + 1.1 * Emu * costh/850); // parte del kaon
	double flujo = 0.14 * (pion +  kaon);            // flujo
	return flujo; 
}

double flujoTot(double *x, double *par) {
	double muonEnergia = x[0]; // energía muon
	// theta = 0
	return muonFlujo(muonEnergia, 0);
}
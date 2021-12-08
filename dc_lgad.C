#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include <TTree.h>
#include <TFile.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

double transfer(double t) {
	constexpr double tauR = 5.6;
	constexpr double tauI = 1.8;
	constexpr double tauA = 47.;
	constexpr double dAI = (tauA - tauI);
	constexpr double dAR = (tauA - tauR);
	constexpr double dIR = (tauI - tauR);
	constexpr double c1 = tauA / (dAI * dAI *dAR);
	constexpr double c2 = 1. / (dAI * tauI * dIR);
	constexpr double c3 = tauR / (dAR * dIR * dIR);
	constexpr double c4 = (tauI * tauI * tauI - tauA * tauI * tauR) / (dAI * dAI * tauI * dIR * dIR);
	const double f1 = -exp(-t / tauA) * c1;
	const double f2 = exp(-t / tauI) * t * c2;
	const double f3 = exp(-t / tauR) * c3;
	const double f4 = exp(-t / tauI) * c4;
	return tauA * tauR * (f1 + f2 + f3 + f4);

}

int main(int argc, char * argv[]) {
	
	TApplication app("app", &argc, argv);

	const double depth = 50.e-4;
	const double left = 175.e-4;
	const double right = 1125.e-4;
	const double width = 650.e-4;

	MediumSilicon si;
	ComponentTcad2d cmp;
//  cmp.EnableDebugging();
  cmp.Initialise("n3397_des.grd", "n3397_des.dat");
	cmp.SetRangeZ(-width, width);

	int nRegions = cmp.GetNumberOfRegions();
	for (int i = 0; i < nRegions; ++i) {
		std::string region;
		bool active;
		cmp.GetRegion(i, region, active);
		if (region == "\"Silicon_1\"") cmp.SetMedium(i, &si);
	}
	
/*
	ComponentAnalyticField wfieldAnalytic;
	wfieldAnalytic.AddPlaneY(0., 0., "top");
	wfieldAnalytic.AddPlaneY(depth, -500., "bot");
	wfieldAnalytic.AddStripOnPlaneY('z', depth, 0., width, "strip");
	wfieldAnalytic.AddReadout("strip");
*/
	ComponentAnalyticField wfieldAnalytic;
	wfieldAnalytic.AddPlaneX(0., 0., "top");
	wfieldAnalytic.AddPlaneX(depth, -100, "bot");
	wfieldAnalytic.AddStripOnPlaneX('z', 0., left, right, "strip");
	wfieldAnalytic.AddReadout("strip");

	ViewField vField;
	constexpr bool plotField = false;
	if (plotField) {
		vField.SetComponent(&cmp);
		vField.PlotContour("v");
	}
  constexpr bool plotWeightingField = false;
  if (plotWeightingField) {
    ViewField* wfieldView = new ViewField();
    wfieldView->SetComponent(&wfieldAnalytic);
    wfieldView->SetArea(0, 0, depth, 1300.e-4);
		wfieldView->PlotContourWeightingField("strip", "v");
  }

	Sensor sensor;
	sensor.AddComponent(&cmp);
	sensor.AddElectrode(&wfieldAnalytic, "strip");
	
	const std::string label = "strip";

	const double endtime = 20;
//	const double endtime = 5;
	const int nSignalBins = 100000;
	double tStep = endtime / nSignalBins;
	sensor.SetTimeWindow(0., tStep, nSignalBins);
	sensor.SetTransferFunction(transfer);

	const double thr1 = -1000 * ElementaryCharge;
	std::cout << "Threshold: " << thr1 << " fC\n";

	AvalancheMC drift;
	drift.SetDistanceSteps(0.1e-4);
	drift.SetSensor(&sensor);
	drift.EnableSignalCalculation();
	drift.EnableInducedChargeCalculation();

	TrackHeed track;
	track.SetSensor(&sensor);
	track.SetParticle("pi");
	track.SetMomentum(420.e6);

	ViewSignal vSignal;
	constexpr bool plotSignal = false;
	vSignal.SetSensor(&sensor);
	
	constexpr bool writeSignal = false;
	constexpr bool writeTime = true;

	ViewDrift vDrift;
	constexpr bool plotDrift = false;
	if (plotDrift) {
		vDrift.SetArea(0., 0., depth, 1300.e-4);
		track.EnablePlotting(&vDrift);
	}

	std::ofstream output_charge, output_time, output_gain;
	std::ofstream log;
	output_charge.open("charge_lgad.dat", std::ios::out);
	output_time.open("time_lgad.dat", std::ios::out);
	output_gain.open("gain.dat", std::ios::out);
	log.open("log_lgad", std::ios::out);

	const unsigned int nEvents = 500;
	for (unsigned int j = 0; j < nEvents; ++j) {
		sensor.ClearSignal();
		double y0 = 0.5 * (left + right);
		y0 += (RndmUniform() - 0.5) * (right - left) * 0.5; 
		double t0 = 0.1;
		track.NewTrack(0., y0, 0., t0, 0.1, 0., 0.);
		double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., dummy = 0.;
		int ne = 0, nh = 0;
		unsigned int nc = 0;
		unsigned int nesum = 0, nesum_real = 0, doublecheck = 0;
		std::cout << "No." << j << "\n";
		while (track.GetCluster(xc, yc, zc, tc, ne, nh, ec, dummy)) {
			++nc;
			nesum += ne;
			nesum_real += ne;
			if (nc % 100 == 0) std::cout << "		Cluster " << nc << "\n";
			drift.DisablePlotting();
			if (plotDrift && RndmUniform() < 0.05) {
				drift.EnablePlotting(&vDrift);
			}
			double xe, ye, ze, te, ee, dxe, dye, dze;
			for (int i = 0; i < ne; ++i) {
				track.GetElectron(i, xe, ye, ze, te, ee, dxe, dye, dze);
//				drift.DriftElectron(xe, ye, ze, te);
				drift.AvalancheElectronHole(xe, ye, ze, te);
				int size = drift.GetNumberOfElectronEndpoints();
				log << "size = " << size << "\n";
				nesum_real += size - 1;
				doublecheck += size;
				if (size > 1) log << "avalanche with " << size << "\n";
			}
			double xh, yh, zh, th;
/*			for (int i = 0; i < nh; ++i) {
				track.GetIon(i, xh, yh, zh, th);
				drift.DriftHole(xh, yh, zh, th);
			}*/
		}
		std::cout << nesum << " electrons, " << nesum * ElementaryCharge << " fC.\n";
//		std::cout << nesum_real << " electrons, " << nesum_real * ElementaryCharge << " fC.\n";
		std::cout << doublecheck << " electrons, " << doublecheck * ElementaryCharge << " fC.\n";
		
		double gain = (double)nesum_real / nesum;
		output_gain << gain << "\n";

		output_charge << nesum_real * ElementaryCharge << "\n";
		sensor.ConvoluteSignals();

		if (writeSignal) {
			char filename[50];
			sprintf(filename, "signal_lgad_%05d.txt", j);
			std::ofstream outfile;
			outfile.open(filename, std::ios::out);
			for (unsigned int i = 0; i < nSignalBins; ++i) {
				const double t = (i + 0.5) * tStep;
				const double f = sensor.GetSignal(label, i);
				const double fe = sensor.GetElectronSignal(label, i);
				const double fh = sensor.GetIonSignal(label, i);
				outfile << t << "  " << f << "	" << fe << "	" << fh << "\n";
			}
			outfile.close();
		}

		if (writeTime) {
			std::cout << "write time" << "\n";
			double min = 999;
			double time = 0;
			double min_time = 0;
			for (unsigned int i = 0; i < nSignalBins; ++i) {
				const double f = sensor.GetSignal(label, i);
				const double t = (i + 0.5) * tStep;
				if (f < min) {
					min = f;
					min_time = t;
				}
			}
			for (unsigned int i = 0; i < nSignalBins; ++i) {
				const double f = sensor.GetSignal(label, i);
				const double t = (i + 0.5) * tStep;
				if (f > 0.5 * min && t < min_time) time = t;
			}
			output_time << time << "	" << min << "\n";
		}
  /*
	AvalancheMicroscopic aval;
  aval.SetSensor(&sensor);

  for (unsigned int i = 0; i < nEvents; ++i) { 
    std::cout << i << "/" << nEvents << "\n";
    // Randomize the initial position. 
    const double x0 = 0;
    const double y0 = 650.e-4;
    const double z0 = 0.; 
    const double t0 = 0.;
    const double e0 = 0.1;
    aval.AvalancheElectron(x0, y0, z0, t0, e0, 0.1, 0., 0.);
    int ne = 0, ni = 0;
    aval.GetAvalancheSize(ne, ni);
    const unsigned int np = aval.GetNumberOfElectronEndpoints();
    double xe1, ye1, ze1, te1, e1;
    double xe2, ye2, ze2, te2, e2;
    double xi1, yi1, zi1, ti1;
    double xi2, yi2, zi2, ti2;
    int status;
    for (unsigned int j = 0; j < np; ++j) {
      aval.GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, 
                                  xe2, ye2, ze2, te2, e2, status);
      drift.DriftHole(xe1, ye1, ze1, te1);
      drift.GetHoleEndpoint(0, xi1, yi1, zi1, ti1, 
                              xi2, yi2, zi2, ti2, status);
    }
		if (!plotSignal) continue;
		vSignal.PlotSignal("strip");
  }
*/
//		if (!plotSignal) continue;
//		vSignal.PlotSignal("strip");
	}
	
	output_charge.close();
	output_time.close();
	output_gain.close();
	log.close();

	if (plotDrift) {
		const bool twod = true;
		vDrift.Plot(twod, true);
	}

//	app.Run(true);

}

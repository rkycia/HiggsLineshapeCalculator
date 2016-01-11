/***********************************************************************
 
	@brief ISR integrator based on Foam - MPI parallel version

	@authors Stanislaw Jadach, Radoslaw A. Kycia

***********************************************************************/

//C headers
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// C++ headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <string>

using namespace std;

// ROOT headers
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TFoam.h>
#include <TFoamIntegrand.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGaxis.h>


#include "TDensity.h"

//MPI support
#include <mpi.h>  


//level of details
//#define DEBUG
#undef DEBUG


#define ELECTRON
//#define MUON

////////////////////////////////////////////////////////////////////////


/// @brief Converts number to string
#define STR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()


////////////////////////////////////////////////////////////////////////

/// @brief Displays simple progress bar
///\arg [in] Now - how much is left
///\arg [in] Total - general number of events
///\arg [in] width  - how many dots will be used in display
int progressBar( double Now, double Total, int width = 40 )
{
  
    double fraction = Now / Total;
    // number of = to display
    int N = round(fraction * width);

    
    int i=0;
    printf("%3.0f%% [",fraction*100);
    // done part
    for ( ; i < N;i++) 
    {
        printf("=");
    }
    // remaining part (spaces)
    for ( ; i < width;i++) 
    {
        printf(" ");
    }
    //print end and info
    printf("] ");
    printf(" %.0f of %.0f", Now, Total );
    //go back at the beginning
    printf("\r");
    fflush(stdout);
    
    return(0);
    
};


/// @brief normalize given histogram using norm histogram
/// @param NorHst  histogram with norm
/// @param Hst histogram to normalize
void HistNorm(TH1D *NorHst, TH1D *Hst)
{
	Hst->ls();
	Double_t Nevt = NorHst->GetBinContent(2);
	Double_t Xsav = NorHst->GetBinContent(1)/Nevt;
	//
	int      nbt  = Hst->GetNbinsX();
	Double_t tmax = Hst->GetXaxis()->GetXmax();
	Double_t tmin = Hst->GetXaxis()->GetXmin();
	Double_t Fact = nbt*Xsav/(tmax-tmin)/Nevt;
	//cout<<"HistNorm: Xsav = "<<Xsav<<"  Nevt =  "<<Nevt<<"  Fact = "<<Fact<<endl;
	Hst->Scale(Fact);
	
};
//----------------------------------------------------------------------

/// @brief Make Monte Carlo integration of Born convolution
/// @param filename filename of root and pdf files to write fistograms
/// @param keyISR   ISR type (a) - 0; (b) - 1; (c) - 2;
/// @param kDim     integration dimension: 2,3 Machine energy spread OFF/ON
/// @param sigE     energy spread if applicable
/// @param NevTot   statistics of integration
void MakeConvBorn( string filename = "histo1", Double_t sigE = 0.0041, long NevTot =  1000000 ) 
{
	
	// total dimension
	Int_t  kDim = 2;   
	
	//allocate density and set up
	TDensity * rho= new TDensity(); 
	//set convolution 2/3 - OFF/ON
	rho->m_kDim = 2;
	//turn off ISR, pure Born
	rho->m_ISROn = 0;
	//set energy spread
	rho->m_sigE = sigE;
	
	
	Double_t MCresult,MCerror;
	//=========================================================
	Int_t  nBin     =       8;   // Number of bins in build-up
	Int_t  OptRej   =       0;   // Wted events for OptRej=0; wt=1 for OptRej=1 (default)
	Int_t  OptDrive =       2;   // (D=2) Option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
	Int_t  EvPerBin =      25;   // Maximum events (equiv.) per bin in buid-up

	#ifdef DEBUG
		Int_t  Chat     =       1;   // Chat level
	#else
		Int_t  Chat     =       0;   // Chat level
	#endif	
		
	//=========================================================
	TRandom *PseRan   = new TRandom3();  // Create random number generator
	TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
	PseRan->SetSeed(4357);
	//=========================================================
	#ifdef DEBUG
		cout<<"*****   Foam version "<<FoamX->GetVersion()<<"    *****"<<endl;
	#else 
		//do nothing
	#endif	
	FoamX->SetkDim(        rho->m_kDim );      // Mandatory!!!
	FoamX->SetnCells(      rho->m_nCells);    // optional
	FoamX->SetnSampl(      rho->m_nSampl);    // optional
	FoamX->SetnBin(        nBin);      // optional
	FoamX->SetOptRej(      OptRej);    // optional
	FoamX->SetOptDrive(    OptDrive);  // optional
	FoamX->SetEvPerBin(    EvPerBin);  // optional
	FoamX->SetChat(        Chat);      // optional
	//===============================
	FoamX->SetRho(rho);
	FoamX->SetPseRan(PseRan);

	#ifdef DEBUG
		cout << "Initialze...." << endl;
	#else 
		//do nothing
	#endif	
	
	FoamX->Initialize(); // Initialize simulator
	
	#ifdef DEBUG
		cout << "...initialze DONE" << endl;
	#else 
		//do nothing
	#endif
	
	double Xnorm, errel;
	FoamX->GetIntNorm(Xnorm,errel);   // universal normalization
 
	//FoamX->Write("FoamX");     // Writing Foam on the disk, TESTING PERSISTENCY!!!

	long nCalls=FoamX->GetnCalls();


	Double_t *MCvect = new Double_t[kDim]; // vector generated in the MC run
	
	//Histograms
	int    nbx=200;       // linear scale
	double Emax = rho->m_MH +3*rho->m_GamH;
	double Emin = rho->m_MH -3*rho->m_GamH;
	
	TH1D * h_Ene = new TH1D("h_Ene", "h_Ene",nbx, Emin, Emax);
	h_Ene->Sumw2();

	//norm and no event - only two bins are used
	TH1D * h_NORM = new TH1D("h_NORM", "h_NORMA",2,0.0,2.0);
	h_NORM->Sumw2();	
	
	
	//event loop
	for(long loop=0; loop<NevTot; loop++)
	{
		//generate MC event
		FoamX->MakeEvent();
		FoamX->GetMCvect( MCvect);
		
		//get point info
		double E, MCwt;
		//double y;
		
		MCwt=FoamX->GetMCwt();
		E = rho->m_E;
		//y = rho->m_y;
    
		//filling histograms
		h_Ene->Fill( E, MCwt );
		
		//  Fill special normalization histogram hNORM
		h_NORM->Fill(0.5, Xnorm);   // 1-st bin = Normal*Nevtot
		h_NORM->Fill(1.5, 1);         // 2-nd bin = Nevtot
		
		#ifdef DEBUG
			//print progress bar
			if( loop % 1000 == 0 ) 
				progressBar( loop, NevTot );
		#else
			//do nothing
		#endif

	}

		#ifdef DEBUG
			cout << "====== Events generated, entering Finalize" << endl;
		#else
			//do nothing
		#endif

	Double_t eps = 0.0005;
	Double_t Effic, WtMax, AveWt, Sigma;
	Double_t IntNorm, Errel;
	FoamX->Finalize( IntNorm, Errel );     // final printout
	FoamX->GetIntegMC( MCresult, MCerror );  // get MC intnegral
	FoamX->GetWtParams( eps, AveWt, WtMax, Sigma ); // get MC wt parameters
	Effic=0; 
	
	if(WtMax>0) 
		Effic=AveWt/WtMax;

	#ifdef DEBUG		
		cout << "================================================================" << endl;
		cout << " MCresult= " << MCresult << " +- " << MCerror << " RelErr= "<< MCerror/MCresult << endl;
		cout << " Dispersion/<wt>= " << Sigma/AveWt << endl;
		cout << "      <wt>/WtMax= " << Effic <<",    for epsilon = "<<eps << endl;
		cout << " nCalls (initialization only) =   " << nCalls << endl;
		cout << "================================================================" << endl;
	#else
		//do nothing
	#endif
	
	
	delete [] MCvect;
	
	//make sigma histograms by normalization 
	TH1D *h_SigEne = (TH1D*)h_Ene->Clone("h_SigEne");
	HistNorm( h_NORM, h_SigEne );
	
	//remove statistics 
	h_SigEne->SetStats(0);
	
	
	//Save pdf
	Float_t  WidPix, HeiPix;
	WidPix = 800; HeiPix =  800;
	TCanvas *cCanv = new TCanvas("cCanv","cCanv", 100,100, WidPix,HeiPix);
	cCanv->SetFillColor(10);
	h_SigEne->Draw();
	cCanv->Update();
	cCanv->Print((filename+string(".pdf")).c_str(), "");
	delete cCanv;
	
	
	//Save ROOT
	TFile RootFile( (filename+string(".root")).c_str(), "RECREATE", "histograms");
	RootFile.ls();
	h_Ene->Write();
	h_SigEne->Write();	
	h_NORM->Write();
	RootFile.Close();

	//cleaning;
	delete FoamX;
	delete PseRan;
	delete rho;
	
	delete h_Ene;
	delete h_SigEne;
	delete h_NORM;
	
	
	
};
//----------------------------------------------------------------------



/// @brief Make Monte Carlo integration of ISR distribution
/// @param filename filename of root and pdf files to write fistograms
/// @param keyISR   ISR type (a) - 0; (b) - 1; (c) - 2;
/// @param kDim     integration dimension: 2,3 Machine energy spread OFF/ON
/// @param sigE     energy spread if applicable
/// @param NevTot   statistics of integration
void MakeISR( string filename = "histo1", Int_t keyISR = 2, Int_t kDim = 3, Double_t sigE = 0.0041, long NevTot =  1000000 ) 
{
	
	//allocate density and set up
	TDensity * rho= new TDensity(); 
	//set convolution 2/3 - OFF/ON
	rho->m_kDim = kDim;
	//set ISR type
	rho->m_KeyISR = keyISR;
	//set energy spread if applicable (kDim = 3)
	rho->m_sigE = sigE;
	
	
	Double_t MCresult,MCerror;
	//=========================================================
	Int_t  nBin     =       8;   // Number of bins in build-up
	Int_t  OptRej   =       0;   // Wted events for OptRej=0; wt=1 for OptRej=1 (default)
	Int_t  OptDrive =       2;   // (D=2) Option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
	Int_t  EvPerBin =      25;   // Maximum events (equiv.) per bin in buid-up

	#ifdef DEBUG
		Int_t  Chat     =       1;   // Chat level
	#else
		Int_t  Chat     =       0;   // Chat level
	#endif
	//=========================================================
	TRandom *PseRan   = new TRandom3();  // Create random number generator
	TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
	PseRan->SetSeed(4357);
	//=========================================================
	#ifdef DEBUG
		cout<<"*****   Foam version "<<FoamX->GetVersion()<<"    *****"<<endl;
	#else 
		//do nothing
	#endif	

	FoamX->SetkDim(        rho->m_kDim );      // Mandatory!!!
	FoamX->SetnCells(      rho->m_nCells);    // optional
	FoamX->SetnSampl(      rho->m_nSampl);    // optional
	FoamX->SetnBin(        nBin);      // optional
	FoamX->SetOptRej(      OptRej);    // optional
	FoamX->SetOptDrive(    OptDrive);  // optional
	FoamX->SetEvPerBin(    EvPerBin);  // optional
	FoamX->SetChat(        Chat);      // optional
	//===============================
	FoamX->SetRho(rho);
	FoamX->SetPseRan(PseRan);

	
	#ifdef DEBUG
		cout << "Initialze...." << endl;
	#else 
		//do nothing
	#endif	
	
	FoamX->Initialize(); // Initialize simulator
	
	#ifdef DEBUG
		cout << "...initialze DONE" << endl;
	#else 
		//do nothing
	#endif
	
	
	double Xnorm, errel;
	FoamX->GetIntNorm(Xnorm,errel);   // universal normalization
 
	//FoamX->Write("FoamX");     // Writing Foam on the disk, TESTING PERSISTENCY!!!

	long nCalls=FoamX->GetnCalls();


	Double_t *MCvect = new Double_t[kDim]; // vector generated in the MC run
	
	//Histograms
	int    nbx=200;       // linear scale
	double Emax = rho->m_MH +3*rho->m_GamH;
	double Emin = rho->m_MH -3*rho->m_GamH;
	
	TH1D * h_Ene = new TH1D("h_Ene", "h_Ene",nbx, Emin, Emax);
	h_Ene->Sumw2();

	//norm and no event - only two bins are used
	TH1D * h_NORM = new TH1D("h_NORM", "h_NORMA",2,0.0,2.0);
	h_NORM->Sumw2();	
	
	
	//event loop
	for(long loop=0; loop<NevTot; loop++)
	{
		//generate MC event
		FoamX->MakeEvent();
		FoamX->GetMCvect( MCvect);
		
		//get point info
		double E, MCwt;
		//double y;
		
		MCwt=FoamX->GetMCwt();
		E = rho->m_E;
		//y = rho->m_y;
    
		//filling histograms
		h_Ene->Fill( E, MCwt );
		
		//  Fill special normalization histogram hNORM
		h_NORM->Fill(0.5, Xnorm);   // 1-st bin = Normal*Nevtot
		h_NORM->Fill(1.5, 1);         // 2-nd bin = Nevtot
	
	
		#ifdef DEBUG
			//print progress bar
			if( loop % 1000 == 0 ) 
				progressBar( loop, NevTot );
		#else
			//do nothing
		#endif
		
	}

	#ifdef DEBUG
		cout << "====== Events generated, entering Finalize" << endl;
	#else
		//do nothing
	#endif

	Double_t eps = 0.0005;
	Double_t Effic, WtMax, AveWt, Sigma;
	Double_t IntNorm, Errel;
	FoamX->Finalize( IntNorm, Errel );     // final printout
	FoamX->GetIntegMC( MCresult, MCerror );  // get MC intnegral
	FoamX->GetWtParams( eps, AveWt, WtMax, Sigma ); // get MC wt parameters
	Effic=0; 
	
	if(WtMax>0) 
		Effic=AveWt/WtMax;
	
	#ifdef DEBUG		
		cout << "================================================================" << endl;
		cout << " MCresult= " << MCresult << " +- " << MCerror << " RelErr= "<< MCerror/MCresult << endl;
		cout << " Dispersion/<wt>= " << Sigma/AveWt << endl;
		cout << "      <wt>/WtMax= " << Effic <<",    for epsilon = "<<eps << endl;
		cout << " nCalls (initialization only) =   " << nCalls << endl;
		cout << "================================================================" << endl;
	#else
		//do nothing
	#endif


	delete [] MCvect;
	
	//make sigma histograms by normalization 
	TH1D *h_SigEne = (TH1D*)h_Ene->Clone("h_SigEne");
	HistNorm( h_NORM, h_SigEne );
	
	//remove statistics 
	h_SigEne->SetStats(0);
	
	
	//Save pdf
	Float_t  WidPix, HeiPix;
	WidPix = 800; HeiPix =  800;
	TCanvas *cCanv = new TCanvas("cCanv","cCanv", 100,100, WidPix,HeiPix);
	cCanv->SetFillColor(10);
	h_SigEne->Draw();
	cCanv->Update();
	cCanv->Print((filename+string(".pdf")).c_str(), "");
	delete cCanv;
	
	
	//Save ROOT
	TFile RootFile( (filename+string(".root")).c_str(), "RECREATE", "histograms");
	RootFile.ls();
	h_Ene->Write();
	h_SigEne->Write();	
	h_NORM->Write();
	RootFile.Close();

	//cleaning;
	delete FoamX;
	delete PseRan;
	delete rho;
	
	delete h_Ene;
	delete h_SigEne;
	delete h_NORM;
	
	
	
};
//----------------------------------------------------------------------

/// @brief Make Born term as histogram of given shape
/// @param filename file to save Born histogram
/// @param pattern  file that contains pattern histogram
void MakeBornSame( string filename = "BornH", string pattern = "./histo-sig0-isr2.root" )
{
	
	#ifdef DEBUG
		cout<< "====MakeBornSame ===="<<endl;
	#endif
	
	// Get acces to MC generator object
    TDensity * Density = new TDensity();
    
	// determine shape of the histogram
	TFile DiskFileA( pattern.c_str() );
    TH1D *h_Ene = (TH1D*)DiskFileA.Get("h_Ene");
    TH1D *h_SigEne = (TH1D*)h_Ene->Clone("h_SigEne");
    int      nbt   = h_SigEne->GetNbinsX();
    Double_t xmax  = h_SigEne->GetXaxis()->GetXmax();
    Double_t xmin  = h_SigEne->GetXaxis()->GetXmin();
    
    
    //make plot
    double svar,sigmaB,x,dxl=(xmax-xmin)/nbt;
    for(int i=1; i<=nbt;i++)
    {
		x = xmin +(i)*dxl;        // RHS of the bin 
		svar=x*x;
		sigmaB = Density->BornH(svar);
		//    sigmaB= 1.0;
		h_SigEne->SetBinContent(i,sigmaB); // Born sigma(Ene)
		h_SigEne->SetBinError(i,0.0);
		
     };
     
    
    //Save pdf
	Float_t  WidPix, HeiPix;
	WidPix = 800; HeiPix =  800;
	TCanvas *cCanv = new TCanvas("cCanv","cCanv", 100,100, WidPix,HeiPix);
	cCanv->SetFillColor(10);
	h_SigEne->Draw();
	cCanv->Update();
	cCanv->Print((filename+string(".pdf")).c_str(), "");
	delete cCanv;
	
	
	//Save ROOT
	TFile RootFile( (filename+string(".root")).c_str(), "RECREATE", "histograms");
	RootFile.ls();
	h_SigEne->Write();	
	RootFile.Close();
    
     
};
//----------------------------------------------------------------------


/// @brief Make Born term as histogram
/// @param filename file to save Born histogram
void MakeBorn( string filename = "BornH" )
{
	
	#ifdef DEBUG
		cout<< "====MakeBorn ===="<<endl;
	#endif
	
	// Get acces to MC generator object
    TDensity * Density = new TDensity();

      
    //make plot
 	int    nbx=200;       // linear scale
	double Emax = Density->m_MH +3*Density->m_GamH;
	double Emin = Density->m_MH -3*Density->m_GamH;
    double dE   = (Emax-Emin)/(nbx+1);
    
    // define shape of the histogram
	TH1D * h_SigEne = new TH1D("h_SigEne", "h_SigEne",nbx,Emin,Emax);
	h_SigEne->Sumw2();
	
    
	//histogram loop
    double svar,sigmaB;
    double E = Emin;
    for(int i=1; i<=nbx; i++)
    {
		E = Emin +(i)*dE; 
		svar=E*E;
		sigmaB = Density->BornH(svar);
		h_SigEne->SetBinContent(i,sigmaB); // Born sigma(Ene)
		h_SigEne->SetBinError(i,0.0);
     
    };
     
    
    //Save pdf
	Float_t  WidPix, HeiPix;
	WidPix = 800; HeiPix =  800;
	TCanvas *cCanv = new TCanvas("cCanv","cCanv", 100,100, WidPix,HeiPix);
	cCanv->SetFillColor(10);
	h_SigEne->Draw();
	cCanv->Update();
	cCanv->Print((filename+string(".pdf")).c_str(), "");
	delete cCanv;
	
	
	//Save ROOT
	TFile RootFile( (filename+string(".root")).c_str(), "RECREATE", "histograms");
	RootFile.ls();
	h_SigEne->Write();	
	RootFile.Close();
    
     
};
//----------------------------------------------------------------------



///@returns bin error for histogram at X-value = E
///@param E - value at which error is returned
///@param histo - pointer to the histogram
double getHistError( double E, TH1 * histo )
{
	return( histo-> GetBinError(histo->GetXaxis()->FindBin(E)) );
};
//----------------------------------------------------------------------


/// @returns TH1D pointer to the Born convolution histogram at given energy spread
/// @warning returned histogram should be deleted when not needed to avoid memory leak
/// @param sigE   energy spread of convolution at which distribution is calculated
/// @param NevTot   statistics of integration
TH1D * getBornHistoAtSigE( Double_t sigE = 0.0042, long NevTot =  1000000 )
{

	//allocate density and set up
	TDensity * rho= new TDensity(); 
	//set convolution ON
	rho->m_kDim = 2;
	//turn off ISR, pure Born
	rho->m_ISROn = 0;
	//set ISR type
	rho->m_KeyISR = 2;
	//set energy spread if applicable (kDim = 3)
	rho->m_sigE = sigE;
	
	//=========================================================
	Int_t  nBin     =       8;   // Number of bins in build-up
	Int_t  OptRej   =       0;   // Wted events for OptRej=0; wt=1 for OptRej=1 (default)
	Int_t  OptDrive =       2;   // (D=2) Option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
	Int_t  EvPerBin =      25;   // Maximum events (equiv.) per bin in buid-up
	
	#ifdef DEBUG
		Int_t  Chat     =       1;   // Chat level
	#else
		Int_t  Chat     =       0;   // Chat level
	#endif
	
	//=========================================================
	TRandom *PseRan   = new TRandom3();  // Create random number generator
	TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
	PseRan->SetSeed(4357);
	//=========================================================
	FoamX->SetkDim(        rho->m_kDim );      // Mandatory!!!
	FoamX->SetnCells(      rho->m_nCells);    // optional
	FoamX->SetnSampl(      rho->m_nSampl);    // optional
	FoamX->SetnBin(        nBin);      // optional
	FoamX->SetOptRej(      OptRej);    // optional
	FoamX->SetOptDrive(    OptDrive);  // optional
	FoamX->SetEvPerBin(    EvPerBin);  // optional
	FoamX->SetChat(        Chat);      // optional
	//===============================
	FoamX->SetRho(rho);
	FoamX->SetPseRan(PseRan);

	
	FoamX->Initialize(); // Initialize simulator
	
	double Xnorm, errel;
	FoamX->GetIntNorm(Xnorm,errel);   // universal normalization
 
	//FoamX->Write("FoamX");     // Writing Foam on the disk, TESTING PERSISTENCY!!!
	
	//Histograms
	int    nbx=200;       // linear scale
	double Emax = rho->m_MH +3*rho->m_GamH;
	double Emin = rho->m_MH -3*rho->m_GamH;
	
	TH1D * h_Ene = new TH1D("h_Ene", "h_Ene",nbx, Emin, Emax);
	h_Ene->Sumw2();

	//norm and no event - only two bins are used
	TH1D * h_NORM = new TH1D("h_NORM", "h_NORMA",2,0.0,2.0);
	h_NORM->Sumw2();	
	
	
	//event loop
	for(long loop=0; loop<NevTot; loop++)
	{
		//generate MC event
		FoamX->MakeEvent();
		
		//get point info
		double E, MCwt;
		//double y;

		MCwt=FoamX->GetMCwt();
		E = rho->m_E;
		//y = rho->m_y;
    
		//filling histograms
		h_Ene->Fill( E, MCwt );
		
		//  Fill special normalization histogram hNORM
		h_NORM->Fill(0.5, Xnorm);   // 1-st bin = Normal*Nevtot
		h_NORM->Fill(1.5, 1);         // 2-nd bin = Nevtot
	}

	
	//make sigma histograms by normalization 
	TH1D *h_SigEne = (TH1D*)h_Ene->Clone("h_SigEne");
	HistNorm( h_NORM, h_SigEne );	


	//cleaning;
	delete FoamX;
	delete PseRan;
	delete rho;
	
	delete h_Ene;
	delete h_NORM;
		
	
	return h_SigEne;	
	
};
//----------------------------------------------------------------------


/// @returns TH1D pointer to ISR (c) histogram at given energy spread
/// @warning returned histogram should be deleted when not needed to avoid memory leak
/// @param sigE   energy spread of convolution at which distribution is calculated
/// @param NevTot   statistics of integration
TH1D * getISRHistoAtSigE( Double_t sigE = 0.0042, long NevTot =  1000000 )
{

	//allocate density and set up
	TDensity * rho= new TDensity(); 
	//set convolution ON
	rho->m_kDim = 3;
	//set ISR type
	rho->m_KeyISR = 2;
	//set energy spread if applicable (kDim = 3)
	rho->m_sigE = sigE;
	
	//=========================================================
	Int_t  nBin     =       8;   // Number of bins in build-up
	Int_t  OptRej   =       0;   // Wted events for OptRej=0; wt=1 for OptRej=1 (default)
	Int_t  OptDrive =       2;   // (D=2) Option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
	Int_t  EvPerBin =      25;   // Maximum events (equiv.) per bin in buid-up
	
	#ifdef DEBUG
		Int_t  Chat     =       1;   // Chat level
	#else
		Int_t  Chat     =       0;   // Chat level
	#endif
	//=========================================================
	TRandom *PseRan   = new TRandom3();  // Create random number generator
	TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
	PseRan->SetSeed(4357);
	//=========================================================
	FoamX->SetkDim(        rho->m_kDim );      // Mandatory!!!
	FoamX->SetnCells(      rho->m_nCells);    // optional
	FoamX->SetnSampl(      rho->m_nSampl);    // optional
	FoamX->SetnBin(        nBin);      // optional
	FoamX->SetOptRej(      OptRej);    // optional
	FoamX->SetOptDrive(    OptDrive);  // optional
	FoamX->SetEvPerBin(    EvPerBin);  // optional
	FoamX->SetChat(        Chat);      // optional
	//===============================
	FoamX->SetRho(rho);
	FoamX->SetPseRan(PseRan);

	
	FoamX->Initialize(); // Initialize simulator
	
	double Xnorm, errel;
	FoamX->GetIntNorm(Xnorm,errel);   // universal normalization
 
	//FoamX->Write("FoamX");     // Writing Foam on the disk, TESTING PERSISTENCY!!!
	
	//Histograms
	int    nbx=200;       // linear scale
	double Emax = rho->m_MH +3*rho->m_GamH;
	double Emin = rho->m_MH -3*rho->m_GamH;
	
	TH1D * h_Ene = new TH1D("h_Ene", "h_Ene",nbx, Emin, Emax);
	h_Ene->Sumw2();

	//norm and no event - only two bins are used
	TH1D * h_NORM = new TH1D("h_NORM", "h_NORMA",2,0.0,2.0);
	h_NORM->Sumw2();	
	
	
	//event loop
	for(long loop=0; loop<NevTot; loop++)
	{
		//generate MC event
		FoamX->MakeEvent();
		
		//get point info
		double E, MCwt;
		//double y;

		MCwt=FoamX->GetMCwt();
		E = rho->m_E;
		//y = rho->m_y;
    
		//filling histograms
		h_Ene->Fill( E, MCwt );
		
		//  Fill special normalization histogram hNORM
		h_NORM->Fill(0.5, Xnorm);   // 1-st bin = Normal*Nevtot
		h_NORM->Fill(1.5, 1);         // 2-nd bin = Nevtot
	}

	
	//make sigma histograms by normalization 
	TH1D *h_SigEne = (TH1D*)h_Ene->Clone("h_SigEne");
	HistNorm( h_NORM, h_SigEne );		

	//cleaning;
	delete FoamX;
	delete PseRan;
	delete rho;
	
	delete h_Ene;
	delete h_NORM;
		
	
	return h_SigEne;	
	
};
//----------------------------------------------------------------------



/// @returns Born convolution value at given energy E and for given energy spread
/// @param E      energy at which distribution is calculated
/// @param sigE   energy spread of convolution at which distribution is calculated
/// @param NevTot   statistics of integration
double getBornatE( Double_t E = 125.09, Double_t sigE = 0.0042, long NevTot =  1000000 )
{

	//allocate density and set up
	TDensity * rho= new TDensity(); 
	//set convolution ON
	rho->m_kDim = 2;
	//turn off ISR, pure Born
	rho->m_ISROn = 0;
	//set ISR type
	rho->m_KeyISR = 2;
	//set energy spread if applicable (kDim = 3)
	rho->m_sigE = sigE;
	
	
	//=========================================================
	Int_t  nBin     =       8;   // Number of bins in build-up
	Int_t  OptRej   =       0;   // Wted events for OptRej=0; wt=1 for OptRej=1 (default)
	Int_t  OptDrive =       2;   // (D=2) Option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
	Int_t  EvPerBin =      25;   // Maximum events (equiv.) per bin in buid-up
	
	#ifdef DEBUG
		Int_t  Chat     =       1;   // Chat level
	#else
		Int_t  Chat     =       0;   // Chat level
	#endif

	//=========================================================
	TRandom *PseRan   = new TRandom3();  // Create random number generator
	TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
	PseRan->SetSeed(4357);
	//=========================================================
	FoamX->SetkDim(        rho->m_kDim );      // Mandatory!!!
	FoamX->SetnCells(      rho->m_nCells);    // optional
	FoamX->SetnSampl(      rho->m_nSampl);    // optional
	FoamX->SetnBin(        nBin);      // optional
	FoamX->SetOptRej(      OptRej);    // optional
	FoamX->SetOptDrive(    OptDrive);  // optional
	FoamX->SetEvPerBin(    EvPerBin);  // optional
	FoamX->SetChat(        Chat);      // optional
	//===============================
	FoamX->SetRho(rho);
	FoamX->SetPseRan(PseRan);

	
	FoamX->Initialize(); // Initialize simulator
	
	double Xnorm, errel;
	FoamX->GetIntNorm(Xnorm,errel);   // universal normalization
 
	//FoamX->Write("FoamX");     // Writing Foam on the disk, TESTING PERSISTENCY!!!
	
	//Histograms
	int    nbx=200;       // linear scale
	double Emax = rho->m_MH +3*rho->m_GamH;
	double Emin = rho->m_MH -3*rho->m_GamH;
	
	TH1D * h_Ene = new TH1D("h_Ene", "h_Ene",nbx, Emin, Emax);
	h_Ene->Sumw2();

	//norm and no event - only two bins are used
	TH1D * h_NORM = new TH1D("h_NORM", "h_NORMA",2,0.0,2.0);
	h_NORM->Sumw2();	
	
	
	//event loop
	for(long loop=0; loop<NevTot; loop++)
	{
		//generate MC event
		FoamX->MakeEvent();
		
		//get point info
		double E, MCwt;
		//double y;

		MCwt=FoamX->GetMCwt();
		E = rho->m_E;
		//y = rho->m_y;
    
		//filling histograms
		h_Ene->Fill( E, MCwt );
		
		//  Fill special normalization histogram hNORM
		h_NORM->Fill(0.5, Xnorm);   // 1-st bin = Normal*Nevtot
		h_NORM->Fill(1.5, 1);         // 2-nd bin = Nevtot
	}

	
	//make sigma histograms by normalization 
	TH1D *h_SigEne = (TH1D*)h_Ene->Clone("h_SigEne");
	HistNorm( h_NORM, h_SigEne );	
	
	
	//get value
	double value = h_SigEne->Interpolate(E);
	

	//cleaning;
	delete FoamX;
	delete PseRan;
	delete rho;
	
	delete h_Ene;
	delete h_SigEne;
	delete h_NORM;
		
	
	return value;	
	
};
//----------------------------------------------------------------------



/// @returns ISR (c) value at given energy E and for given energy spread
/// @param E      energy at which distribution is calculated
/// @param sigE   energy spread of convolution at which distribution is calculated
/// @param NevTot   statistics of integration
double getISRatE( Double_t E = 125.09, Double_t sigE = 0.0042, long NevTot =  1000000 )
{

	//allocate density and set up
	TDensity * rho= new TDensity(); 
	//set convolution ON
	rho->m_kDim = 3;
	//set ISR type
	rho->m_KeyISR = 2;
	//set energy spread if applicable (kDim = 3)
	rho->m_sigE = sigE;
	
	//=========================================================
	Int_t  nBin     =       8;   // Number of bins in build-up
	Int_t  OptRej   =       0;   // Wted events for OptRej=0; wt=1 for OptRej=1 (default)
	Int_t  OptDrive =       2;   // (D=2) Option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
	Int_t  EvPerBin =      25;   // Maximum events (equiv.) per bin in buid-up
	
	#ifdef DEBUG
		Int_t  Chat     =       1;   // Chat level
	#else
		Int_t  Chat     =       0;   // Chat level
	#endif

	//=========================================================
	TRandom *PseRan   = new TRandom3();  // Create random number generator
	TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
	PseRan->SetSeed(4357);
	//=========================================================
	FoamX->SetkDim(        rho->m_kDim );      // Mandatory!!!
	FoamX->SetnCells(      rho->m_nCells);    // optional
	FoamX->SetnSampl(      rho->m_nSampl);    // optional
	FoamX->SetnBin(        nBin);      // optional
	FoamX->SetOptRej(      OptRej);    // optional
	FoamX->SetOptDrive(    OptDrive);  // optional
	FoamX->SetEvPerBin(    EvPerBin);  // optional
	FoamX->SetChat(        Chat);      // optional
	//===============================
	FoamX->SetRho(rho);
	FoamX->SetPseRan(PseRan);

	
	FoamX->Initialize(); // Initialize simulator
	
	double Xnorm, errel;
	FoamX->GetIntNorm(Xnorm,errel);   // universal normalization
 
	//FoamX->Write("FoamX");     // Writing Foam on the disk, TESTING PERSISTENCY!!!
	
	//Histograms
	int    nbx=200;       // linear scale
	double Emax = rho->m_MH +3*rho->m_GamH;
	double Emin = rho->m_MH -3*rho->m_GamH;
	
	TH1D * h_Ene = new TH1D("h_Ene", "h_Ene",nbx, Emin, Emax);
	h_Ene->Sumw2();

	//norm and no event - only two bins are used
	TH1D * h_NORM = new TH1D("h_NORM", "h_NORMA",2,0.0,2.0);
	h_NORM->Sumw2();	
	
	
	//event loop
	for(long loop=0; loop<NevTot; loop++)
	{
		//generate MC event
		FoamX->MakeEvent();
		
		//get point info
		double E, MCwt;
		//double y;

		MCwt=FoamX->GetMCwt();
		E = rho->m_E;
		//y = rho->m_y;
    
		//filling histograms
		h_Ene->Fill( E, MCwt );
		
		//  Fill special normalization histogram hNORM
		h_NORM->Fill(0.5, Xnorm);   // 1-st bin = Normal*Nevtot
		h_NORM->Fill(1.5, 1);         // 2-nd bin = Nevtot
	}

	
	//make sigma histograms by normalization 
	TH1D *h_SigEne = (TH1D*)h_Ene->Clone("h_SigEne");
	HistNorm( h_NORM, h_SigEne );	
	
	
	//get value
	double value = h_SigEne->Interpolate(E);
	

	//cleaning;
	delete FoamX;
	delete PseRan;
	delete rho;
	
	delete h_Ene;
	delete h_SigEne;
	delete h_NORM;
		
	
	return value;	
	
};
//----------------------------------------------------------------------



////////////////////////////////////////////////////////////////////////
//Plotting functions
////////////////////////////////////////////////////////////////////////

/// @brief Read root histograms from disk and plot ISR ratios for ISR (a), (b) and (c)
void plotISRabc( void )
{
	
	//asociate root files - files must exist !!!
	TFile DiskFileA(  "./histo-sig0-isr2.root");
	TFile DiskFileB(  "./histo-sig04-isr2.root"); 
	TFile DiskFileC(  "./histo-sig08-isr2.root"); 
	TFile DiskFileISR0(  "./histo-sig0-isr0.root"); 
	TFile DiskFileISR1(  "./histo-sig0-isr1.root"); 
	TFile DiskFileISR2(  "./histo-sig0-isr2.root"); 
	TFile DiskFileBornH(  "./BornH.root"); 
	
	//load histograms
	//TH1D *h_Ene    = (TH1D*)DiskFileISR0.Get("h_Ene");
	TH1D *h_SigEne   = (TH1D*)DiskFileISR0.Get("h_SigEne");
	//TH1D *h_EneB   = (TH1D*)DiskFileISR1.Get("h_Ene");
	TH1D *h_SigEneB  = (TH1D*)DiskFileISR1.Get("h_SigEne");
	//TH1D *h_EneC   = (TH1D*)DiskFileISR2.Get("h_Ene");
	TH1D *h_SigEneC  = (TH1D*)DiskFileISR2.Get("h_SigEne");
	TH1D *h_Born     = (TH1D*)DiskFileBornH.Get("h_SigEne");

	Double_t Emax  = h_SigEne->GetXaxis()->GetXmax();
	Double_t Emin  = h_SigEne->GetXaxis()->GetXmin();
	
	//collect data from Density object
	TDensity* Density = new TDensity();
	double MH   = Density->m_MH;
	double GamH = Density->m_GamH;
  
	cout << "Creating histograms ..." << endl;
  
	// Caption in latex
	TLatex *CaptTa = new TLatex();
	// CaptTa->SetNDC();
	CaptTa->SetTextAlign(11);
	CaptTa->SetTextSize(0.040);

	TLine *line = new TLine(126,0,126,2);
	Float_t  WidPix, HeiPix;
	WidPix = 800; HeiPix =  1000;
	TCanvas *cCanV3b = new TCanvas("cCanVISR","cCanVISR", 200, 50, WidPix,HeiPix);
	cCanV3b->SetFillColor(10);
	cCanV3b->Draw();
	cCanV3b->Divide(0, 2);
	
	//set y logscale
	//cCanV3b->SetLogy();

	
	//make histogram
	TH1D *Hst;
	cCanV3b->cd(1);
	Hst = h_Born;
 
	Hst->SetStats(0);
	Hst->SetTitle(0);


	Hst->SetStats(0);
	Hst->SetTitle(0);
	Hst->GetYaxis()->CenterTitle();
	Hst->GetYaxis()->SetLabelSize(0.06);
	Hst->GetYaxis()->SetTitleSize(0.085);
	
	#if defined( ELECTRON )
		Hst->GetYaxis()->SetTitle( "#sigma(s) [fb]" );
	#elif defined( MUON )
		Hst->GetYaxis()->SetTitle("#sigma(s) [pb]");
	#endif 	
	
	
	Hst->GetXaxis()->CenterTitle();
	Hst->GetXaxis()->SetTitleOffset(0.90);
	Hst->GetYaxis()->SetTitleOffset(0.60);
	Hst->GetXaxis()->SetTitleSize(0.070);
	Hst->GetXaxis()->SetLabelSize(0.07);
	//Hst->GetXaxis()->SetTitle( "#sqrt{s}   [GeV]" );
	Hst->GetXaxis()->SetTitle( "" );

	Hst->SetLineColor(kBlack);
	Hst->DrawCopy("l");

	h_SigEne->SetLineColor(kBlue);
	h_SigEne->DrawCopy("lsame");

	h_SigEneB->SetLineColor(kRed);
	h_SigEneB->DrawCopy("lsame");

	h_SigEneC->SetLineColor(kMagenta);
	h_SigEneC->DrawCopy("lsame");

	//place captions on the plot
	CaptTa->DrawLatex(MH, Hst->Interpolate(MH) ,    "  Born");
	CaptTa->SetTextColor(kBlue);
	CaptTa->DrawLatex(MH, h_SigEne->Interpolate(MH)-0.1 ,  "   (a)");
	CaptTa->SetTextColor(kRed);
	CaptTa->DrawLatex(MH+0.001, h_SigEneB->Interpolate(MH) , "   (b)");
	CaptTa->SetTextColor(kMagenta);
	CaptTa->DrawLatex(MH, h_SigEneC->Interpolate(MH) , "   (c)");

	//Make ratio histograms
	cCanV3b->cd(2);
	TH1D *h_SigRat = (TH1D*)h_SigEne->Clone("h_SigRat");
	h_SigRat->Divide( h_Born );

	TH1D *h_SigRatB = (TH1D*)h_SigEneB->Clone("h_SigRatB");
	h_SigRatB->Divide( h_Born );

	TH1D *h_SigRatC = (TH1D*)h_SigEneC->Clone("h_SigRatC");
	h_SigRatC->Divide( h_Born );


	h_SigRat->SetStats(0);
	h_SigRat->SetTitle(0);
	h_SigRat->GetYaxis()->CenterTitle();
	h_SigRat->GetYaxis()->SetLabelSize(0.06);
	h_SigRat->GetYaxis()->SetTitleSize(0.085);
	h_SigRat->GetYaxis()->SetTitle( "#sigma/#sigma_{Born}" );
	h_SigRat->GetXaxis()->CenterTitle();
	h_SigRat->GetXaxis()->SetTitleOffset(0.50);
	h_SigRat->GetYaxis()->SetTitleOffset(0.60);
	h_SigRat->GetXaxis()->SetTitleSize(0.085);
	h_SigRat->GetXaxis()->SetLabelSize(0.07);
	h_SigRat->GetXaxis()->SetTitle( "#sqrt{s}   [GeV]" );


	h_SigRat->GetXaxis()->SetLabelOffset(999);
	h_SigRat->GetXaxis()->SetLabelSize(0);
	


	h_SigRat->SetMinimum(0);
	h_SigRat->SetMaximum(2);
	h_SigRat->SetLineColor(kBlue);
	h_SigRat->DrawCopy("l");

	h_SigRatB->SetLineColor(kRed);
	h_SigRatB->DrawCopy("lsame");

	h_SigRatC->SetLineColor(kMagenta);
	h_SigRatC->DrawCopy("lsame");


	TH1D* h_Hone = new TH1D("h_Hone","line",1, Emin,Emax);
	h_Hone->SetBinContent(1,1);
	h_Hone->DrawCopy("lsame");

	CaptTa->SetTextAlign(11);
	CaptTa->SetTextColor(kBlue);
	CaptTa->DrawLatex(MH-2.3*GamH, 1.10*h_SigRat->Interpolate(MH-2.3*GamH)-0.12  ,  "(a)");
  
	CaptTa->SetTextColor(kRed);
	CaptTa->DrawLatex(MH-2.3*GamH+0.001, 1.10*h_SigRatB->Interpolate(MH-2.3*GamH) ,  "(b)");
	
	CaptTa->SetTextColor(kMagenta);
	CaptTa->DrawLatex(MH-2.3*GamH, 1.05*h_SigRatC->Interpolate(MH-2.3*GamH) ,  "(c)");

	TLine *line2 = new TLine(MH,0,MH,2);
	line2->Draw();

	cCanV3b->Update();
	//cCanV3b->cd();
	
	cCanV3b->Print("ISRabc.eps", "");
  
	cout<< "...creating histograms DONE"<<endl;
  
  
	// print suppresion table
	cout << "Suppresion factors..." << endl;
	
	//select histograms
	TH1D *h_Enea  = h_SigEne;
	TH1D *h_Eneb  = h_SigEneB;
	TH1D *h_Enec  = h_SigEneC;

	const double E0   = MH;
	const double E0pG = MH + GamH;
	const double E0mG = MH - GamH;
	
  	cout << "E0 = " << E0  << endl;
	cout << "Gamma = " << GamH << endl;
	cout << "E0+Gamma = " << E0pG  << endl;
	cout << "E0-Gamma = " << E0mG  << endl;
	
	const double sig0MH   = Density->BornH(E0*E0);
	const double sig0MHpG = Density->BornH(E0pG*E0pG);
	const double sig0MHmG = Density->BornH(E0mG*E0mG);
	
	cout << setw(10) <<" type" << "\t" << setw(10) << " Born " << "\t" << setw(30) << "(a)" << "\t" << setw(30) << "(b)" << "\t" << setw(30) << "(c)" << "\t" << setw(10) << "(a)/B" << "\t" << setw(10) <<"(b)/B" << "\t" << setw(10) << "(c)/B" << endl;
	
	cout << setw(10) << "E0" << "\t" << setw(10) << sig0MH  << "\t" << setw(14) << h_Enea->Interpolate(E0) << "+-" << setw(14) << getHistError( E0, h_Enea ) << "\t" << setw(14) << h_Eneb->Interpolate(E0) << "+-" << setw(14) << getHistError( E0, h_Eneb ) << "\t" << setw(14) << h_Enec->Interpolate(E0) << "+-" << setw(14) << getHistError( E0, h_Enec ) << "\t" << setw(10) << h_Enea->Interpolate(E0)/sig0MH << "\t" << setw(10) << h_Eneb->Interpolate(E0)/sig0MH << "\t" << setw(10) << h_Enec->Interpolate(E0)/sig0MH << endl;
	
	cout << setw(10) << "E0+Gamma" << "\t" << setw(9) << sig0MHpG << "\t" << setw(14) << h_Enea->Interpolate(E0pG) << "+-" << setw(14) << getHistError( E0pG, h_Enea ) << "\t" << setw(14) << h_Eneb->Interpolate(E0pG) << "+-" << setw(14) << getHistError( E0pG, h_Eneb ) << "\t" << setw(14) << h_Enec->Interpolate(E0pG) << "+-" << setw(14) << getHistError( E0pG, h_Enec ) << "\t" << setw(10) << h_Enea->Interpolate(E0pG)/sig0MHpG << "\t" << setw(10) << h_Eneb->Interpolate(E0pG)/sig0MHpG << "\t" << setw(10) << h_Enec->Interpolate(E0pG)/sig0MHpG  << endl;
	
	cout << setw(10) << "E0-Gamma" << "\t" << setw(10) << sig0MHmG  << "\t" << setw(14) << h_Enea->Interpolate(E0mG) << "+-" << setw(14) << getHistError( E0mG, h_Enea ) << "\t" << setw(14) << h_Eneb->Interpolate(E0mG) << "+-" << setw(14) << getHistError( E0mG, h_Eneb ) << "\t" << setw(14) << h_Enec->Interpolate(E0mG) << "+-" << setw(14) << getHistError( E0mG, h_Enec ) << "\t" << setw(10) << h_Enea->Interpolate(E0mG)/sig0MHmG << "\t" << setw(10) << h_Eneb->Interpolate(E0mG)/sig0MHmG << "\t"<< setw(10) << h_Enec->Interpolate(E0mG)/sig0MHmG << endl;


  
	cout << "...suppresion factors DONE" << endl;
	
	
	//cleaning
	delete Density;
	delete h_Hone;
	delete cCanV3b;
	delete line;
	delete line2;
	delete CaptTa;

	
};
//----------------------------------------------------------------------


/// @brief Make plots of ISR (c) for three values of spread
void plotISR123()
{

	//asociate root files - files must exist !!!
	TFile DiskFileA(  "./histo-sig0-isr2.root");
	TFile DiskFileB(  "./histo-sig04-isr2.root"); 
	TFile DiskFileC(  "./histo-sig08-isr2.root"); 
	TFile DiskFileISR0(  "./histo-sig0-isr0.root"); 
	TFile DiskFileISR1(  "./histo-sig0-isr1.root"); 
	TFile DiskFileISR2(  "./histo-sig0-isr2.root"); 
	TFile DiskFileBornH(  "./BornH.root");   

	//load histograms
	//TH1D *h_Ene    = (TH1D*)DiskFileA.Get("h_Ene");
	TH1D *h_SigEne   = (TH1D*)DiskFileA.Get("h_SigEne");
	//TH1D *h_EneB   = (TH1D*)DiskFileB.Get("h_Ene");
	TH1D *h_SigEneB  = (TH1D*)DiskFileB.Get("h_SigEne");
	//TH1D *h_EneC   = (TH1D*)DiskFileC.Get("h_Ene");
	TH1D *h_SigEneC  = (TH1D*)DiskFileC.Get("h_SigEne");
	TH1D *h_Born     = (TH1D*)DiskFileBornH.Get("h_SigEne");

	Double_t Emax  = h_SigEne->GetXaxis()->GetXmax();
	Double_t Emin  = h_SigEne->GetXaxis()->GetXmin();

  	//collect data from Density object
	TDensity* Density = new TDensity();
	double MH  = Density->m_MH;
	double GamH= Density->m_GamH;
  
	
	cout << "Creating histograms ..." << endl;
  
	
	// Caption in latex
	TLatex *CaptTa = new TLatex();
	CaptTa->SetTextAlign(11);
	CaptTa->SetTextSize(0.040);

	TLine *line = new TLine(126,0,126,2);
	Float_t  WidPix, HeiPix;
	WidPix = 800; HeiPix =  1000;
	TCanvas *cCanV3b = new TCanvas("cCanV3b","cCanV3b", 200, 50, WidPix,HeiPix);
	cCanV3b->SetFillColor(10);
	cCanV3b->Draw();
	cCanV3b->Divide(0, 2);
  
  
	TH1D *Hst;
  
	cCanV3b->cd(1);
	Hst = h_Born;
	Hst->SetStats(0);
	Hst->SetTitle(0);

	Hst->SetStats(0);
	Hst->SetTitle(0);
	Hst->GetYaxis()->CenterTitle();
	Hst->GetYaxis()->SetLabelSize(0.06);
	Hst->GetYaxis()->SetTitleSize(0.085);
	
	#if defined( ELECTRON )
		Hst->GetYaxis()->SetTitle( "#sigma(s) [fb]" );
	#elif defined( MUON )
		Hst->GetYaxis()->SetTitle("#sigma(s) [pb]");
	#endif 
	
	//Hst->GetYaxis()->SetTitle( "#sigma(s) [fb]" );
	Hst->GetXaxis()->CenterTitle();
	Hst->GetXaxis()->SetTitleOffset(0.90);
	Hst->GetYaxis()->SetTitleOffset(0.60);
	Hst->GetXaxis()->SetTitleSize(0.070);
	Hst->GetXaxis()->SetLabelSize(0.07);
	//Hst->GetXaxis()->SetTitle( "#sqrt{s}   [GeV]" );
	Hst->GetXaxis()->SetTitle( "" );

	Hst->SetLineColor(kBlack);
	//Hst->SetLineWidth(2);
	Hst->DrawCopy("l");

	h_SigEne->SetLineColor(kBlue);
	h_SigEne->DrawCopy("lsame");

	h_SigEneB->SetLineColor(kRed);
	h_SigEneB->DrawCopy("lsame");

	h_SigEneC->SetLineColor(kMagenta);
	h_SigEneC->DrawCopy("lsame");

	CaptTa->DrawLatex(MH, Hst->Interpolate(MH) ,    "  Born");
	CaptTa->SetTextColor(kBlue);
	
	CaptTa->DrawLatex(MH, h_SigEne->Interpolate(MH) ,  "  (1)");
	CaptTa->SetTextColor(kRed);
	
	CaptTa->DrawLatex(MH, h_SigEneB->Interpolate(MH) , "   (2)");
	CaptTa->SetTextColor(kMagenta);
	
	CaptTa->DrawLatex(MH, h_SigEneC->Interpolate(MH) , "   (3)");

	//make ratio plot
	cCanV3b->cd(2);
	TH1D *h_SigRat = (TH1D*)h_SigEne->Clone("h_SigRat");
	h_SigRat->Divide( h_Born );

	TH1D *h_SigRatB = (TH1D*)h_SigEneB->Clone("h_SigRatB");
	h_SigRatB->Divide( h_Born );

	TH1D *h_SigRatC = (TH1D*)h_SigEneC->Clone("h_SigRatC");
	h_SigRatC->Divide( h_Born );

	h_SigRat->SetStats(0);
	h_SigRat->SetTitle(0);
	h_SigRat->GetYaxis()->CenterTitle();
	h_SigRat->GetYaxis()->SetLabelSize(0.06);
	h_SigRat->GetYaxis()->SetTitleSize(0.085);
	h_SigRat->GetYaxis()->SetTitle( "#sigma/#sigma_{Born}" );
	h_SigRat->GetXaxis()->CenterTitle();
	h_SigRat->GetXaxis()->SetTitleOffset(0.50);
	h_SigRat->GetYaxis()->SetTitleOffset(0.60);
	h_SigRat->GetXaxis()->SetTitleSize(0.085);
	h_SigRat->GetXaxis()->SetLabelSize(0.07);
	h_SigRat->GetXaxis()->SetTitle( "#sqrt{s}   [GeV]" );


	h_SigRat->GetXaxis()->SetLabelOffset(999);
	h_SigRat->GetXaxis()->SetLabelSize(0);
	

	h_SigRat->SetMinimum(0);
	h_SigRat->SetMaximum(2);
	h_SigRat->SetLineColor(kBlue);
	h_SigRat->DrawCopy("l");

	h_SigRatB->SetLineColor(kRed);
	h_SigRatB->DrawCopy("lsame");

	h_SigRatC->SetLineColor(kMagenta);
	h_SigRatC->DrawCopy("lsame");

	TH1D *h_Hone = new TH1D("h_Hone","line",1, Emin,Emax);
	h_Hone->SetBinContent(1,1);
	h_Hone->DrawCopy("lsame");

	CaptTa->SetTextAlign(11);
	CaptTa->SetTextColor(kBlue);
	CaptTa->DrawLatex(MH-2.3*GamH, 1.10*h_SigRat->Interpolate(MH-2.3*GamH)  ,  "(1)");
	CaptTa->SetTextColor(kRed);
	CaptTa->DrawLatex(MH-2.3*GamH, 1.10*h_SigRatB->Interpolate(MH-2.3*GamH) ,  "(2)");
	CaptTa->SetTextColor(kMagenta);
	CaptTa->DrawLatex(MH-2.3*GamH, 1.05*h_SigRatC->Interpolate(MH-2.3*GamH) ,  "(3)");

	TLine *line2 = new TLine(MH,0,MH,2);
	line2->Draw();

	cCanV3b->Update();
	cCanV3b->cd();
	cCanV3b->Print("ISR123.eps", "");
  
	cout<< "...creating histograms DONE"<<endl;
	
	// print suppresion table
	cout << "Suppresion factors..." << endl;
	
	//select histograms
	TH1D *h_Enec  = (TH1D*)DiskFileC.Get("h_SigEne");
	TH1D *h_Eneb  = (TH1D*)DiskFileB.Get("h_SigEne");
	TH1D *h_Enea  = (TH1D*)DiskFileA.Get("h_SigEne");

	const double E0 = MH;
	const double E0pG = MH + GamH;
	const double E0mG = MH - GamH;
	
  	cout << "E0 = " << E0  << endl;
	cout << "Gamma = " << GamH << endl;
	cout << "E0+Gamma = " << E0pG  << endl;
	cout << "E0-Gamma = " << E0mG  << endl;
	
	const double sig0MH   = Density->BornH(E0*E0);
	const double sig0MHpG = Density->BornH(E0pG*E0pG);
	const double sig0MHmG = Density->BornH(E0mG*E0mG);
  
  
	// print suppresion table
	cout << setw(10) <<" type" << "\t" << setw(10) << " Born " << "\t" << setw(30) << "(1)" << "\t" << setw(30) << "(2)" << "\t" << setw(30) << "(3)" << "\t" << setw(10) << "(1)/B" << "\t" << setw(10) <<"(2)/B" << "\t" << setw(10) << "(3)/B" << endl;
	
	cout << setw(10) << "E0" << "\t" << setw(10) << sig0MH  << "\t" << setw(14) << h_Enea->Interpolate(E0) << "+-" << setw(14) << getHistError( E0, h_Enea ) << "\t" << setw(14) << h_Eneb->Interpolate(E0) << "+-" << setw(14) << getHistError( E0, h_Eneb ) << "\t" << setw(14) << h_Enec->Interpolate(E0) << "+-" << setw(14) << getHistError( E0, h_Enec ) << "\t" << setw(10) << h_Enea->Interpolate(E0)/sig0MH << "\t" << setw(10) << h_Eneb->Interpolate(E0)/sig0MH << "\t" << setw(10) << h_Enec->Interpolate(E0)/sig0MH << endl;
	
	cout << setw(10) << "E0+Gamma" << "\t" << setw(9) << sig0MHpG << "\t" << setw(14) << h_Enea->Interpolate(E0pG) << "+-" << setw(14) << getHistError( E0pG, h_Enea ) << "\t" << setw(14) << h_Eneb->Interpolate(E0pG) << "+-" << setw(14) << getHistError( E0pG, h_Eneb ) << "\t" << setw(14) << h_Enec->Interpolate(E0pG) << "+-" << setw(14) << getHistError( E0pG, h_Enec ) << "\t" << setw(10) << h_Enea->Interpolate(E0pG)/sig0MHpG << "\t" << setw(10) << h_Eneb->Interpolate(E0pG)/sig0MHpG << "\t" << setw(10) << h_Enec->Interpolate(E0pG)/sig0MHpG  << endl;
	
	cout << setw(10) << "E0-Gamma" << "\t" << setw(10) << sig0MHmG  << "\t" << setw(14) << h_Enea->Interpolate(E0mG) << "+-" << setw(14) << getHistError( E0mG, h_Enea ) << "\t" << setw(14) << h_Eneb->Interpolate(E0mG) << "+-" << setw(14) << getHistError( E0mG, h_Eneb ) << "\t" << setw(14) << h_Enec->Interpolate(E0mG) << "+-" << setw(14) << getHistError( E0mG, h_Enec ) << "\t" << setw(10) << h_Enea->Interpolate(E0mG)/sig0MHmG << "\t" << setw(10) << h_Eneb->Interpolate(E0mG)/sig0MHmG << "\t"<< setw(10) << h_Enec->Interpolate(E0mG)/sig0MHmG << endl;
	
	
	cout << "...suppresion factors DONE" << endl;
	
	
	//cleaning
	delete Density;
	delete h_Hone;
	delete cCanV3b;
	delete line;
	delete line2;
	delete CaptTa;
  
  
 
};
//----------------------------------------------------------------------


/// @brief Plot ISR (c) for different values of energy spread
void plotISRdelta( void )
{
	//asociate root files - files must exist !!!
	TFile DiskFileISR0(  "./histo-sig0-isr2.root");   
	TFile DiskFileISR04(  "./histo-sig04-isr2.root"); 
	TFile DiskFileISR08(  "./histo-sig08-isr2.root"); 
	TFile DiskFileISR15(  "./histo-sig15-isr2.root"); 
	TFile DiskFileISR30(  "./histo-sig30-isr2.root"); 
	TFile DiskFileISR100( "./histo-sig100-isr2.root"); 

	//load histograms	
	TH1D *h_SigEne0    = (TH1D*)DiskFileISR0.Get("h_SigEne");
	TH1D *h_SigEne04   = (TH1D*)DiskFileISR04.Get("h_SigEne");
	TH1D *h_SigEne08   = (TH1D*)DiskFileISR08.Get("h_SigEne");
	TH1D *h_SigEne15   = (TH1D*)DiskFileISR15.Get("h_SigEne");
	TH1D *h_SigEne30   = (TH1D*)DiskFileISR30.Get("h_SigEne");
	TH1D *h_SigEne100  = (TH1D*)DiskFileISR100.Get("h_SigEne");
	
	//collect data from Density object
	TDensity* Density = new TDensity();
	double MH  = Density->m_MH;
	//double GamH= Density->m_GamH;
	
	cout << "Creating histograms ..." << endl;
	
	//set up histograms
	TH1D * h_Hist = h_SigEne0;
	h_Hist->SetStats(0);
	h_Hist->SetStats(0);
	h_Hist->SetTitle(0);
	h_Hist->GetYaxis()->CenterTitle();
	h_Hist->GetYaxis()->SetTitleSize(0.07);
	#if defined( ELECTRON )
		h_Hist->GetYaxis()->SetTitle( "#sigma(s) [fb]" );
	#elif defined( MUON )
		h_Hist->GetYaxis()->SetTitle("#sigma(s) [pb]");
	#endif 
	//h_Hist->GetYaxis()->SetTitle( "#sigma(s) [fb]" );
	h_Hist->GetXaxis()->CenterTitle();
	h_Hist->GetXaxis()->SetTitleSize(0.07);
	h_Hist->GetXaxis()->SetTitle( "#sqrt{s}   [GeV]" );
	h_Hist->GetXaxis()->SetLabelSize(0.055);
	h_Hist->GetYaxis()->SetLabelSize(0.055);
	h_Hist->GetXaxis()->SetTitleOffset(0.72);
	h_Hist->GetYaxis()->SetTitleOffset(0.72);
	
	
	h_SigEne0->SetLineStyle(1);
	h_SigEne0->SetLineColor(1);
	h_SigEne04->SetLineStyle(2);
	h_SigEne04->SetLineColor(2);
	h_SigEne08->SetLineStyle(3);
	h_SigEne08->SetLineColor(3);
	h_SigEne15->SetLineStyle(4);
	h_SigEne15->SetLineColor(4);
	h_SigEne30->SetLineStyle(6);
	h_SigEne30->SetLineColor(6);
	h_SigEne100->SetLineStyle(7);
	h_SigEne100->SetLineColor(7);
	
	
	//legend
	TLegend * leg = new TLegend(0.6,0.6,0.9,0.9);
	
	leg->AddEntry(h_SigEne0,(string("(0) #delta = 0") ).c_str(),"l");
	leg->AddEntry(h_SigEne04,(string("(1) #delta = ")+STR(4.2) ).c_str(),"l");
	leg->AddEntry(h_SigEne08,(string("(2) #delta = ")+STR(8)   ).c_str(),"l");
	leg->AddEntry(h_SigEne15,(string("(3) #delta = ")+STR(15)  ).c_str(),"l");  
	leg->AddEntry(h_SigEne30,(string("(4) #delta = ")+STR(30)  ).c_str(),"l");
	leg->AddEntry(h_SigEne100,(string("(5) #delta = ")+STR(100)).c_str(),"l");
	
	leg->SetHeader(" Energy spread [MeV]");
	//leg->SetTextFont( 62 );
	//leg->SetTextSize(12);

	// Caption in latex
	TLatex *CaptTa = new TLatex();
	CaptTa->SetTextAlign(11);
	CaptTa->SetTextSize(0.030);

	TCanvas* canv1 = new TCanvas("canv","plot");
	
	//set logscale y
	//gPad->SetLogy();
	
	canv1->cd(1);
	
	h_SigEne0->Draw("l");
	h_SigEne04->Draw("lsame");
	h_SigEne08->Draw("lsame");
	h_SigEne15->Draw("lsame");
	h_SigEne30->Draw("lsame");
	h_SigEne100->Draw("lsame");
	leg->Draw("same");
	
	CaptTa->SetTextAlign(11);
	CaptTa->SetTextColor( 1 );
	CaptTa->DrawLatex(MH, 1.02*h_SigEne0->Interpolate(MH)  ,   "(0)");
	CaptTa->SetTextColor( 2 );
	CaptTa->DrawLatex(MH, 1.10*h_SigEne04->Interpolate(MH) ,  "(1)");
	CaptTa->SetTextColor( 3 );
	CaptTa->DrawLatex(MH, 1.10*h_SigEne08->Interpolate(MH) ,  "(2)");
	CaptTa->SetTextColor( 4 );
	CaptTa->DrawLatex(MH, 1.10*h_SigEne15->Interpolate(MH) ,  "(3)");
	CaptTa->SetTextColor( 6 );
	CaptTa->DrawLatex(MH, 1.10*h_SigEne30->Interpolate(MH) ,  "(4)");
	CaptTa->SetTextColor( 7 );
	CaptTa->DrawLatex(MH, 1.10*h_SigEne100->Interpolate(MH) , "(5)");
	
	canv1->Update();
	
	canv1->Print( "ISRspread.eps", "" );		
  
    cout << "...creating histograms DONE" << endl;
	
	
	//cleaning
	delete canv1;
	delete leg;
	delete CaptTa;
	delete Density;
	
	
};
//----------------------------------------------------------------------


/// @brief Plot Born convolution for different values of energy spread
void plotBorndelta( void )
{
	//asociate root files - files must exist !!!
	TFile DiskFileBornH(  "./BornH.root");   
	TFile DiskFileBorn04(  "./histo-sig04-born.root"); 
	TFile DiskFileBorn08(  "./histo-sig08-born.root"); 
	TFile DiskFileBorn15(  "./histo-sig15-born.root"); 
	TFile DiskFileBorn30(  "./histo-sig30-born.root"); 
	TFile DiskFileBorn100( "./histo-sig100-born.root"); 

	//load histograms	
	TH1D *h_SigBorn    = (TH1D*)DiskFileBornH.Get("h_SigEne");
	TH1D *h_SigEne04   = (TH1D*)DiskFileBorn04.Get("h_SigEne");
	TH1D *h_SigEne08   = (TH1D*)DiskFileBorn08.Get("h_SigEne");
	TH1D *h_SigEne15   = (TH1D*)DiskFileBorn15.Get("h_SigEne");
	TH1D *h_SigEne30   = (TH1D*)DiskFileBorn30.Get("h_SigEne");
	TH1D *h_SigEne100  = (TH1D*)DiskFileBorn100.Get("h_SigEne");
	
	
	
	//collect data from Density object
	TDensity* Density = new TDensity();
	double MH  = Density->m_MH;	
	
	cout << "Creating histograms ..." << endl;
	
	
	TH1D * h_Hist = h_SigBorn;
	h_Hist->SetStats(0);
	h_Hist->SetStats(0);
	h_Hist->SetTitle(0);
	h_Hist->GetYaxis()->CenterTitle();
	h_Hist->GetYaxis()->SetTitleSize(0.07);
	
	#if defined( ELECTRON )
		h_Hist->GetYaxis()->SetTitle( "#sigma(s) [fb]" );
	#elif defined( MUON )
		h_Hist->GetYaxis()->SetTitle("#sigma(s) [pb]");
	#endif 
	
	//h_Hist->GetYaxis()->SetTitle( "#sigma(s) [fb]" );
	h_Hist->GetXaxis()->CenterTitle();
	h_Hist->GetXaxis()->SetTitleSize(0.07);
	h_Hist->GetXaxis()->SetTitle( "#sqrt{s}   [GeV]" );
	h_Hist->GetXaxis()->SetLabelSize(0.055);
	h_Hist->GetYaxis()->SetLabelSize(0.055);
	h_Hist->GetXaxis()->SetTitleOffset(0.72);
	h_Hist->GetYaxis()->SetTitleOffset(0.72);
	
	h_SigBorn->SetLineStyle(1);
	h_SigBorn->SetLineColor(1);
	h_SigEne04->SetLineStyle(2);
	h_SigEne04->SetLineColor(2);
	h_SigEne08->SetLineStyle(3);
	h_SigEne08->SetLineColor(3);
	h_SigEne15->SetLineStyle(4);
	h_SigEne15->SetLineColor(4);
	h_SigEne30->SetLineStyle(6);
	h_SigEne30->SetLineColor(6);
	h_SigEne100->SetLineStyle(7);
	h_SigEne100->SetLineColor(7);
	
	
	//legend
	TLegend * leg = new TLegend(0.6,0.6,0.9,0.9);
	
	leg->AddEntry(h_SigBorn,(string("(Born) #delta = 0") ).c_str(),"l");
	leg->AddEntry(h_SigEne04,(string("(1) #delta = ")+STR(4.2) ).c_str(),"l");
	leg->AddEntry(h_SigEne08,(string("(2) #delta = ")+STR(8)   ).c_str(),"l");
	leg->AddEntry(h_SigEne15,(string("(3) #delta = ")+STR(15)  ).c_str(),"l");  
	leg->AddEntry(h_SigEne30,(string("(4) #delta = ")+STR(30)  ).c_str(),"l");
	leg->AddEntry(h_SigEne100,(string("(5) #delta = ")+STR(100)).c_str(),"l");
	
	leg->SetHeader(" Energy spread [MeV]");
	//leg->SetTextFont( 62 );
	//leg->SetTextSize(12);

	// Caption in latex
	TLatex *CaptTa = new TLatex();
	CaptTa->SetTextAlign(11);
	CaptTa->SetTextSize(0.030);


	TCanvas* canv1 = new TCanvas("canv","plot");	
	//set logscale y
	//gPad->SetLogy();
	
	canv1->cd(1);
	
	h_SigBorn->Draw("l");
	h_SigEne04->Draw("lsame");
	h_SigEne08->Draw("lsame");
	h_SigEne15->Draw("lsame");
	h_SigEne30->Draw("lsame");
	h_SigEne100->Draw("lsame");
	leg->Draw("same");
	
	CaptTa->SetTextAlign(11);
	CaptTa->SetTextColor( 1 );
	CaptTa->DrawLatex(MH, 1.02*h_SigBorn->Interpolate(MH)  ,   "(Born)");
	CaptTa->SetTextColor( 2 );
	CaptTa->DrawLatex(MH, 1.10*h_SigEne04->Interpolate(MH) ,  "(1)");
	CaptTa->SetTextColor( 3 );
	CaptTa->DrawLatex(MH, 1.10*h_SigEne08->Interpolate(MH) ,  "(2)");
	CaptTa->SetTextColor( 4 );
	CaptTa->DrawLatex(MH, 1.10*h_SigEne15->Interpolate(MH) ,  "(3)");
	CaptTa->SetTextColor( 6 );
	CaptTa->DrawLatex(MH, 1.10*h_SigEne30->Interpolate(MH) ,  "(4)");
	CaptTa->SetTextColor( 7 );
	CaptTa->DrawLatex(MH, 1.10*h_SigEne100->Interpolate(MH) , "(5)");
	
	canv1->Update();
	
	canv1->Print( "Bornspread.eps", "" );		
  
    cout << "...creating histograms DONE" << endl;
	
	
	//cleaning
	delete canv1;
	delete leg;
	delete CaptTa;
	delete Density;
	
	
	
};
//----------------------------------------------------------------------


//----------------------------------------------------------------------
/// @brief Calculate convolution of integral with beam dispersion function for ISR and Born
/// @param nbins number of bins in histogram
/// @param NevTot   statistics of integration
int makeISRsigEDistribution( int nbins = 100, long NevTot =  1000000 )
{
	
	//sigma range [MeV]
	double sigE_min = 0.0;
	double sigE_max = 20.0;
	
	int nbt=nbins;
	double dsigE=(sigE_max-sigE_min)/nbt;

	//allocate histograms
		//for E=MH
		TH1D * h_BornMH  = new TH1D("h_BornMH","h_BornMH",nbt,sigE_min,sigE_max);
		h_BornMH->Sumw2();
		
		TH1D * h_ISRMH  = new TH1D("h_ISRMH","h_ISRMH",nbt,sigE_min,sigE_max);
		h_ISRMH->Sumw2();
	
		//for E=MH+GammaH
		TH1D * h_BornMHpGH  = new TH1D("h_BornMHpGH","h_BornMHpGH",nbt,sigE_min,sigE_max);
		h_BornMHpGH->Sumw2();
		
		TH1D * h_ISRMHpGH  = new TH1D("h_ISRMHpGH","h_ISRMHpGH",nbt,sigE_min,sigE_max);
		h_ISRMHpGH->Sumw2();
	
		//for E=MH+GammaH
		TH1D * h_BornMHmGH  = new TH1D("h_BornMHmGH","h_BornMHmGH",nbt,sigE_min,sigE_max);
		h_BornMHmGH->Sumw2(); 
		
		TH1D * h_ISRMHmGH  = new TH1D("h_ISRMHmGH","h_ISRMHmGH",nbt,sigE_min,sigE_max);
		h_ISRMHmGH->Sumw2();
	
	
	//get Higgs mass and width
	TDensity * Density = new TDensity();
	double MH   = Density->m_MH;
	double GamH = Density->m_GamH;

	double sigE;


	//loop over Gaussian spread/sigma
	for(int i=1;i<=nbt;i++)
	{
		sigE=sigE_min+(i-0.5)*dsigE;
		sigE *= 0.001; //MeV to GeV	
		
		//Born histograms
			//do MC integration for given sigE
			TH1D * h_Born = getBornHistoAtSigE( sigE, NevTot );
			//extract values
			double valBornMH = h_Born->Interpolate( MH );
			double valBornMHpGH = h_Born->Interpolate( MH + GamH );
			double valBornMHmGH = h_Born->Interpolate( MH - GamH );
			//fill histograms
			h_BornMH->SetBinContent( i, valBornMH );
			h_BornMHpGH->SetBinContent( i, valBornMHpGH );
			h_BornMHmGH->SetBinContent( i, valBornMHmGH );
		
		//ISR (c) histograms
			//do MC integration for given sigE
			TH1D * h_ISR  = getISRHistoAtSigE( sigE, NevTot );
			//extract values
			double valISRMH  = h_ISR->Interpolate( MH );
			double valISRMHpGH  = h_ISR->Interpolate( MH + GamH );
			double valISRMHmGH  = h_ISR->Interpolate( MH - GamH );
			//fill histograms
			h_ISRMH->SetBinContent( i,valISRMH );
			h_ISRMHpGH->SetBinContent( i,valISRMHpGH );
	 		h_ISRMHmGH->SetBinContent( i,valISRMHmGH );
	 
	 
		cout << i << " of " << nbt <<  ", sigE = " << sigE << ", Born at MH = " << valBornMH << ", ISR at MH = " << valISRMH << endl;
		cout << i << " of " << nbt <<  ", sigE = " << sigE << ", Born at MH + GammaH = " << valBornMHpGH << ", ISR at MH + GammaH = " << valISRMHpGH << endl;
		cout << i << " of " << nbt <<  ", sigE = " << sigE << ", Born at MH - GammaH = " << valBornMHmGH << ", ISR at MH - GammaH = " << valISRMHmGH << endl;
	
		//cleaning
		delete h_Born;
		delete h_ISR;
	
	}
	
	//approximation plots

		//Born normalization
		double BornNormMH = Density->BornH( MH*MH );
		double BornNormMHpGH = Density->BornH( (MH + GamH)*(MH + GamH) );
		double BornNormMHmGH = Density->BornH( (MH - GamH)*(MH - GamH) );
		cout << "BornNormMH = " << BornNormMH << endl;
		cout << "BornNormMHpGH = " << BornNormMHpGH << endl;
		cout << "BornNormMHmGH = " << BornNormMHmGH << endl;
		
		

		// ISR3 normalization create the histogram and get vaues - it assures that the function is self-contained but last longer
			//create temp plot for normalization
				MakeISR( string("tmp-histo-sig0-isr2"), 2, 2, 0.0, NevTot ); 
			//load plot from the file
				TFile DiskFileISR( "./tmp-histo-sig0-isr2.root"); 	
				TH1D *h_SigISR    = (TH1D*)DiskFileISR.Get("h_SigEne");
			//extract normalization values
				double ISRNormMH = h_SigISR->Interpolate( MH );
				double ISRNormMHpGH = h_SigISR->Interpolate( MH + GamH );
				double ISRNormMHmGH = h_SigISR->Interpolate( MH - GamH );
				cout << "ISRNormMH = " << ISRNormMH << endl;	
				cout << "ISRNormMHpGH = " << ISRNormMHpGH << endl;	
				cout << "ISRNormMHmGH = " << ISRNormMHmGH << endl;	
		
		//approx plot at MH
		TH1D * h_approxMH = (TH1D*)h_BornMH->Clone("h_approxMH");
		h_approxMH->Divide( h_ISRMH );
		h_approxMH->Scale( ISRNormMH / BornNormMH );
	
		//approx plot at MH + GammaH
		TH1D * h_approxMHpGH = (TH1D*)h_BornMHpGH->Clone("h_approxMHpGH");
		h_approxMHpGH->Divide( h_ISRMHpGH );
		h_approxMHpGH->Scale( ISRNormMHpGH / BornNormMHpGH );
	
		//approx plot at MH - GammaH
		TH1D * h_approxMHmGH = (TH1D*)h_BornMHmGH->Clone("h_approxMHmGH");
		h_approxMHmGH->Divide( h_ISRMHmGH );
		h_approxMHmGH->Scale( ISRNormMHmGH / BornNormMHmGH );
	
	
	//helper function/closure that set up and save histograms
	struct {
        void operator() ( TH1D * hist, string filename ) const 
        {
			//setup histograms
				TH1D * h_Hist = hist;
				h_Hist->SetStats(0);
				h_Hist->SetStats(0);
				h_Hist->SetTitle(0);
				h_Hist->GetYaxis()->CenterTitle();
				h_Hist->GetYaxis()->SetTitleSize(0.05);
				h_Hist->GetYaxis()->SetTitle("#sigma_{B}^{conv}(E,#delta)/#sigma_{B}(E) #left( #sigma_{(c)}^{conv}(E,#delta)/#sigma_{(c)}^{conv}(E,0) #right)^{-1}");
				h_Hist->GetYaxis()->SetLabelSize(0.05);
				h_Hist->GetXaxis()->CenterTitle();
				h_Hist->GetXaxis()->SetTitleSize(0.05);
				h_Hist->GetXaxis()->SetTitle("#delta [MeV]");
				h_Hist->GetXaxis()->SetLabelSize(0.05);
			
				h_Hist->GetYaxis()->SetRangeUser(0., 2.);
				
			//smooth histograms
				hist->Smooth(5);
		
				
			//save
				TCanvas* canv1 = new TCanvas("canv","plot");
				//canv1->SetGrid(bGrid);
	
				canv1->cd(1);
				hist->Draw("h");
				//hist->Draw("L");
				canv1->Update();
	
				canv1->Print( filename.c_str(), "" );

				delete canv1;
			
        }
    } saveApproxHistograms;	
	
	
	//save approximation histograms
	saveApproxHistograms( h_approxMH, string( "approxVoigtMH.eps" ) );
	saveApproxHistograms( h_approxMHpGH, string( "approxVoigtMHpGH.eps" ) );
	saveApproxHistograms( h_approxMHmGH, string( "approxVoigtMHmGH.eps" ) );
	
	
	
	
	//Voigt distribution plots
		// Born normalization at E=MH
		h_BornMH->Scale( 1.0/BornNormMH );
		h_ISRMH->Scale( 1.0/BornNormMH );

		// Born normalization at E=MH + GammaH
		h_BornMHpGH->Scale( 1.0/BornNormMHpGH );
		h_ISRMHpGH->Scale( 1.0/BornNormMHpGH );

		// Born normalization at E=MH - GammaH
		h_BornMHmGH->Scale( 1.0/BornNormMHmGH );
		h_ISRMHmGH->Scale( 1.0/BornNormMHmGH );



	
	//helper function/closure that set up and save histograms
	struct {
        void operator() ( TH1D * h_ISROFF, TH1D * h_ISRON, string filename ) const 
        {
			//setup histograms
				TH1D * h_Hist = h_ISROFF;
				h_Hist->SetStats(0);
				h_Hist->SetStats(0);
				h_Hist->SetTitle(0);
				h_Hist->GetYaxis()->CenterTitle();
				h_Hist->GetYaxis()->SetTitleSize(0.05);
				h_Hist->GetYaxis()->SetTitle( "" );
				h_Hist->GetYaxis()->SetLabelSize(0.05);
				h_Hist->GetXaxis()->CenterTitle();
				h_Hist->GetXaxis()->SetTitleSize(0.05);
				h_Hist->GetXaxis()->SetTitle("#delta [MeV]");
				h_Hist->GetXaxis()->SetLabelSize(0.05);
			
				h_Hist->GetYaxis()->SetRangeUser(0., 2.);
			
				h_ISROFF->SetLineStyle(1);
				h_ISROFF->SetLineColor(1);
				h_ISRON->SetLineStyle(2);
				h_ISRON->SetLineColor(2);
				
				//smooth histograms
					h_ISROFF->Smooth(5);
					h_ISRON->Smooth(5);	
				
				
				//legend
				TLegend * leg = new TLegend(0.6,0.6,0.9,0.9);
				//leg->SetTextFont( 62 );
				//leg->SetTextSize(12);
	
				leg->AddEntry(h_ISROFF,(string("ISR OFF") ).c_str(),"l");
				leg->AddEntry(h_ISRON,(string("ISR ON") ).c_str(),"l");
	

				TCanvas* canv1 = new TCanvas("canv","plot");
				
				canv1->cd(1);
				h_ISROFF->Draw("h");
				h_ISRON->Draw("hsame");
				//h_ISROFF->Draw("L");
				//h_ISRON->Draw("Lsame");
				leg->Draw();
				canv1->Update();
	
				canv1->Print( filename.c_str(), "" );

				delete canv1;
			
        }
    } saveHistograms;	
	

  
	saveHistograms( h_BornMH, h_ISRMH, string("VoigtMH.eps") );
	saveHistograms( h_BornMHpGH, h_ISRMHpGH, string("VoigtMHpGH.eps") );
	saveHistograms( h_BornMHmGH, h_ISRMHmGH, string("VoigtMHmGH.eps") );
  
  
  
	//ROOT file for histograms
	TFile * rootFile = new TFile( "Voigt.root","RECREATE", "Histograms");
	h_BornMH->Write();
	h_ISRMH->Write();
	h_BornMHpGH->Write();
	h_ISRMHpGH->Write();
	h_BornMHmGH->Write();
	h_ISRMHmGH->Write();
	h_approxMH->Write();
	h_approxMHpGH->Write();
	h_approxMHmGH->Write();
	rootFile->Write();
	rootFile->Close();
	delete rootFile;

	//control printouts
	cout << "-----------------------------------------------------------" << endl;
	cout << "BornConv/Born( delta = 0.5 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(0.5) << endl;
	cout << "BornConv/Born( delta = 1.0 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(1.0) << endl;
	cout << "BornConv/Born( delta = 1.5 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(1.5) << endl;
	cout << "BornConv/Born( delta = 2.0 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(2.0) << endl;
	cout << "BornConv/Born( delta = 2.5 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(2.5) << endl;
	cout << "-----------------------------------------------------------" << endl;
	cout << "ISRConv/Born( delta = 0.5 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(0.5) << endl;
	cout << "ISRConv/Born( delta = 1.0 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(1.0) << endl;
	cout << "ISRConv/Born( delta = 1.5 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(1.5) << endl;
	cout << "ISRConv/Born( delta = 2.0 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(2.0) << endl;
	cout << "ISRConv/Born( delta = 2.5 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(2.5) << endl;
	cout << "-----------------------------------------------------------" << endl;

	//cleaning
	delete h_BornMH;
	delete h_ISRMH;
	delete h_BornMHpGH;
	delete h_ISRMHpGH;
	delete h_BornMHmGH;
	delete h_ISRMHmGH;
	delete h_approxMH;
	delete h_approxMHpGH;
	delete h_approxMHmGH;
	delete Density;

	return 0;

};
//----------------------------------------------------------------------


/// @brief Adjust plots made by makeISRsigEDistribution()
/// @warning you should invoke makeISRsigEDistribution() first
/// @param nbins number of bins in histogram
/// @param NevTot   statistics of integration
int plotISRsigEDistribution( int nbins = 100, long NevTot =  1000000 )
{
	
	//asociate root files - files must exist !!!
	TFile DiskFileVoigt(  "./Voigt.root");   
	TFile DiskFileISR( "./tmp-histo-sig0-isr2.root"); 	

	//load histograms	
	TH1D *h_BornMH    = (TH1D*) DiskFileVoigt.Get("h_BornMH");
	TH1D *h_BornMHpGH    = (TH1D*) DiskFileVoigt.Get("h_BornMHpGH");
	TH1D *h_BornMHmGH    = (TH1D*) DiskFileVoigt.Get("h_BornMHmGH");
	
	TH1D *h_ISRMH    = (TH1D*) DiskFileVoigt.Get("h_ISRMH");
	TH1D *h_ISRMHpGH    = (TH1D*) DiskFileVoigt.Get("h_ISRMHpGH");
	TH1D *h_ISRMHmGH    = (TH1D*) DiskFileVoigt.Get("h_ISRMHmGH");
	
	TH1D *h_approxMH    = (TH1D*) DiskFileVoigt.Get("h_approxMH");
	TH1D *h_approxMHpGH    = (TH1D*) DiskFileVoigt.Get("h_approxMHpGH");
	TH1D *h_approxMHmGH    = (TH1D*) DiskFileVoigt.Get("h_approxMHmGH");
	
	//TH1D *h_SigISR    = (TH1D*)DiskFileISR.Get("h_SigEne");
	
	
	//collect data from Density object
	TDensity* Density = new TDensity();
	double MH   = Density->m_MH;
	//double GamH = Density->m_GamH;
	
	
	cout << "Creating histograms ..." << endl;
	
	
	//helper function/closure that set up and save histograms
	struct {
        void operator() ( TH1D * hist, string filename ) const 
        {
			//setup histograms
				TH1D * h_Hist = hist;
				h_Hist->SetStats(0);
				h_Hist->SetStats(0);
				h_Hist->SetTitle(0);
				h_Hist->GetYaxis()->CenterTitle();
				h_Hist->GetYaxis()->SetTitleSize(0.05);
				h_Hist->GetYaxis()->SetTitle("#sigma_{B}^{conv}(E,#delta)/#sigma_{B}(E) #left( #sigma_{(c)}^{conv}(E,#delta)/#sigma_{(c)}^{conv}(E,0) #right)^{-1}");
				h_Hist->GetYaxis()->SetLabelSize(0.05);
				h_Hist->GetXaxis()->CenterTitle();
				h_Hist->GetXaxis()->SetTitleSize(0.05);
				h_Hist->GetXaxis()->SetTitle("#delta [MeV]");
				h_Hist->GetXaxis()->SetLabelSize(0.05);
				
				h_Hist->GetYaxis()->SetRangeUser(0., 2.);
			
				
			//smooth histograms
				hist->Smooth(5);
		
				
			//save
				TCanvas* canv1 = new TCanvas("canv","plot");
				//canv1->SetGrid(bGrid);
	
				canv1->cd(1);
				hist->Draw("l");
				canv1->Update();
	
				canv1->Print( filename.c_str(), "" );

				delete canv1;
			
        }
    } saveApproxHistograms;	
	
	
	//save approximation histograms
	saveApproxHistograms( h_approxMH, string( "approxVoigtMH.eps" ) );
	saveApproxHistograms( h_approxMHpGH, string( "approxVoigtMHpGH.eps" ) );
	saveApproxHistograms( h_approxMHmGH, string( "approxVoigtMHmGH.eps" ) );
	

	
	//helper function/closure that set up and save histograms
	struct {
        void operator() ( TH1D * h_ISROFF, TH1D * h_ISRON, string filename ) const 
        {
			//setup histograms
				TH1D * h_Hist = h_ISROFF;
				h_Hist->SetStats(0);
				h_Hist->SetStats(0);
				h_Hist->SetTitle(0);
				h_Hist->GetYaxis()->CenterTitle();
				h_Hist->GetYaxis()->SetTitleSize(0.05);
				h_Hist->GetYaxis()->SetTitle( "" );
				h_Hist->GetYaxis()->SetLabelSize(0.05);
				h_Hist->GetXaxis()->CenterTitle();
				h_Hist->GetXaxis()->SetTitleSize(0.05);
				h_Hist->GetXaxis()->SetTitle("#delta [MeV]");
				h_Hist->GetXaxis()->SetLabelSize(0.05);
				
				h_Hist->GetYaxis()->SetRangeUser(0., 2.);
			
				h_ISROFF->SetLineStyle(1);
				h_ISROFF->SetLineColor(1);
				h_ISRON->SetLineStyle(2);
				h_ISRON->SetLineColor(2);
				
				//smooth histograms
					h_ISROFF->Smooth(5);
					h_ISRON->Smooth(5);	
				
				
				//legend
				TLegend * leg = new TLegend(0.6,0.6,0.9,0.9);
				//leg->SetTextFont( 62 );
				//leg->SetTextSize(12);
	
				leg->AddEntry(h_ISROFF,(string("ISR OFF") ).c_str(),"l");
				leg->AddEntry(h_ISRON,(string("ISR ON") ).c_str(),"l");
	

				TCanvas* canv1 = new TCanvas("canv","plot");
				
				canv1->cd(1);
				h_ISROFF->Draw("l");
				h_ISRON->Draw("lsame");
				leg->Draw();
				canv1->Update();
	
				canv1->Print( filename.c_str(), "" );

				delete canv1;
			
        }
    } saveHistograms;	
	

  
	saveHistograms( h_BornMH, h_ISRMH, string("VoigtMH.eps") );
	saveHistograms( h_BornMHpGH, h_ISRMHpGH, string("VoigtMHpGH.eps") );
	saveHistograms( h_BornMHmGH, h_ISRMHmGH, string("VoigtMHmGH.eps") );
  
  
  
	//ROOT file for histograms
	TFile * rootFile = new TFile( "VoigtFinal.root","RECREATE", "Histograms");
	h_BornMH->Write();
	h_ISRMH->Write();
	h_BornMHpGH->Write();
	h_ISRMHpGH->Write();
	h_BornMHmGH->Write();
	h_ISRMHmGH->Write();
	h_approxMH->Write();
	h_approxMHpGH->Write();
	h_approxMHmGH->Write();
	rootFile->Write();
	rootFile->Close();
	delete rootFile;

	//control printouts
	cout << "-----------------------------------------------------------" << endl;
	cout << "BornConv/Born( delta = 0.5 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(0.5) << endl;
	cout << "BornConv/Born( delta = 1.0 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(1.0) << endl;
	cout << "BornConv/Born( delta = 1.5 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(1.5) << endl;
	cout << "BornConv/Born( delta = 2.0 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(2.0) << endl;
	cout << "BornConv/Born( delta = 2.5 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(2.5) << endl;
	cout << "-----------------------------------------------------------" << endl;
	cout << "ISRConv/Born( delta = 0.5 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(0.5) << endl;
	cout << "ISRConv/Born( delta = 1.0 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(1.0) << endl;
	cout << "ISRConv/Born( delta = 1.5 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(1.5) << endl;
	cout << "ISRConv/Born( delta = 2.0 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(2.0) << endl;
	cout << "ISRConv/Born( delta = 2.5 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(2.5) << endl;
	cout << "-----------------------------------------------------------" << endl;

	//cleaning
	delete h_BornMH;
	delete h_ISRMH;
	delete h_BornMHpGH;
	delete h_ISRMHpGH;
	delete h_BornMHmGH;
	delete h_ISRMHmGH;
	delete Density;

	return 0;

};
//----------------------------------------------------------------------
	

/*
/// gives starting index for given process
/// @param N  dimension of array
/// @param workers number of processes
/// @param rank numer of given process
int startIndex( int N, int workers, int rank)
{
	//integer part of division
	int count = N / workers;
	//reminder part of the division
	int remainder = N % workers;
	
	int start;

	if (rank < remainder) 
	{
		start = rank * (count + 1);
	} 
	else 
	{
		start = rank * count + remainder;
	}
	
	return( start );
	
};

/// gives ending index for given process
/// @param N  dimension of array
/// @param workers number of processes
/// @param rank numer of given process
int stopIndex( int N, int workers, int rank)
{
	//integer part of division
	int count = N / workers;
	//reminder part of the division
	int remainder = N % workers;
	
	int start, stop;

	if (rank < remainder) 
	{
		start = rank * (count + 1);
		stop =  start + count;
	} 
	else 
	{
		start = rank * count + remainder;
		stop =  start + count - 1;
	}
	
	return( stop );
	
};
*/


////////////////////////////////////////////////////////////////////////
int main()
{

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);	

	//Monte Carlo statistics, the higher the number, the less the error; error \approx \sqrt{NeVTot}
	long NevTot = 1000000;  //fair
	//long NevTot = 10000000;  //good
	//long NevTot = 100000000; //production

    // Get the number of processes
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // Get the rank of the process
    int iproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
	

	//sigma range [MeV]
	double sigE_min = 0.0;
	double sigE_max = 20.0;
	
	//number of bins in energy spread
	//WARNING - the number of bins have to be multiplicity of the number of processes
	int nbins = 6 * 10;
	
	int nbt = nbins;
	double dsigE=(sigE_max-sigE_min)/nbt;
	
	
	//get Higgs mass and width
	TDensity * Density = new TDensity();
	double MH   = Density->m_MH;
	double GamH = Density->m_GamH;

	//tables for storing E values
	if( nbt % nproc )
		perror("Dimension of table(nbt) should divide by number of processors");
  
  
	int n = nbt / nproc;
	
	
	double * localBornMH = (double * ) malloc(  n * sizeof(double) );
	if( !localBornMH  )
		perror( "Error in creating localBornMH" );
	
	double * globalBornMH = (double * ) malloc( nbt * sizeof(double) );
	if( !globalBornMH  )
		perror( "Error in creating globalBornMH" );
	
	double * localBornMHpGH = (double * ) malloc(  n * sizeof(double) );
	if( !localBornMHpGH  )
		perror( "Error in creating localBornMHpGH" );
	
	double * globalBornMHpGH = (double * ) malloc( nbt * sizeof(double) );
	if( !globalBornMHpGH  )
		perror( "Error in creating globalBornMHpGH" );
	
	double * localBornMHmGH = (double * ) malloc(  n * sizeof(double) );
	if( !localBornMHmGH  )
		perror( "Error in creating localBornMHmGH" );
	
	double * globalBornMHmGH = (double * ) malloc( nbt * sizeof(double) );
	if( !globalBornMHmGH  )
		perror( "Error in creating globalBornMHmGH" );
	
	
	double * localISRMH = (double * ) malloc(  n * sizeof(double) );
	if( !localISRMH  )
		perror( "Error in creating localISRMH" );
	
	double * globalISRMH = (double * ) malloc( nbt * sizeof(double) );
	if( !globalISRMH  )
		perror( "Error in creating globalISRMH" );
	
	double * localISRMHpGH = (double * ) malloc(  n * sizeof(double) );
	if( !localISRMHpGH  )
		perror( "Error in creating localISRMHpGH" );
	
	double * globalISRMHpGH = (double * ) malloc( nbt * sizeof(double) );
	if( !globalISRMHpGH  )
		perror( "Error in creating ISRMHpGH" );
	
	double * localISRMHmGH = (double * ) malloc(  n * sizeof(double) );
	if( !localISRMHmGH  )
		perror( "Error in creating localISRMHmGH" );
	
	double * globalISRMHmGH = (double * ) malloc( nbt * sizeof(double) );
	if( !globalISRMHmGH  )
		perror( "Error in creating ISRMHmGH" );
	


	//loop over Gaussian spread/sigma
	
	//int iDim = 0; 
	int ifirst = n * iproc;
	//int ilast  = ifirst + n;
	int iDim = ifirst; 
	
	double sigE;
	
	//for( int i = ifirst;  i < ilast; i++ )
	for( int i = 0;  i < n; i++ )
	{
		
		sigE=sigE_min+(iDim-0.5)*dsigE;
		sigE *= 0.001; //MeV to GeV	
		
		//Born histograms
			//do MC integration for given sigE
			TH1D * h_Born = getBornHistoAtSigE( sigE, NevTot );
			//extract values
				localBornMH[i]    = h_Born->Interpolate( MH);
				localBornMHpGH[i] = h_Born->Interpolate( MH + GamH );
				localBornMHmGH[i] = h_Born->Interpolate( MH - GamH );		
		
		//ISR (c) histograms
			//do MC integration for given sigE
			TH1D * h_ISR  = getISRHistoAtSigE( sigE, NevTot );
			//extract values
			localISRMH[i]    = h_ISR->Interpolate( MH );
			localISRMHpGH[i] = h_ISR->Interpolate( MH + GamH );
			localISRMHmGH[i] = h_ISR->Interpolate( MH - GamH );


			
		//cleaning
			delete h_Born;
			delete h_ISR;	
				
		iDim++;
	}	
	
	//gather data to process 0
		MPI_Gather( localBornMH, n, MPI_DOUBLE, globalBornMH, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather( localBornMHpGH, n, MPI_DOUBLE, globalBornMHpGH, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather( localBornMHmGH, n, MPI_DOUBLE, globalBornMHmGH, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
		MPI_Gather( localISRMH, n, MPI_DOUBLE, globalISRMH, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather( localISRMHpGH, n, MPI_DOUBLE, globalISRMHpGH, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather( localISRMHmGH, n, MPI_DOUBLE, globalISRMHmGH, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	
	//fill histograms in 0 process
	if( iproc == 0 )
	{
		
		//allocate histograms
			//for E=MH
			TH1D * h_BornMH  = new TH1D("h_BornMH","h_BornMH",nbt,sigE_min,sigE_max);
			h_BornMH->Sumw2();
		
			TH1D * h_ISRMH  = new TH1D("h_ISRMH","h_ISRMH",nbt,sigE_min,sigE_max);
			h_ISRMH->Sumw2();
	
			//for E=MH+GammaH
			TH1D * h_BornMHpGH  = new TH1D("h_BornMHpGH","h_BornMHpGH",nbt,sigE_min,sigE_max);
			h_BornMHpGH->Sumw2();
		
			TH1D * h_ISRMHpGH  = new TH1D("h_ISRMHpGH","h_ISRMHpGH",nbt,sigE_min,sigE_max);
			h_ISRMHpGH->Sumw2();
	
			//for E=MH+GammaH
			TH1D * h_BornMHmGH  = new TH1D("h_BornMHmGH","h_BornMHmGH",nbt,sigE_min,sigE_max);
			h_BornMHmGH->Sumw2(); 
		
			TH1D * h_ISRMHmGH  = new TH1D("h_ISRMHmGH","h_ISRMHmGH",nbt,sigE_min,sigE_max);
			h_ISRMHmGH->Sumw2();
		
		
		
		
			
		//double sigE;
		
		for(int i=1;i<=nbt;i++)
		{
			sigE=sigE_min+(i-0.5)*dsigE;
			sigE *= 0.001; //MeV to GeV	
		
		//Born histograms
			h_BornMH->SetBinContent( i, globalBornMH[i] );
			h_BornMHpGH->SetBinContent( i, globalBornMHpGH[i] );
			h_BornMHmGH->SetBinContent( i, globalBornMHmGH[i] );
		
		//ISR (c) histograms			
			h_ISRMH->SetBinContent( i, globalISRMH[i] );
			h_ISRMHpGH->SetBinContent( i, globalISRMHpGH[i] );
	 		h_ISRMHmGH->SetBinContent( i, globalISRMHmGH[i] );
	 
	 /*
		cout << i << " of " << nbt <<  ", sigE = " << sigE << ", Born at MH = " << valBornMH << ", ISR at MH = " << valISRMH << endl;
		cout << i << " of " << nbt <<  ", sigE = " << sigE << ", Born at MH + GammaH = " << valBornMHpGH << ", ISR at MH + GammaH = " << valISRMHpGH << endl;
		cout << i << " of " << nbt <<  ", sigE = " << sigE << ", Born at MH - GammaH = " << valBornMHmGH << ", ISR at MH - GammaH = " << valISRMHmGH << endl;
	*/
			
		
		};
		
	
		//approximation plots

			//Born normalization
				double BornNormMH = Density->BornH( MH*MH );
				double BornNormMHpGH = Density->BornH( (MH + GamH)*(MH + GamH) );
				double BornNormMHmGH = Density->BornH( (MH - GamH)*(MH - GamH) );
				cout << "BornNormMH = " << BornNormMH << endl;
				cout << "BornNormMHpGH = " << BornNormMHpGH << endl;
				cout << "BornNormMHmGH = " << BornNormMHmGH << endl;
		
		

			// ISR3 normalization create the histogram and get vaues - it assures that the function is self-contained but last longer
				//create temp plot for normalization
					MakeISR( string("tmp-histo-sig0-isr2"), 2, 2, 0.0, NevTot ); 
				//load plot from the file
					TFile DiskFileISR( "./tmp-histo-sig0-isr2.root"); 	
					TH1D *h_SigISR    = (TH1D*)DiskFileISR.Get("h_SigEne");
				//extract normalization values
					double ISRNormMH = h_SigISR->Interpolate( MH );
					double ISRNormMHpGH = h_SigISR->Interpolate( MH + GamH );
					double ISRNormMHmGH = h_SigISR->Interpolate( MH - GamH );
					cout << "ISRNormMH = " << ISRNormMH << endl;	
					cout << "ISRNormMHpGH = " << ISRNormMHpGH << endl;	
					cout << "ISRNormMHmGH = " << ISRNormMHmGH << endl;	
		
			//approx plot at MH
			TH1D * h_approxMH = (TH1D*)h_BornMH->Clone("h_approxMH");
			h_approxMH->Divide( h_ISRMH );
			h_approxMH->Scale( ISRNormMH / BornNormMH );
	
			//approx plot at MH + GammaH
			TH1D * h_approxMHpGH = (TH1D*)h_BornMHpGH->Clone("h_approxMHpGH");
			h_approxMHpGH->Divide( h_ISRMHpGH );
			h_approxMHpGH->Scale( ISRNormMHpGH / BornNormMHpGH );
	
			//approx plot at MH - GammaH
			TH1D * h_approxMHmGH = (TH1D*)h_BornMHmGH->Clone("h_approxMHmGH");
			h_approxMHmGH->Divide( h_ISRMHmGH );
			h_approxMHmGH->Scale( ISRNormMHmGH / BornNormMHmGH );
	
	
		//helper function/closure that set up and save histograms
		struct {
			void operator() ( TH1D * hist, string filename ) const 
			{
				//setup histograms
					TH1D * h_Hist = hist;
					h_Hist->SetStats(0);
					h_Hist->SetStats(0);
					h_Hist->SetTitle(0);
					h_Hist->GetYaxis()->CenterTitle();
					h_Hist->GetYaxis()->SetTitleSize(0.05);
					h_Hist->GetYaxis()->SetTitle("#sigma_{B}^{conv}(E,#delta)/#sigma_{B}(E) #left( #sigma_{(c)}^{conv}(E,#delta)/#sigma_{(c)}^{conv}(E,0) #right)^{-1}");
					h_Hist->GetYaxis()->SetLabelSize(0.05);
					h_Hist->GetXaxis()->CenterTitle();
					h_Hist->GetXaxis()->SetTitleSize(0.05);
					h_Hist->GetXaxis()->SetTitle("#delta [MeV]");
					h_Hist->GetXaxis()->SetLabelSize(0.05);
			
					h_Hist->GetYaxis()->SetRangeUser(0., 2.);
				
				//smooth histograms
					hist->Smooth(5);
		
				
				//save
					TCanvas* canv1 = new TCanvas("canv","plot");
					//canv1->SetGrid(bGrid);
	
					canv1->cd(1);
					hist->Draw("h");
					//hist->Draw("L");
					canv1->Update();
	
					canv1->Print( filename.c_str(), "" );

					delete canv1;
			
			}
		} saveApproxHistograms;	
	
	
		//save approximation histograms
		saveApproxHistograms( h_approxMH, string( "approxVoigtMH.eps" ) );
		saveApproxHistograms( h_approxMHpGH, string( "approxVoigtMHpGH.eps" ) );
		saveApproxHistograms( h_approxMHmGH, string( "approxVoigtMHmGH.eps" ) );
	
	
	
	
		//Voigt distribution plots
			// Born normalization at E=MH
			h_BornMH->Scale( 1.0/BornNormMH );
			h_ISRMH->Scale( 1.0/BornNormMH );

			// Born normalization at E=MH + GammaH
			h_BornMHpGH->Scale( 1.0/BornNormMHpGH );
			h_ISRMHpGH->Scale( 1.0/BornNormMHpGH );

			// Born normalization at E=MH - GammaH
			h_BornMHmGH->Scale( 1.0/BornNormMHmGH );
			h_ISRMHmGH->Scale( 1.0/BornNormMHmGH );



	
		//helper function/closure that set up and save histograms
		struct {
        void operator() ( TH1D * h_ISROFF, TH1D * h_ISRON, string filename ) const 
			{
				//setup histograms
					TH1D * h_Hist = h_ISROFF;
					h_Hist->SetStats(0);
					h_Hist->SetStats(0);
					h_Hist->SetTitle(0);
					h_Hist->GetYaxis()->CenterTitle();
					h_Hist->GetYaxis()->SetTitleSize(0.05);
					h_Hist->GetYaxis()->SetTitle( "" );
					h_Hist->GetYaxis()->SetLabelSize(0.05);
					h_Hist->GetXaxis()->CenterTitle();
					h_Hist->GetXaxis()->SetTitleSize(0.05);
					h_Hist->GetXaxis()->SetTitle("#delta [MeV]");
					h_Hist->GetXaxis()->SetLabelSize(0.05);
			
					h_Hist->GetYaxis()->SetRangeUser(0., 2.);
			
					h_ISROFF->SetLineStyle(1);
					h_ISROFF->SetLineColor(1);
					h_ISRON->SetLineStyle(2);
					h_ISRON->SetLineColor(2);
				
					//smooth histograms
						h_ISROFF->Smooth(5);
						h_ISRON->Smooth(5);	
				
				
					//legend
						TLegend * leg = new TLegend(0.6,0.6,0.9,0.9);
						//leg->SetTextFont( 62 );
						//leg->SetTextSize(12);
		
						leg->AddEntry(h_ISROFF,(string("ISR OFF") ).c_str(),"l");
						leg->AddEntry(h_ISRON,(string("ISR ON") ).c_str(),"l");
	

					TCanvas* canv1 = new TCanvas("canv","plot");
				
					canv1->cd(1);
					h_ISROFF->Draw("h");
					h_ISRON->Draw("hsame");
					//h_ISROFF->Draw("L");
					//h_ISRON->Draw("Lsame");
					leg->Draw();
					canv1->Update();
	
					canv1->Print( filename.c_str(), "" );

					delete canv1;
			
			}
		} saveHistograms;	
	

  
		saveHistograms( h_BornMH, h_ISRMH, string("VoigtMH.eps") );
		saveHistograms( h_BornMHpGH, h_ISRMHpGH, string("VoigtMHpGH.eps") );
		saveHistograms( h_BornMHmGH, h_ISRMHmGH, string("VoigtMHmGH.eps") );
  
  
  
		//ROOT file for histograms
		TFile * rootFile = new TFile( "Voigt.root","RECREATE", "Histograms");
		h_BornMH->Write();
		h_ISRMH->Write();
		h_BornMHpGH->Write();
		h_ISRMHpGH->Write();
		h_BornMHmGH->Write();
		h_ISRMHmGH->Write();
		h_approxMH->Write();
		h_approxMHpGH->Write();
		h_approxMHmGH->Write();
		rootFile->Write();
		rootFile->Close();
		delete rootFile;

		//control printouts
		cout << "-----------------------------------------------------------" << endl;
		cout << "BornConv/Born( delta = 0.5 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(0.5) << endl;
		cout << "BornConv/Born( delta = 1.0 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(1.0) << endl;
		cout << "BornConv/Born( delta = 1.5 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(1.5) << endl;
		cout << "BornConv/Born( delta = 2.0 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(2.0) << endl;
		cout << "BornConv/Born( delta = 2.5 )( E = "  << MH  << ") = " << h_BornMH->Interpolate(2.5) << endl;
		cout << "-----------------------------------------------------------" << endl;
		cout << "ISRConv/Born( delta = 0.5 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(0.5) << endl;
		cout << "ISRConv/Born( delta = 1.0 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(1.0) << endl;
		cout << "ISRConv/Born( delta = 1.5 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(1.5) << endl;
		cout << "ISRConv/Born( delta = 2.0 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(2.0) << endl;
		cout << "ISRConv/Born( delta = 2.5 )( E = "  << MH  << ") = " << h_ISRMH->Interpolate(2.5) << endl;
		cout << "-----------------------------------------------------------" << endl;

		//cleaning
		delete h_BornMH;
		delete h_ISRMH;
		delete h_BornMHpGH;
		delete h_ISRMHpGH;
		delete h_BornMHmGH;
		delete h_ISRMHmGH;
		delete h_approxMH;
		delete h_approxMHpGH;
		delete h_approxMHmGH;
		delete Density;

	}



	//synchronization barier
	MPI_Barrier( MPI_COMM_WORLD ) ;


	if( iproc == 0 )
	{
		cout << "Making plots ...." << endl;
	
			plotISRsigEDistribution();
	
		cout << "... making plots DONE" << endl;
	}

	//synchronization barier
	MPI_Barrier( MPI_COMM_WORLD ) ;

	cout << "Goodbye from process " << iproc << endl;
	
	
	if( iproc == 0 )
	{
		cout << "***** End of Demonstration Program  *****" << endl;
	}
	
	//cleaning
		free( localBornMH );
		free( globalBornMH );
	
		free( localBornMHpGH );
		free( globalBornMHpGH );
	
		free( localBornMHmGH );
		free( globalBornMHmGH );
	
		free( localISRMH );
		free( globalISRMH );
	
		free( localISRMHpGH );
		free( globalISRMHpGH );
	
		free( localISRMHmGH );
		free( globalISRMHmGH );
		
	
	// Finalize the MPI environment.
		MPI_Finalize();
	
	return 0;
	
};

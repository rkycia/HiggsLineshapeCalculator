#ifndef TDensity_H_
#define TDensity_H_


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


#define ELECTRON
//#define MUON


///Class density for Foam
class TDensity: public TFoamIntegrand 
{
public:
	
	///Construcotr
	TDensity();
	
	///Destructor
	virtual ~TDensity(){;};
	
	Double_t Density(int nDim, Double_t *Xarg);
	
	Double_t DensityISR(int nDim, Double_t *Xarg);
	
	Double_t DensityBorn(int nDim, Double_t *Xarg);
	
	/// ON/OFF ISR corrections nonzero/0, default ON
	int    m_ISROn;		
	/// =2 for Bremss, =3 for energy spread
	int    m_kDim;      
	/// Type of ISR/QED switch a-0, b-1, c-2
	int    m_KeyISR;    
	/// No. of cells, optional, default=10000
	int    m_nCells;    
	/// No. of MC evts/cell in exploration, default=10000
	int    m_nSampl;    
	
	/// machine spread of sqrt(s)
	double m_sigE;      
		
	/// Energy used to calculate weight		
	double m_E;			
	/// y value used to calculate weight
	double m_y;			
	
	/// GeV^2 to nanobarn conversion constant
	double m_gnanob;     
	/// Ludolphian number
	double m_pi;         
	///  Euler-Mascheroni constant
	double m_ceuler;     
	/// inverse of QED coupling constant
	double m_alfinv;     
	/// QED coubling constant divided by pi
	double m_alfpi;      
	/// electron mass
	double m_amel;   
	
	/// Higgs mass
	double m_MH;  
	/// Higgs width            
	double m_GamH;
	
	///Low bound of energy integration
	double m_Emin;		
	///High bound of energy integration
	double m_Emax;		
	///Current energy value used to weight calculation
	double m_Ene;		
	
	
	/// @returns Born term for given s = E * E
	/// @param s = E * E
	double BornH( double s );
  
	/// @returns ISR distribution
	/// @param svar s in CM
	/// @param vv   v varaible
	double RhoISR( double svar, double vv );
  
	/// @returns gamma from the paper
	/// @param svar = E*E
	double gamISR( double svar){ return  2*m_alfpi*( log(svar/sqr(m_amel)) -1);}
	
	/// @returns YFS formfactor
	/// @param svar = E*E
	double FFact( double svar)
	{
		//YFS formfactor
		double beti  = gamISR(svar);
		return  exp(-m_ceuler*beti)/TMath::Gamma(1+beti) * exp( beti/4 +m_alfpi *(-0.5  +sqr(m_pi)/3.0) );
	}
  
	///@returns square of x
	double sqr( double x ){ return x*x; };

	
private:
	
	//nothing
    
};
 
#endif

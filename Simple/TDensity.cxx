#include "TDensity.h"


TDensity::TDensity()
{
	
	// all defaults defined here can be changed by the user
	// before calling TFoam::Initialize()

	m_gnanob  = 389.37966e3;
	m_pi      = 3.1415926535897932;
	m_ceuler  = 0.57721566;
	m_alfinv  = 137.035;
	m_alfpi   = 1/m_alfinv/m_pi;
	
	#if defined( ELECTRON )
		m_amel    = 0.510999e-3;
	#elif defined( MUON )
		m_amel    = 105.658e-3;
	#endif	
	
	m_MH      =  125.09;     // Higgs mass
	m_GamH    =  0.0042;    // Higgs width
	
	m_kDim    = 2;        // =2,3 energy spread off/on
	m_sigE    = 0.005;    // machine spread of sqrt(s)

	m_kDim    =    3;         // No. of dim. for Foam, =2,3 Machine energy spread OFF/ON
	m_nCells  = 10000;         // No. of cells, optional, default=2000
	m_nSampl  = 10000;         // No. of MC evts/cell in exploration, default=200

	m_ISROn	  = 1;	// ON ISR
	m_KeyISR  = 2;  // Type of ISR/QED switch
	
	//bounds of integration over E +- 10 sigma
	m_Emin=m_MH - 10.0 *m_GamH;
	m_Emax=m_MH + 10.0 *m_GamH;

};
//----------------------------------------------------------------------


double TDensity::BornH( double s) 
{
	// SqChiZ= s**2/((s-MZ2)**2+(GammZ*s/MZ)**2) *RaZ**2 ! variable width
	// double bornZ= s/( sqr(s-sqr(m_MZ))+ sqr(m_GamZ*s/m_MZ)); // variable width


	#if defined( ELECTRON )
		double BrHee = 5.3e-9; //e e branching ratio
	#elif defined( MUON )
	  double BrHee = 2.19e-4; //mu mu branching ratio
	#endif 	
	
	double bornH= 1/( sqr(s-sqr(m_MH))+ sqr(m_GamH*m_MH) );
	//double bornH= 1/( sqr(s-sqr(m_MH))+ sqr(m_GamH*s/m_MH) );
	bornH *= 4*m_pi*BrHee *sqr(m_GamH);  // Higgs total
	
	#if defined( ELECTRON )
		bornH *= m_gnanob*1e6; // femtobarn
	#elif defined( MUON )
		bornH *= m_gnanob*1e3; // picobarn
	#endif
	
	return bornH;
	
};
//----------------------------------------------------------------------


double TDensity::RhoISR(double svar, double vv)
{

	double alf1   = m_alfpi;
	double beti   = gamISR(svar);
	///
	double gamfac = exp(-m_ceuler*beti)/TMath::Gamma(1+beti);
	double delb   = beti/4 +alf1*(-0.5  +sqr(m_pi)/3.0);
	double ffact  = gamfac*exp(delb);


	double rho,dels,delh;
	if(       m_KeyISR == 0)
	{
		// zero   order exponentiated
		dels = 0;
		delh = 0;
		//rho  = beti* exp( log(vv)*(beti-1) ) *(1 +dels +delh);
		rho  = ffact*beti* exp( log(vv)*(beti-1) ) *(1 +dels +delh);
	}
	else if( m_KeyISR == 1)
	{
		// first  order
		dels = beti/2;   // NLO part =0 as for vector boson
		delh = vv*(-1 +vv/2);
		rho = ffact*beti* exp( log(vv)*(beti-1) ) *(1 +dels +delh);
	}	
	else if( m_KeyISR == 2)
	{
		// second order without NLO part
		dels = beti/2 +sqr(beti)/8;
		delh = vv*(-1+vv/2.0)+beti*0.5*(-0.25*(4.0-6.0*vv+3.0*vv*vv)*log(1-vv)-vv);
		rho = ffact*beti* exp( log(vv)*(beti-1) ) *(1 +dels +delh);
	}
	else
	{
	  cout<<"+++++TDensity::RhoISR: Wrong KeyISR = " << m_KeyISR<<endl;
	  exit(5);
	}
	
	return rho;

};
//----------------------------------------------------------------------


///Density function to integrate by Foam. It represents convoluted Born only.
/// @warning Requires nDim = 2.
/// @param nDim dimenison of integration
/// @param Xarg vector of random numebrs of dimension nDim
double TDensity::DensityBorn(int nDim, double *Xarg)
{
	assert( nDim == 2 && "DensityBorn() convolution with nDim != 2" );
	
	double Dist=1;

	// Machine energy
	m_Ene = m_Emin +(m_Emax-m_Emin)*Xarg[0]; // Mapping
	Dist *= (m_Emax-m_Emin);                 // Jacobian
	m_E = m_Ene;

	// Machine energy spread
	double DelE=0;
	DelE = -10*m_sigE +20*m_sigE*Xarg[1];
	Dist *= 20*m_sigE;
	Dist *= 1/(m_sigE*sqrt(2*m_pi))*exp(-0.5*(DelE/m_sigE)*(DelE/m_sigE) );
	
	// Born at reduced/smeared energy
	double Ene = m_Ene +DelE;
	double svar = Ene*Ene;
	Dist *= BornH( svar );

	//cout << "Density = " << Dist << endl;

	return Dist;

};
//----------------------------------------------------------------------



///Density function to integrate by Foam. It calculates ISR and convoluted ISR.
/// @param nDim dimenison of integration
/// @param Xarg vector of random numebrs of dimension nDim
double TDensity::DensityISR(int nDim, double *Xarg)
{
	double Dist=1;

	// Machine energy
	m_Ene = m_Emin +(m_Emax-m_Emin)*Xarg[0]; // Mapping
	Dist *= (m_Emax-m_Emin);                 // Jacobian
	m_E = m_Ene;

	// ISR energy loss
	double eps = 1e-32;
	double ymin =log(eps);
	double ymax =0;
	m_y = ymin+(ymax-ymin)*Xarg[1]; // Mapping
	Dist *= (ymax-ymin);            // Jacobian
	double v=exp(m_y);              // Mapping
	Dist *= v;                      // Jacobian

	// Machine energy spread
	double DelE=0;
	if(m_kDim == 3)
	{
		DelE = -10*m_sigE +20*m_sigE*Xarg[2];
		Dist *= 20*m_sigE;
		Dist *= 1/(m_sigE*sqrt(2*m_pi))*exp(-0.5*(DelE/m_sigE)*(DelE/m_sigE) );
	}

	// Born at reduced/smeared energy
	double Ene = m_Ene +DelE;
	double svar = Ene*Ene;
	Dist *= BornH(svar*(1-v));

	Dist *= RhoISR(svar,v);    // ISR loss functiom

	//cout << "Density = " << Dist << endl;

	return Dist;

};
//----------------------------------------------------------------------


double TDensity::Density(int nDim, double *Xarg)
{
	double Dist=1;

	if( m_ISROn > 0 )
	{
		Dist = DensityISR( nDim, Xarg );
	}
	else if( m_ISROn == 0)
	{
		Dist = DensityBorn( nDim, Xarg );
	}
	else 
	{
		cout << "TDensity::Density: wrong m_ISROn value." << endl;
		exit( 1 );
	};	
	
	
	//cout << "Density = " << Dist << endl;

	return Dist;

};
//----------------------------------------------------------------------

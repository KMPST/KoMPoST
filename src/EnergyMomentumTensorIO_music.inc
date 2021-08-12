/*
 * Copyright (c) 2018, Aleksi Kurkela, Aleksas Mazeliauskas, Jean-Francois
 * Paquet, Soeren Schlichting and Derek Teaney
 * All rights reserved.
 *
 * KoMPoST is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/KMPST/KoMPoST/
 */
#include <fstream>
#include "EventInput.h"
// routines for finding the hydro fields from TMUNU
#include "EnergyMomentumDecomposition.inc"
// routines for calculating Navier-Stokes estimate of PIMUNU
#include "ComputeSigma.inc"

// EnergyMomentumTensoIO_music.inc defines read in/out routines loading/writing out
// different components of EnergyMomentumTensorMap data.
// It should be replaced by your own favourite read in/out routine.


// SAVE TMUNU TO FILE //
// File format
//
// for iy in Ns:
//     for ix in Ns:
//         1:x 2:y  
//         3:T^{tautau} 4:T^{xx} 5:T^{yy} 6:tau^2 T^{etaeta} 
//         7:-T^{tau x} 8:-T^{tauy} 9:-tau T^{taueta} 
//         10:-T^{xy} 11:-tau T^{yeta} 12:-tau T^{etax} 
//         13-28:CellData[0 -- 15] -- additonal cell data. See KineticEvolution.cpp for definitions
//
// where ix and iy are cell indecies with grid parameters given by EventInput
// with grid spacing afm in [fm] and grid size [Ns x Ns]
// x, y given in [fm]
//       
// T^{munu} writen out in units [Gev/fm^3]
void EnergyMomentumTensorMapSave(EnergyMomentumTensorMap *TOut,std::string fname){
    
    using namespace EventInput ;

    std::ofstream OutStream;
    OutStream.open(fname.c_str());
    
    std::cerr << "#SAVING ENERGY-MOMENTUM TENSOR AT " << TOut->tau << " fm/c TO " << fname << std::endl;
    
    // Optionally CONVERT  E [GeV^4] TO E [GeV/fm^3] //
    double ScaleFactor= 1./(M_HBARC*M_HBARC*M_HBARC) ;
    
    for(int yS=ySTART;yS<=yEND;yS++){
        for(int xS=xSTART;xS<=xEND;xS++){
        
            // GET RELEVANT COMPONENTS //
            double T00=TOut->Get(0,0,xS,yS);
            double TXX=TOut->Get(1,1,xS,yS);
            double TYY=TOut->Get(2,2,xS,yS);
            double TZZ=TOut->Get(3,3,xS,yS);
            
            double T0X=TOut->Get(0,1,xS,yS);
            double T0Y=TOut->Get(0,2,xS,yS);
            double T0Z=TOut->Get(0,3,xS,yS);
            
            double TXY=TOut->Get(1,2,xS,yS);
            double TYZ=TOut->Get(2,3,xS,yS);
            double TZX=TOut->Get(3,1,xS,yS);

            double x = (xS-Ns/2)*afm ;
            double y = (yS-Ns/2)*afm ;

            // CREATE OUTPUT //
            OutStream << x << " " << y << " ";
            OutStream << T00*ScaleFactor << " " << TXX*ScaleFactor << " " << TYY*ScaleFactor << " " << TZZ*ScaleFactor << " ";
            OutStream << T0X*ScaleFactor << " " << T0Y*ScaleFactor << " " << T0Z*ScaleFactor << " ";
            OutStream << TXY*ScaleFactor << " " << TYZ*ScaleFactor << " " << TZX*ScaleFactor << " ";

            // Additional CellData. See KineticEvolution.cpp for definitions 
            for (int id = 0 ; id < EnergyMomentumTensorMap::NCellData; id++) {
                OutStream << TOut->GetCellData(id,xS,yS) << " " ;
            }
            

            OutStream << std::endl;
            
        }
        OutStream << std::endl;
    }
    
    // CLOSE //
    OutStream.close();
    
    // COMMANDLINE OUTPUT //
    std::cerr << "#DONE" << std::endl;
    
}




// LOAD TMUNU FROM FILE //
// File format
//
// for iy in Ns:
//     for ix in Ns:
//         1:ix 2:iy  
//         3:T^{tautau} 4:T^{xx} 5:T^{yy} 6:tau^2 T^{etaeta} 
//         7:-T^{tau x} 8:-T^{tauy} 9:-tau T^{taueta} 
//         10:-T^{xy} 11:-tau T^{yeta} 12:-tau T^{etax}
//
// where ix and iy are cell indecies with grid parameters given by EventInput
// with grid spacing afm in [fm] and grid size [Ns x Ns]
//       
// T^{munu} given in units Gev/fm^3, but saved in units of [GeV^4]
void EnergyMomentumTensorMapLoad(EnergyMomentumTensorMap *TIn, std::string fname){
    using namespace EventInput ;
        
    // OPEN INPUT FILE //
    std::ifstream InStream;
    InStream.open(fname.c_str());
    
    // COMMANDLINE OUTPUT //
    std::cerr << "#GETTING ENERGY-MOMENTUM TENSOR FROM FILE " << fname << std::endl;
    
    // INPUT LINE STRING //
    std::string InputLine;
    
    // SPATIAL COORDINATES //
    int xS,yS;
    
    // COMPONENTS OF ENERGY MOMENTUM TENSOR//
    double T00,TXX,TYY,TZZ,T0X,T0Y,T0Z,TXY,TYZ,TZX;
    
    // CONVERT E [GeV/fm^3] TO  E [GeV^4] //
    double ScaleFactor=(M_HBARC*M_HBARC*M_HBARC);
    
    int InputCount=0;
    bool printfirst=true; 
    // GET DATA LINE BY LINE //
    // Get header first
    getline(InStream,InputLine);
    while(InStream.good()){
        
        // READ LINE //
        getline(InStream,InputLine);
        
        if(!(InputLine.empty())){
            
            std::stringstream InputValues(InputLine);
            
            // GET COMPONENTS OF ENERGY-MOMENTUM TENSOR //
            InputValues >> xS; InputValues >> yS;
            
            InputValues >> T00; InputValues >> TXX; InputValues >> TYY; InputValues >> TZZ;
            InputValues >> T0X; InputValues >> T0Y; InputValues >> T0Z;
            InputValues >> TXY; InputValues >> TYZ; InputValues >> TZX;
            
            // SAVE ENERGY-MOMENTUM TENSOR TO MEMORY //
            const double small_num=1000*GSL_DBL_MIN;
            if (T00 < small_num) {
		    if (printfirst) {
			    std::cerr << "Warning: encountered small energy density T^{\\tau\\tau}<" << small_num << " at (ix,iy)=" << xS << "," << yS << ")!! Using T^{\\mu\\nu}=diag(" << 10*small_num << "," << 5*small_num << "," << 5*small_num << ",0) instead.\n"; 
			    std::cerr << "Further warnings of small T^{\\tau\\tau} will be suppressed." << std::endl;
			    printfirst=false;
		    }
                    TIn->Set(xS,yS,10*small_num,5*small_num,5*small_num,0,0,0,0,0,0,0);
            }
            else {
                    TIn->Set(xS,yS,ScaleFactor*T00,ScaleFactor*TXX,ScaleFactor*TYY,ScaleFactor*TZZ,ScaleFactor*T0X,ScaleFactor*T0Y,ScaleFactor*T0Z,ScaleFactor*TXY,ScaleFactor*TYZ,ScaleFactor*TZX);
            }

            // COUNTER //
            InputCount++;
            
            
        }
        
    }
    
    // CHECK //
    if(InputCount!=Ns*Ns){
        std::cerr << "#ERROR -- COULD NOT LOAD INPUT FILE CORRECTLY" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // CLOSE //
    InStream.close();
    
    // COMMANDLINE OUTPUT //
    std::cerr << "#DONE" << std::endl;
    
    
}
    
    
//Take Tmunu from EnergyMomentumTensorMap class and return T^{munu} (hyperbolic coordinates)
void get_Tmunu_raised_hyperbolic(class EnergyMomentumTensorMap * Tmunu, int ix, int iy, double Tmunu_raised[4][4]) {

	//tau
	const double tau_in_fm=Tmunu->tau;
	const double tau_in_GeV=tau_in_fm/(M_HBARC);
	const double tau_inv_in_GeV=1./tau_in_GeV;

	const double Ttautau=Tmunu->Get(0,0,ix,iy);
	const double Ttaux=-1*Tmunu->Get(0,1,ix,iy);
	const double Ttauy=-1*Tmunu->Get(0,2,ix,iy);
	const double Ttaueta=-1*tau_inv_in_GeV*Tmunu->Get(0,3,ix,iy);
	const double Txx=Tmunu->Get(1,1,ix,iy);
	const double Txy=-1*Tmunu->Get(1,2,ix,iy);
	const double Txeta=-1*tau_inv_in_GeV*Tmunu->Get(1,3,ix,iy);
	const double Tyy=Tmunu->Get(2,2,ix,iy);
	const double Tyeta=-1*tau_inv_in_GeV*Tmunu->Get(2,3,ix,iy);
	const double Tetaeta=tau_inv_in_GeV*tau_inv_in_GeV*Tmunu->Get(3,3,ix,iy);

	Tmunu_raised[0][0]=Ttautau;
	Tmunu_raised[0][1]=Ttaux;
	Tmunu_raised[1][0]=Ttaux;
	Tmunu_raised[0][2]=Ttauy;
	Tmunu_raised[2][0]=Ttauy;
	Tmunu_raised[0][3]=Ttaueta;
	Tmunu_raised[3][0]=Ttaueta;
	Tmunu_raised[1][1]=Txx;
	Tmunu_raised[1][2]=Txy;
	Tmunu_raised[2][1]=Txy;
	Tmunu_raised[1][3]=Txeta;
	Tmunu_raised[3][1]=Txeta;
	Tmunu_raised[2][2]=Tyy;
	Tmunu_raised[2][3]=Tyeta;
	Tmunu_raised[3][2]=Tyeta;
	Tmunu_raised[3][3]=Tetaeta;
}

// Get conformal entropy from energy density (given in [GeV^4
double get_conformal_entropy(double energy_density) {
	const double p=energy_density/3.0;
	const double T=pow(p/(16*M_PI*M_PI/90),1./4);	
	const double entropy=(energy_density+p)/T;
	return entropy;
};



// Write initial condition file that can be read by MUSIC
// for ix in Ns:
//     for iy in Ns:
//         1:ieta 2:ix 3:iy  
//         4:e [GeV/fm^3] 5:u^tau 6:u^x 7:u^y 8:u^eta [fm^{-1}] 
//         9: pit^{tautau} 10:pi^{taux} 11:pi^{tauy} 12: pi^{taueta} 
//         13:pi^{xx} 14:pi^{xy} 15:pi^{xeta} 16:pi^{yy} 17:pi^{yeta} 18:pi^{etaeta}
//
// where ieta=0, and ix and iy are cell indecies with grid parameters given by EventInput
// with grid spacing afm in [fm] and grid size [Ns x Ns]
//       
//\pi^{\mu\nu} is in units of fermi.
// Note that the loop order is different than in the input file
void write_initial_conditions_MUSIC(std::string outfile_name, bool use_sigmamunu_NavierStokes, double dx, double dy, int Nx, int Ny, EnergyMomentumTensorMap *TOut_Full, EnergyMomentumTensorMap *TOut_BG) {
  //void write_initial_conditions_MUSIC(std::string outfile_name, double tau_in_fm, double dx, double dy, int Nx, int Ny, double ** energy_density, double *** flow, double **** shear_tensor, double ** bulk_pressure) {

  std::cerr << "#SAVING " << outfile_name.c_str() << std::endl;

  std::fstream outfile;
  outfile.open(outfile_name.c_str(),std::fstream::out);

  //From GeV^4 to fm^-4
  const double shear_tensor_normalisation=1./(M_HBARC*M_HBARC*M_HBARC*M_HBARC);
  const double shear_tensor_normalisation_one_eta=1./(M_HBARC*M_HBARC*M_HBARC*M_HBARC*M_HBARC);
  const double shear_tensor_normalisation_two_eta=1./(M_HBARC*M_HBARC*M_HBARC*M_HBARC*M_HBARC*M_HBARC);

  //
  const double tau_in_fm=TOut_Full->tau;
  const double tau_in_GeV=tau_in_fm/M_HBARC;

  //>> dummy >>dummy >>dummy >>dummy >>neta >> dummy >> nx >> dummy >> ny >> dummy >> deta >> dummy >> dx >> dummy >> dy;
  //# dummy 1 etamax= 1 xmax= 200 ymax= 200 deta= 0 dx= 0.17 dy= 0.17
  std::stringstream header;
  header << "# tau_in_fm " << tau_in_fm << " etamax= 1 xmax= " << Nx << " ymax= " << Ny << " deta= 0 dx= " << dx << " dy= " << dy << "\n";
  outfile << header.str().c_str();

  // Careful, the order of x and y is not the same as in input file
  for(int ix=0;ix<Nx;ix++) {
    double xpos=dx*(ix-Nx/2.0);
    for(int iy=0;iy<Ny;iy++) {
      double ypos=dy*(iy-Ny/2.0);

      // Get Tmunu in a convenient format
      void get_Tmunu_raised_hyperbolic(class EnergyMomentumTensorMap * Tmunu, int ix, int iy, double Tmunu_raised[4][4]);
      double Tmunu_Full_raised[4][4];
      get_Tmunu_raised_hyperbolic(TOut_Full, ix, iy, Tmunu_Full_raised);
      double Tmunu_BG_raised[4][4];
      get_Tmunu_raised_hyperbolic(TOut_BG, ix, iy, Tmunu_BG_raised);


      // Get the energy density and flow
      //void get_energy_flow_linear_improved(double Tmunu[4][4], double Tmunu_bg[4][4], double tau_in_GeV, double &e, double flow[]);
      //double energy_density;
      //double flow[4];
      //get_energy_flow_linear_improved(Tmunu_Full_raised,Tmunu_BG_raised,tau_in_GeV,energy_density,flow);

      void get_energy_flow_from_Tmunu(double Tmunu_raised[4][4], double tau_in_GeV, double * energy_density, double flow[], enum error_codes * error_code);
      double energy_density;
      double flow[4];
      enum error_codes error;
      get_energy_flow_from_Tmunu(Tmunu_Full_raised,tau_in_GeV,&energy_density,flow,&error);
      if (error != decomp_ok) {
        std::cout << "No solution to eigenvalue problem. Using \"linear improved\" energy and flow extraction. Energy density=" << energy_density/(M_HBARC*M_HBARC*M_HBARC) << "GeV/fm^3\n";
        TOut_BG->SetCellData(15, ix, iy, 1);
      }
      else {
        TOut_BG->SetCellData(15, ix, iy, 0);
      }


      // Get the shear stress tensor and the bulk pressure
      void get_shear_and_bulk_from_Tmunu(double Tmunu_raised[4][4], double tau_in_GeV, double energy_density, double * flow, double pimunu_raised[4][4], double & bulk_pressure);
      double shear_tensor_raised[4][4];
      double bulk_pressure;
      get_shear_and_bulk_from_Tmunu(Tmunu_Full_raised, tau_in_GeV,energy_density, flow, shear_tensor_raised, bulk_pressure);
      if (use_sigmamunu_NavierStokes) {
        void get_sigma_proper(class EnergyMomentumTensorMap *TOut_Full, class EnergyMomentumTensorMap *TOut_BG, const int ix, const int iy, double sigma[4][4]);
        double get_conformal_entropy(double energy_density);
        const double entropy_density=get_conformal_entropy(energy_density);
        get_sigma_proper(TOut_Full, TOut_BG, ix, iy, shear_tensor_raised);
        for(int mu=0;mu<4;mu++) for(int nu=0;nu<4;nu++) shear_tensor_raised[mu][nu]*=2*entropy_density*KoMPoSTParameters::EtaOverS;
      }

      // Output information

      //profile >> dummy1 >> dummy2 >> dummy3 >> density >> utau >> ux >> uy  >> dummy  >> dummy  >> dummy  >> dummy;
      //energy density should be in GeV/fm^3
      //utau should be in fm^-1 (T^{\tau/x/y} in fm^-4, T^{\tau/x/y \eta} in fm^-5, T^{\eta\eta} in fm^-6)
      outfile << "0 " << xpos << " " << ypos << " " << energy_density/(M_HBARC*M_HBARC*M_HBARC) << " " << flow[0] << " " << flow[1] << " " << flow[2] << " " << flow[3]/M_HBARC << " ";
      //\pi^{\mu\nu} should be in fermi
      //order: pitautau >> pitaux >> pitauy >> pitaueta >> pixx >> pixy >> pixeta >> piyy >> piyeta >> pietaeta
      outfile << shear_tensor_raised[0][0]*shear_tensor_normalisation << " " << shear_tensor_raised[0][1]*shear_tensor_normalisation << " " << shear_tensor_raised[0][2]*shear_tensor_normalisation << " " << shear_tensor_raised[0][3]*shear_tensor_normalisation_one_eta << " " << shear_tensor_raised[1][1]*shear_tensor_normalisation << " "<< shear_tensor_raised[1][2]*shear_tensor_normalisation << " " << shear_tensor_raised[1][3]*shear_tensor_normalisation_one_eta << " " << shear_tensor_raised[2][2]*shear_tensor_normalisation << " " << shear_tensor_raised[2][3]*shear_tensor_normalisation_one_eta << " " << shear_tensor_raised[3][3]*shear_tensor_normalisation_two_eta << " ";
      outfile << "\n";
    }
  }

  outfile.close();

}

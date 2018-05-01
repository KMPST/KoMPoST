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
#ifndef EnergyMomentumTensorMap_h
#define EnergyMomentumTensorMap_h

#include "EventInput.h"

// ENERGYMOMENTUMTESORMAP -- data holder of input/output energy-momentum tensor
class EnergyMomentumTensorMap {
    
public:

    // Number of data cells allocated per cell
    static const int NCellData = 16;
    // Number of fluid velocity components
    static const int NUi=3 ;
    
    //////////////////////////////////////////
    // COMPONENTS OF ENERGY MOMENTUM TENSOR //
    //////////////////////////////////////////
    
    // T00= T^{tau tau}
    // TXX= T^{xx}
    // TYY= T^{yy}
    // TZZ= 1/t^2 T_{zz} = t^2 T^{zz}
    
    // T0X=    (F_{ty}F^{y}_{x}+F_{tz}F^{z}_{x})  =      T_{tx} = -T^{tx}
    // T0Y=    (F_{tz}F^{z}_{y}+F_{tx}F^{x}_{y})  =      T_{ty} = -T^{ty}
    // T0Z= 1/t(F_{tx}F^{x}_{z}+F_{ty}F^{y}_{z})  =  1/t T_{tz} = -t T^{tz}
    
    // TXY=    -(F_{xt}F^{t}_{y}+F_{xz}F^{z}_{y}) =     -T_{xy} = -T^{xy}
    // TYZ= -1/t(F_{xt}F^{t}_{z}+F_{xy}F^{y}_{z}) = -1/t T_{yz} = -t T^{yz}
    // TZX= -1/t(F_{zt}F^{t}_{x}+F_{zy}F^{y}_{x}) = -1/t T_{zx} = -t T^{zx}
    
public:
    
    // EVOLUTION TIME in fm//
    double tau;
    // ENERGY MOMENTUM TENSOR  -- T_{mu nu} in GeV**4 //
    double *T;
    // Additional data about the corresponding fluid state //
    double *CellData;
    // Energy density from reconstructed stress tensor
    double *Ed ;
    // Flow velocity  from reconstructed stress tensor
    double *Ui ;
    
    // INDEXING //
    int Index2D(int x,int y){
        using EventInput::Ns;
        return x+Ns*y;
    }
    // Index of mu nu component of stress tensor
    int Index(int mu,int nu,int xS,int yS){
        return mu+4*(nu+4*(Index2D(xS,yS)));
    }
    // Index of the id data slot of the CellData
    int IndexCellData(int id,int xS,int yS){
        return id + NCellData*(Index2D(xS,yS));
    }
    // Index of the id data slot of the CellData
    int IndexUi(int i,int xS,int yS){
        return i + NUi*(Index2D(xS,yS));
    }
    
    // GET VALUES //
    double Get(int mu,int nu,int xS,int yS){
        return T[Index(mu,nu,xS,yS)];
    }
    // Get Data //
    double GetCellData(int id,int xS,int yS){
        return CellData[IndexCellData(id,xS,yS)];
    }
    // Get Energy Density //
    double GetEd(int xS,int yS){
        return Ed[Index2D(xS,yS)];
    }
    // Get Flow Velocity //
    double GetUi(int i,int xS,int yS){
        return Ui[IndexUi(i,xS,yS)];
    }

    // SET COMPONENT VALUES //
    void SetComponent(int mu,int nu,int xS,int yS,double Value){
        T[Index(mu,nu,xS,yS)]=Value;
    }
    // SET CELLDATA //
    void SetCellData(int id,int xS,int yS,double Value){
        CellData[IndexCellData(id,xS,yS)]=Value;
    }
    // SET Ed //
    void SetEd(int xS,int yS,double Value){
        Ed[Index2D(xS,yS)]=Value;
    }
    // SET Ui //
    void SetUi(int i,int xS,int yS,double Value){
        Ui[IndexUi(i,xS,yS)]=Value;
    }
    // Set the entire stress tensor 
    void Set(int xS,int yS,double T00,double TXX,double TYY,double TZZ,double T0X,double T0Y,double T0Z,double TXY,double TYZ,double TZX){
        T[Index(0,0,xS,yS)]=T00;    T[Index(0,1,xS,yS)]=T0X;    T[Index(0,2,xS,yS)]=T0Y;    T[Index(0,3,xS,yS)]=T0Z;
        T[Index(1,0,xS,yS)]=T0X;    T[Index(1,1,xS,yS)]=TXX;    T[Index(1,2,xS,yS)]=TXY;    T[Index(1,3,xS,yS)]=TZX;
        T[Index(2,0,xS,yS)]=T0Y;    T[Index(2,1,xS,yS)]=TXY;    T[Index(2,2,xS,yS)]=TYY;    T[Index(2,3,xS,yS)]=TYZ;
        T[Index(3,0,xS,yS)]=T0Z;    T[Index(3,1,xS,yS)]=TZX;    T[Index(3,2,xS,yS)]=TYZ;    T[Index(3,3,xS,yS)]=TZZ;
    }

    // Set the stress tensor from a one dimensional array, consisting of the
    // raised components of T^{munu} stored in order
    void SetRaised(int xS, int yS, double *tmunu_raised)  {
        double u = tau / M_HBARC ;
        double measure[16] = {1., -1., -1., -u ,
                      -1.,  1., -1., -u ,
                      -1., -1.,  1., -u ,
                      -u , -u , -u , u*u } ;
        for (int id = 0 ; id < 16 ; id++) {
            T[id + 16*Index2D(xS,yS) ] = tmunu_raised[id]*measure[id] ;
        }
    }

    // Get the components of the stress tensor as a one dimnensional array
    void GetRaised(int xS, int yS, double *tmunu_raised) {
        double u = M_HBARC / tau ;
        double measure[16] = {1., -1., -1., -u ,
                      -1.,  1., -1., -u ,
                      -1., -1.,  1., -u ,
                      -u , -u , -u , u*u } ;
        for (int id = 0 ; id < 16 ; id++) {
            tmunu_raised[id] = T[id + 16*Index2D(xS,yS)] * measure[id] ;
        }
    }

    
    // RESET //
    void Reset(){
        using EventInput::Ns;
        // SET ALL ENTRIES TO ZERO //
        for(int yS=0;yS<Ns;yS++){
            for(int xS=0;xS<Ns;xS++){
                for(int nu=0;nu<4;nu++){
                    for(int mu=0;mu<4;mu++){
                        T[Index(mu,nu,xS,yS)]=0.0;
                    }
                }
                for(int id=0;id<NCellData;id++){
                    CellData[IndexCellData(id,xS,yS)]=0.0;
                }
                Ed[Index2D(xS,yS)]=0.0;
                for(int k=0;k<NUi;k++){
                    Ui[IndexUi(k,xS,yS)]=0.0;
                }
            }
        }
    }
    // CONSTRUCTOR //
    EnergyMomentumTensorMap(double tau){
        
        using EventInput::Ns;

        // SET TIME //
        this->tau=tau;
        
        // ALLOCATE MEMORY for T//
        this->T=new double[4*4*Ns*Ns];

        // ALLOCATE MEMORY for Additional CellData //
        this->CellData = new double[NCellData*Ns*Ns];
        
        // ALLOCATE MEMORY for Energy Density //
        this->Ed = new double[Ns*Ns];

        // Allocate memory for velocity
        this->Ui = new double[NUi*Ns*Ns];
        
        // RESET //
        this->Reset();
    }
    
    // DE-STRUCTOR //
    ~EnergyMomentumTensorMap(){
        delete T;
        delete CellData;
        delete Ed;
        delete Ui;
    }
};
#endif

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
#include <iostream>
#include <fstream>
#include <cmath>
#include "EventInput.h"
#include "BackgroundEvolution.h"
#include "ScalingVariable.h"

#ifndef M_HBARC
#define M_HBARC  0.197326979
#endif


namespace BackgroundEvolution {


    /////////////////////////////////////////
    // KINETIC THEORY BACKGROUND EVOLUTION //
    /////////////////////////////////////////

    namespace KineticTheory{
        //Number of DOF
        const int NuG=16;

        // MAXIMUM AND MINIMUM SCALING VARIABLES //
        double sMin=0.0; double sMax=512.0;


        double EnergyScalingCurveValue(double x){


            double C2,F0,b,c,xSwitch,xRange;

            // DEFINE YANG-MILLS CONSTANTS //
            if(KoMPoSTParameters::KineticTheory=="EKT"){

                C2=1.02;

                F0 =0.3037314623;
                b = -0.02260000191;
                c = 0.00084856491;

                xSwitch=8.16814089933;
                xRange=3.14159265359;

            }

            // DEFINE RTA CONSTANTS //
            if(KoMPoSTParameters::KineticTheory=="RTA"){

                C2=1.17265;

                F0=0.32862;
                b=-0.0296451;
                c=0.000972941;

                xSwitch=10.5228;
                xRange=6.00391;

            }


            // KoMPoST PARAMETRIZATION //
            double eHD=1.0/(1.0+8./3./x+8./9.*(5+C2)/x/x);
            double eFS=sqrt(tanh(x*pow(F0*(1+b*x+c*x*x),2)));

            double Switch=0.5*(1.0+tanh((x*x-xSwitch*xSwitch)/x/xRange));

            return Switch*eHD+(1.0-Switch)*eFS;

        }

        double PressureScalingCurveValue(double x){
            double dx=0.001*x;
            return EnergyScalingCurveValue(x)-2*x*(EnergyScalingCurveValue(x+dx)-EnergyScalingCurveValue(x-dx))/(2*dx);
        }

        // DETERMINE SCALING FACTOR K [GeV] //
        double DetermineScalingFactor(double EIn,double tIn, ScalingVariable &sv, double &ScalingVarIn){

            // CONVERT t [fm/c] TO t [GeV^-1] //
            double tInGeV=tIn/M_HBARC;

            double EtaOverS=sv.GetEtaOverS0();
            double KLow=std::pow(sMin*EtaOverS,3.0/2.0)/(tInGeV);
            double KHigh=std::pow(sMax*EtaOverS,3.0/2.0)/(tInGeV);

            while(std::abs(KHigh-KLow)>0.001*(KHigh+KLow)/2.0){

                // CHECK MID-POINT //
                double KValue=0.5*(KLow+KHigh);
                double sValue=sv.ScalingVar(tIn, KValue) ;

                // COMPUTE EKT ENERGY //
                double ERec=NuG*(M_PI*M_PI)/30.0*(KValue*KValue*KValue*KValue)/(KValue*tInGeV*sqrt(EtaOverS))*EnergyScalingCurveValue(sValue)/sqrt(sValue);

                // CHOOSE NEW INTERVAL //
                if(ERec>EIn){
                    KHigh=KValue;
                }
                else{
                    KLow=KValue;
                }

            }

            // PERFORM LINEAR INTPEROLATION //
            double sLow=sv.ScalingVar(tIn,KLow);
            double ELow=NuG*(M_PI*M_PI)/30.0*(KLow*KLow*KLow*KLow)/(KLow*tInGeV*sqrt(EtaOverS))*EnergyScalingCurveValue(sLow)/sqrt(sLow);

            double sHigh=sv.ScalingVar(tIn,KHigh);
            double EHigh=NuG*(M_PI*M_PI)/30.0*(KHigh*KHigh*KHigh*KHigh)/(KHigh*tInGeV*sqrt(EtaOverS))*EnergyScalingCurveValue(sHigh)/sqrt(sHigh);

            double KValue = KLow + (EIn-ELow)/(EHigh-ELow)*(KHigh-KLow);
            ScalingVarIn = sv.ScalingVar(tIn,KValue);

            return KValue;

        }

        // PERFORM PROPAGATION OF ENERGY MOMENTUM TENSOR //
        void Propagate(double tOut,double KValue,double ScalingVarOut,double EtaOverS,double &T00,double &TXX,double &TYY,double &TZZ){

            // CONVERT t [fm/c] TO t [GeV^-1] //
            double tOutGeV=tOut/M_HBARC;

            // SET ENERGY-MOMENTUM TENSOR VALUES //
            T00=NuG*(M_PI*M_PI)/30.0*(KValue*KValue*KValue*KValue)/(KValue*tOutGeV*sqrt(EtaOverS))*EnergyScalingCurveValue(ScalingVarOut)/sqrt(ScalingVarOut);
            TZZ=NuG*(M_PI*M_PI)/90.0*(KValue*KValue*KValue*KValue)/(KValue*tOutGeV*sqrt(EtaOverS))*PressureScalingCurveValue(ScalingVarOut)/sqrt(ScalingVarOut);
            TXX=0.5*(T00-TZZ);
            TYY=0.5*(T00-TZZ);

        }

        // CHECK MATCHING EFFICIENCY OF T00,TXX,TYY,TZZ //
        void CheckMatchingEfficiency(int xS,int yS,double T00In,double TXXIn,double TYYIn,double TZZIn,double tIn,ScalingVariable &sv){

            // CONVERT t [fm/c] TO t [GeV^-1] //
            double tInGeV=tIn/M_HBARC;

            // DETERMINE SCALING FACTOR AND SCALED OUTPUT TIME //
            double EtaOverS = sv.GetEtaOverS0() ;
            double sInValue ;
            double KValue=BackgroundEvolution::KineticTheory::DetermineScalingFactor(T00In,tIn,sv,sInValue);

            // GET RECONSTRUCTED ENERGY-MOMENTUM TENSOR VALUES //
            double T00Rec=NuG*(M_PI*M_PI)/30.0*(KValue*KValue*KValue*KValue)/(KValue*tInGeV*sqrt(EtaOverS))*EnergyScalingCurveValue(sInValue)/sqrt(sInValue);
            double TZZRec=NuG*(M_PI*M_PI)/90.0*(KValue*KValue*KValue*KValue)/(KValue*tInGeV*sqrt(EtaOverS))*PressureScalingCurveValue(sInValue)/sqrt(sInValue);
            double TXXRec=0.5*(T00Rec-TZZRec);
            double TYYRec=0.5*(T00Rec-TZZRec);

            // CREATE OUTPUT //
            std::cout << xS << " " << yS << " " << T00In << " " << T00Rec << " " << TXXIn << " " << TXXRec << " " << TYYIn << " " << TYYRec << " " << TZZIn << " " << TZZRec << std::endl;

        }


    }


    /////////////////////////////////////////
    // FREE STREAMING BACKGROUND EVOLUTION //
    /////////////////////////////////////////

    namespace FreeStreaming{

        // PERFORM PROPAGATION OF ENERGY MOMENTUM TENSOR //
        void Propagate(double T00In,double TXXIn,double TYYIn,double TZZIn,double tIn,double tOut,double &T00,double &TXX,double &TYY,double &TZZ){

            // SET ENERGY-MOMENTUM TENSOR VALUES //
            T00=T00In*(tIn/tOut); TXX=0.5*T00; TYY=0.5*T00; TZZ=0.0;

        }

        // CHECK MATCHING EFFICIENCY OF T00,TXX,TYY,TZZ //
        void CheckMatchingEfficiency(int xS,int yS,double T00In,double TXXIn,double TYYIn,double TZZIn,double tIn){

            // GET RECONSTRUCTED ENERGY-MOMENTUM TENSOR VALUES //
            double T00Rec=T00In;
            double TZZRec=0.0;
            double TXXRec=0.5*(T00Rec-TZZRec);
            double TYYRec=0.5*(T00Rec-TZZRec);

            // CREATE OUTPUT //
            std::cout << xS << " " << yS << " " << T00In << " " << T00Rec << " " << TXXIn << " " << TXXRec << " " << TYYIn << " " << TYYRec << " " << TZZIn << " " << TZZRec << std::endl;

        }
    }

}

#ifndef TOOLS
#define TOOLS
#include <string>
#include <iostream>
#include "TMath.h"

inline float dPhi(float phi1, float phi2){
  return TMath::ACos(TMath::Cos(phi1-phi2));
}

inline float dR(float eta1, float phi1, float eta2, float phi2){
  return TMath::Power( TMath::Power(eta1-eta2,2) + TMath::Power( dPhi(phi1, phi2), 2), 0.5);
}
#endif

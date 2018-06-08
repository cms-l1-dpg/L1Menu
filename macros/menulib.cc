/* automatically generated from L1Menu_Collisions2018_v1_0_0 with menu2lib.py */
/* https://gitlab.cern.ch/cms-l1t-utm/scripts */

#include <algorithm>
#include <map>
#include <string>
#include <sstream>

#include "menulib.hh"

//
// common functions for algorithm implementations
//
std::pair<double, double>
get_missing_et(L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower,
               const int max_eta,
               const double threshold)
{
  // https://github.com/cms-sw/cmssw/blob/CMSSW_9_0_X/L1Trigger/L1TCalorimeter/src/CaloTools.cc#L13=L15
  const int64_t cos_coeff[72] = {1023, 1019, 1007, 988, 961, 927, 886, 838, 784, 723, 658, 587, 512, 432, 350, 265, 178, 89, 0, -89, -178, -265, -350, -432, -512, -587, -658, -723, -784, -838, -886, -927, -961, -988, -1007, -1019, -1023, -1019, -1007, -988, -961, -927, -886, -838, -784, -723, -658, -587, -512, -432, -350, -265, -178, -89, 0, 89, 178, 265, 350, 432, 511, 587, 658, 723, 784, 838, 886, 927, 961, 988, 1007, 1019};

  const int64_t sin_coeff[72] = {0, 89, 178, 265, 350, 432, 512, 587, 658, 723, 784, 838, 886, 927, 961, 988, 1007, 1019, 1023, 1019, 1007, 988, 961, 927, 886, 838, 784, 723, 658, 587, 512, 432, 350, 265, 178, 89, 0, -89, -178, -265, -350, -432, -512, -587, -658, -723, -784, -838, -886, -927, -961, -988, -1007, -1019, -1023, -1019, -1007, -988, -961, -927, -886, -838, -784, -723, -658, -587, -512, -432, -350, -265, -178, -89};

  if (not calo_tower) return std::make_pair(-1., -9999.);

  double ex = 0.;
  double ey = 0.;

  for (int ii = 0; ii < calo_tower->nTower; ii++)
  {
    if (abs(calo_tower->ieta.at(ii)) <= max_eta)
    {
      const double et = calo_tower->iet.at(ii) * 0.5;
      if (et > threshold)
      {
        const int index = calo_tower->iphi.at(ii) - 1;
        ex += (et*cos_coeff[index]/1024.);
        ey += (et*sin_coeff[index]/1024.);
      }
    }
  }

  return std::make_pair(sqrt(ex*ex + ey*ey), atan2(-ey, -ex));
}


double
get_total_ht(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade,
             const int max_eta,
             const double threshold)
{
  double sum = 0.;

  for (int ii = 0; ii < upgrade->nJets; ii++)
  {
    if (upgrade->jetBx.at(ii) != 0) continue;

    if (abs(upgrade->jetIEta.at(ii)) <= 2*max_eta)
    {
      const double et = upgrade->jetEt.at(ii);
      if (et > threshold)
      {
        sum += et;
      }
    }
  }

  return sum;
} 


double
get_transverse_mass(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade,
                    const double threshold_eg,
                    const double threshold_met)
{
  double mt = -1.;

  if (upgrade->nEGs == 0) return mt;

  // leading-eg
  int id_leading_eg = -1;
  for (int ii = 0; ii < upgrade->nEGs; ii++)
  {
    if (upgrade->egBx.at(ii) != 0) continue;
    if (id_leading_eg < 0)
    {
      id_leading_eg = ii;
      break;
    }
  }

  if (id_leading_eg < 0) return mt;

  const double eg_et = upgrade->egEt.at(id_leading_eg);
  const double eg_phi = upgrade->egPhi.at(id_leading_eg);

  if (eg_et < threshold_eg) return mt;


  // missing-Et
  int id_missing_et = -1;
  for (int ii = 0; ii < upgrade->nSums; ii++)
  {
    if (upgrade->sumBx.at(ii) != 0) continue;
    if (upgrade->sumType.at(ii) == L1Analysis::kMissingEt)
    {
      id_missing_et = ii;
      break;
    }
  }

  if (id_missing_et < 0) return mt;

  const double met_et = upgrade->sumEt.at(id_missing_et);
  const double met_phi = upgrade->sumPhi.at(id_missing_et);

  if (met_et < threshold_met) return mt;


  // mt
  double delta_phi = eg_phi - met_phi;
  while (delta_phi >= M_PI) delta_phi -= 2.*M_PI;
  while (delta_phi < -M_PI) delta_phi += 2.*M_PI;

  mt = sqrt(2.*eg_et*met_et*(1. - cos(delta_phi)));
  return mt;
}


// utility methods
void
getCombination(int N,
               int K,
               std::vector<std::vector<int> >& combination)
{
  std::string bitmask(K, 1);
  bitmask.resize(N, 0);

  do
  {
    std::vector<int> set;
    for (int ii = 0; ii < N; ++ii)
    {
      if (bitmask[ii]) set.push_back(ii);
    }
    combination.push_back(set);
  }
  while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}


void
getPermutation(int N,
               std::vector<std::vector<int> >& permutation)
{
  std::vector<int> indicies(N);
  for (int ii = 0; ii < N; ii++) indicies.at(ii) = ii;

  do
  {
    std::vector<int> set;
    for (int ii = 0; ii < N; ++ii)
    {
      set.push_back(indicies.at(ii));
    }
    permutation.push_back(set);
  }
  while (std::next_permutation(indicies.begin(), indicies.end()));
}




//
// NB: tmEventSetup.XxxWithOverlapRemoval was removed between utm-overlapRemoval-xsd330 and utm_0.6.5
//
// generate conditions
    





bool
CaloCaloCorrelation_1626121879521236767
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 80)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 48));
            
          if (not etaWindow1) continue;
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 48));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.6 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
CaloCaloCorrelation_18379087122140179561
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(ii) >= 44)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
      nobj1++;
              if (nobj1 > 12) break;
              
                                      // TAU70: ET >= 140 at BX = 0
      if (not (data->tauIEt.at(jj) >= 140)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(jj)) and (data->tauIEta.at(jj) <= 48));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_TAU[deltaIEta];
  
    int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_EG_TAU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.09 * POW10[6]);
  maximum = (long long)(139.24 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
CaloCaloCorrelation_3813196582576312175
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET32: ET >= 64 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 64)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 52));
            
          if (not etaWindow1) continue;
                                // JET32: ET >= 64 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 64)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 52));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.6 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
CaloCaloCorrelation_3813196582576378703
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 80)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 52));
            
          if (not etaWindow1) continue;
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 52));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.6 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
CaloCaloCorrelation_7041035331702023693
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 200)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 52));
            
          if (not etaWindow1) continue;
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 200)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 52));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.6 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
CaloCaloCorrelation_7041035331710545453
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET112: ET >= 224 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 224)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 52));
            
          if (not etaWindow1) continue;
                                // JET112: ET >= 224 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 224)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 52));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.6 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
CaloCaloCorrelation_9154320878437213226
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(ii) >= 52)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->jetBx.size(); jj++)
    {
      if (not (data->jetBx.at(jj) == 0)) continue;
      nobj1++;
              
                                      // JET34: ET >= 68 at BX = 0
      if (not (data->jetIEt.at(jj) >= 68)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(jj)) and (data->jetIEta.at(jj) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_JET[deltaIEta];
  
    int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_EG_JET[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.09 * POW10[6]);
  maximum = (long long)(139.24 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
CaloCaloCorrelation_9154320878437278762
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(ii) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->jetBx.size(); jj++)
    {
      if (not (data->jetBx.at(jj) == 0)) continue;
      nobj1++;
              
                                      // JET34: ET >= 68 at BX = 0
      if (not (data->jetIEt.at(jj) >= 68)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(jj)) and (data->jetIEta.at(jj) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_JET[deltaIEta];
  
    int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_EG_JET[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.09 * POW10[6]);
  maximum = (long long)(139.24 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
CaloCaloCorrelation_9154320878441210922
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(ii) >= 60)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->jetBx.size(); jj++)
    {
      if (not (data->jetBx.at(jj) == 0)) continue;
      nobj1++;
              
                                      // JET34: ET >= 68 at BX = 0
      if (not (data->jetIEt.at(jj) >= 68)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(jj)) and (data->jetIEta.at(jj) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_JET[deltaIEta];
  
    int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_EG_JET[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.09 * POW10[6]);
  maximum = (long long)(139.24 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
CaloCaloCorrelation_980305370107485378
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(ii) >= 44)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
      nobj1++;
              if (nobj1 > 12) break;
              
                                      // TAU26: ET >= 52 at BX = 0
      if (not (data->tauIEt.at(jj) >= 52)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(jj)) and (data->tauIEta.at(jj) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(jj)) & 1)) continue;

          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_TAU[deltaIEta];
  
    int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_EG_TAU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.09 * POW10[6]);
  maximum = (long long)(139.24 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
CaloCaloCorrelation_980305438827093186
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(ii) >= 48)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
      nobj1++;
              if (nobj1 > 12) break;
              
                                      // TAU27: ET >= 54 at BX = 0
      if (not (data->tauIEt.at(jj) >= 54)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(jj)) and (data->tauIEta.at(jj) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(jj)) & 1)) continue;

          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_TAU[deltaIEta];
  
    int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_EG_TAU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.09 * POW10[6]);
  maximum = (long long)(139.24 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      






bool
CaloMuonCorrelation_12420459208926927530
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET32: ET >= 64 at BX = 0
      if (not (data->jetIEt.at(ii) >= 64)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 52));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU10: ET >= 21 at BX = 0
      if (not (data->muonIEt.at(jj) >= 21)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000625
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(jj)) and (data->muonIEtaAtVtx.at(jj) <= 211));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_14690273421121723050
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(ii) >= 80)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 52));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(jj) >= 25)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000625
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(jj)) and (data->muonIEtaAtVtx.at(jj) <= 211));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_15675017008716733697
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(ii) >= 70)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_17980860017930427617
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET16: ET >= 32 at BX = 0
      if (not (data->jetIEt.at(ii) >= 32)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_2018052583404954969
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(ii) >= 240)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_2018052583404955481
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(ii) >= 240)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.80
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.641 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_4145801962648263985
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(ii) >= 120)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_4145801962648264017
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(ii) >= 160)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_8802031396140112104
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(ii) >= 180)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(jj) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.80
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.641 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_8802031396140113640
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(ii) >= 180)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.80
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.641 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_9314160007219977660
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(ii) >= 80)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 48));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(jj) >= 25)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000625
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(jj)) and (data->muonIEtaAtVtx.at(jj) <= 211));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_13299746526184647851
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_2355036583129339571
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_2931778810409473715
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_56167094600470587
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG20: ET >= 40 at BX = 0
      if (not (data->egIEt.at(idx) >= 40)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_56237463344648251
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_56237497704386619
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_56343050820653115
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG25: ET >= 50 at BX = 0
      if (not (data->egIEt.at(idx) >= 50)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_7286147382174463995
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_7286147403649300475
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG15: ET >= 30 at BX = 0
      if (not (data->egIEt.at(idx) >= 30)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_7286147425124136955
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_7286147489548646395
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG15: ET >= 30 at BX = 0
      if (not (data->egIEt.at(idx) >= 30)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG15: ET >= 30 at BX = 0
      if (not (data->egIEt.at(idx) >= 30)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_7286147532498319355
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG17: ET >= 34 at BX = 0
      if (not (data->egIEt.at(idx) >= 34)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG17: ET >= 34 at BX = 0
      if (not (data->egIEt.at(idx) >= 34)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_7286147931930277883
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG20: ET >= 40 at BX = 0
      if (not (data->egIEt.at(idx) >= 40)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_7286147940520212475
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_7286147987764852731
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG25: ET >= 50 at BX = 0
      if (not (data->egIEt.at(idx) >= 50)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_7286148022124591099
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG25: ET >= 50 at BX = 0
      if (not (data->egIEt.at(idx) >= 50)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG14: ET >= 28 at BX = 0
      if (not (data->egIEt.at(idx) >= 28)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_7286148030714525691
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG27: ET >= 54 at BX = 0
      if (not (data->egIEt.at(idx) >= 54)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG14: ET >= 28 at BX = 0
      if (not (data->egIEt.at(idx) >= 28)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_14043965433308858268
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_16307690244847013269
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx) >= 200)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx) >= 60)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_16379747838884941845
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET110: ET >= 220 at BX = 0
      if (not (data->jetIEt.at(idx) >= 220)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_16382562588652064149
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET115: ET >= 230 at BX = 0
      if (not (data->jetIEt.at(idx) >= 230)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_16451805432922886165
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx) >= 90)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_17548339888472803228
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET75: ET >= 150 at BX = 0
      if (not (data->jetIEt.at(idx) >= 150)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET65: ET >= 130 at BX = 0
      if (not (data->jetIEt.at(idx) >= 130)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_209751802956826525
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx) >= 160)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET70: ET >= 140 at BX = 0
      if (not (data->jetIEt.at(idx) >= 140)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_254798794346809245
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET85: ET >= 170 at BX = 0
      if (not (data->jetIEt.at(idx) >= 170)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET75: ET >= 150 at BX = 0
      if (not (data->jetIEt.at(idx) >= 150)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_4162612533456677351
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(idx) >= 180)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx) >= 60)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_7222390773192019523
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx) >= 200)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx) >= 200)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_7222953723145440851
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_7223798148075572843
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET150: ET >= 300 at BX = 0
      if (not (data->jetIEt.at(idx) >= 300)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET150: ET >= 300 at BX = 0
      if (not (data->jetIEt.at(idx) >= 300)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_12394181525097886766
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU10: ET >= 21 at BX = 0
      if (not (data->muonIEt.at(idx) >= 21)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU10: ET >= 21 at BX = 0
      if (not (data->muonIEt.at(idx) >= 21)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_13627348644483379947
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064375
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 184));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064375
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 184));
            
          if (not etaWindow1) continue;
            
      // charge correlation
      bool equal = true;
      bool invalid = false;
      for (size_t mm = 0; mm < 2 -1; mm++)
      {
        int idx0 = candidates.at(set.at(indicies.at(mm)));
        int idx1 = candidates.at(set.at(indicies.at(mm+1)));
        if ((data->muonChg.at(idx0) == 0) or (data->muonChg.at(idx1) == 0))
        {
          invalid = true;
          break;
        }
        if (data->muonChg.at(idx0) != data->muonChg.at(idx1))
        {
          equal = false;
          break;
        }
      }
      if (invalid) continue;

      // charge correlation: "os"
      if (equal) continue;

      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_14585777620730815295
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_14585778097856672575
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_14585778268989477695
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_14585786515326686015
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_14585789264105755455
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx) >= 9)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx) >= 9)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_14585792012884824895
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_14585803008001102655
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU9: ET >= 19 at BX = 0
      if (not (data->muonIEt.at(idx) >= 19)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU9: ET >= 19 at BX = 0
      if (not (data->muonIEt.at(idx) >= 19)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_14617142003772573591
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      // charge correlation
      bool equal = true;
      bool invalid = false;
      for (size_t mm = 0; mm < 2 -1; mm++)
      {
        int idx0 = candidates.at(set.at(indicies.at(mm)));
        int idx1 = candidates.at(set.at(indicies.at(mm+1)));
        if ((data->muonChg.at(idx0) == 0) or (data->muonChg.at(idx1) == 0))
        {
          invalid = true;
          break;
        }
        if (data->muonChg.at(idx0) != data->muonChg.at(idx1))
        {
          equal = false;
          break;
        }
      }
      if (invalid) continue;

      // charge correlation: "os"
      if (equal) continue;

      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_16323903523977050720
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU18: ET >= 37 at BX = 0
      if (not (data->muonIEt.at(idx) >= 37)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU18: ET >= 37 at BX = 0
      if (not (data->muonIEt.at(idx) >= 37)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_16961157256621881348
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(idx) >= 25)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_16961159554147985412
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU15: ET >= 31 at BX = 0
      if (not (data->muonIEt.at(idx) >= 31)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_16961163303935834116
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU15: ET >= 31 at BX = 0
      if (not (data->muonIEt.at(idx) >= 31)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx) >= 15)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_16961163952194496516
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU15: ET >= 31 at BX = 0
      if (not (data->muonIEt.at(idx) >= 31)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx) >= 15)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_2011765979326275391
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      // charge correlation
      bool equal = true;
      bool invalid = false;
      for (size_t mm = 0; mm < 2 -1; mm++)
      {
        int idx0 = candidates.at(set.at(indicies.at(mm)));
        int idx1 = candidates.at(set.at(indicies.at(mm+1)));
        if ((data->muonChg.at(idx0) == 0) or (data->muonChg.at(idx1) == 0))
        {
          invalid = true;
          break;
        }
        if (data->muonChg.at(idx0) != data->muonChg.at(idx1))
        {
          equal = false;
          break;
        }
      }
      if (invalid) continue;

      // charge correlation: "os"
      if (equal) continue;

      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_3139255731352238604
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      // charge correlation
      bool equal = true;
      bool invalid = false;
      for (size_t mm = 0; mm < 2 -1; mm++)
      {
        int idx0 = candidates.at(set.at(indicies.at(mm)));
        int idx1 = candidates.at(set.at(indicies.at(mm+1)));
        if ((data->muonChg.at(idx0) == 0) or (data->muonChg.at(idx1) == 0))
        {
          invalid = true;
          break;
        }
        if (data->muonChg.at(idx0) != data->muonChg.at(idx1))
        {
          equal = false;
          break;
        }
      }
      if (invalid) continue;

      // charge correlation: "os"
      if (equal) continue;

      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_3229327723899648524
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx) >= 9)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx) >= 9)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      // charge correlation
      bool equal = true;
      bool invalid = false;
      for (size_t mm = 0; mm < 2 -1; mm++)
      {
        int idx0 = candidates.at(set.at(indicies.at(mm)));
        int idx1 = candidates.at(set.at(indicies.at(mm+1)));
        if ((data->muonChg.at(idx0) == 0) or (data->muonChg.at(idx1) == 0))
        {
          invalid = true;
          break;
        }
        if (data->muonChg.at(idx0) != data->muonChg.at(idx1))
        {
          equal = false;
          break;
        }
      }
      if (invalid) continue;

      // charge correlation: "os"
      if (equal) continue;

      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_3947425258490794706
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_14812192849374407856
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU32: ET >= 64 at BX = 0
      if (not (data->tauIEt.at(idx) >= 64)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU32: ET >= 64 at BX = 0
      if (not (data->tauIEt.at(idx) >= 64)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_17539608616528615651
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU70: ET >= 140 at BX = 0
      if (not (data->tauIEt.at(idx) >= 140)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU70: ET >= 140 at BX = 0
      if (not (data->tauIEt.at(idx) >= 140)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_5588820814667115697
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU36: ET >= 72 at BX = 0
      if (not (data->tauIEt.at(idx) >= 72)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU36: ET >= 72 at BX = 0
      if (not (data->tauIEt.at(idx) >= 72)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_977134795165985969
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU34: ET >= 68 at BX = 0
      if (not (data->tauIEt.at(idx) >= 68)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU34: ET >= 68 at BX = 0
      if (not (data->tauIEt.at(idx) >= 68)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

                          

  



bool
InvariantMassOvRm_10944688471548334117
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
        // remove overlap -- reference: TAU40
  std::vector<int> reference;
    size_t nref = 0;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nref++;
          if (nref > 12) break;
          
                              // TAU40: ET >= 80 at BX = 0
      if (not (data->tauIEt.at(ii) >= 80)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(ii)) & 1)) continue;

          
    reference.push_back(ii);
  }
  if (not reference.size()) return false;

    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                     int iEta = -9999999; unsigned int deltaIEta = 9999999;
                          
  // remove overlap -- target: JET80
       // 0.00 <= DeltaR <= 0.20
  long long minDeltaR2 = std::numeric_limits<long long>::max();
  const long long cutDeltaR2Min = (long long)(0.0 * POW10[6]);
  const long long cutDeltaR2Max = (long long)(0.041 * POW10[6]);
         
  // compute minimum distance to reference objects
  for (size_t _jj = 0; _jj < reference.size(); _jj++)
  {
    const int index = reference.at(_jj);
                iEta = data->jetIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(index));
      unsigned int deltaEta = LUT_DETA_JET_TAU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_TAU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
      if (deltaR2 < minDeltaR2) minDeltaR2 = deltaR2;
                  }

  // skip if needed
      if ((cutDeltaR2Min <= minDeltaR2) and (minDeltaR2 <= cutDeltaR2Max)) continue;

        
        candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 160)) continue;

          
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
            // 420.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(88200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                          

  



bool
InvariantMassOvRm_10967205787862279205
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
        // remove overlap -- reference: TAU45
  std::vector<int> reference;
    size_t nref = 0;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nref++;
          if (nref > 12) break;
          
                              // TAU45: ET >= 90 at BX = 0
      if (not (data->tauIEt.at(ii) >= 90)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(ii)) & 1)) continue;

          
    reference.push_back(ii);
  }
  if (not reference.size()) return false;

    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                     int iEta = -9999999; unsigned int deltaIEta = 9999999;
                          
  // remove overlap -- target: JET35
       // 0.00 <= DeltaR <= 0.20
  long long minDeltaR2 = std::numeric_limits<long long>::max();
  const long long cutDeltaR2Min = (long long)(0.0 * POW10[6]);
  const long long cutDeltaR2Max = (long long)(0.041 * POW10[6]);
         
  // compute minimum distance to reference objects
  for (size_t _jj = 0; _jj < reference.size(); _jj++)
  {
    const int index = reference.at(_jj);
                iEta = data->jetIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(index));
      unsigned int deltaEta = LUT_DETA_JET_TAU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_TAU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
      if (deltaR2 < minDeltaR2) minDeltaR2 = deltaR2;
                  }

  // skip if needed
      if ((cutDeltaR2Min <= minDeltaR2) and (minDeltaR2 <= cutDeltaR2Max)) continue;

        
        candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 70)) continue;

          
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 70)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
            // 450.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(101250.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_12340992865530516744
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU15: ET >= 31 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 31)) continue;

          
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 15)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 1.0 <= mass <= 151982.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.5 * POW10[6]);
  maximum = (long long)(11549264162.0 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_12401534428203007157
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1, etaWindow2;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 90)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 90)) continue;

                        // -2.697 <= eta <= 2.697
              etaWindow1 = ((-62 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 61));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_13689376201502793133
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000625
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 211));
            
          if (not etaWindow1) continue;
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000625
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 211));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 8.0 <= mass <= 14.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(32.0 * POW10[6]);
  maximum = (long long)(98.0 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_14031536345113029771
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 120)) continue;

          
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -2.697 <= eta <= 2.697
              etaWindow1 = ((-62 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 61));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_14031539093892099211
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 120)) continue;

          
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 90)) continue;

                        // -2.697 <= eta <= 2.697
              etaWindow1 = ((-62 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 61));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_15191958030943548804
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU2p5: ET >= 6 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 6)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 5.0 <= mass <= 17.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(12.5 * POW10[6]);
  maximum = (long long)(144.5 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_15192153509407276420
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 11)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU2p5: ET >= 6 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 6)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 5.0 <= mass <= 17.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(12.5 * POW10[6]);
  maximum = (long long)(144.5 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_16026464453068411612
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1, etaWindow2;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 80)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_16026640374949827292
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1, etaWindow2;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 90)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 90)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_16981538589298500419
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // EG7p5: ET >= 15 at BX = 0
      if (not (data->egIEt.at(idx0) >= 15)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx0)) and (data->egIEta.at(idx0) <= 48));
            
          if (not etaWindow1) continue;
                                // EG7p5: ET >= 15 at BX = 0
      if (not (data->egIEt.at(idx1) >= 15)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx1)) and (data->egIEta.at(idx1) <= 48));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 20.0
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
  
    int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_EG_EG[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_EG_EG[deltaIPhi];
  long long pt0 = LUT_EG_ET[data->egIEt.at(idx0)];
  long long pt1 = LUT_EG_ET[data->egIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.0 * POW10[5]);
  maximum = (long long)(200.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_2342552854377181621
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064375
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 184));
            
          if (not etaWindow1) continue;
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064375
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 184));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 7.0 <= mass <= 18.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(24.5 * POW10[6]);
  maximum = (long long)(162.0 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_2443380592745462540
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // EG3: ET >= 6 at BX = 0
      if (not (data->egIEt.at(idx0) >= 6)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx0)) and (data->egIEta.at(idx0) <= 48));
            
          if (not etaWindow1) continue;
                                // EG3: ET >= 6 at BX = 0
      if (not (data->egIEt.at(idx1) >= 6)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx1)) and (data->egIEta.at(idx1) <= 48));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 20.0
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
  
    int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_EG_EG[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_EG_EG[deltaIPhi];
  long long pt0 = LUT_EG_ET[data->egIEt.at(idx0)];
  long long pt1 = LUT_EG_ET[data->egIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.0 * POW10[5]);
  maximum = (long long)(200.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_2940638391876117895
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

          
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_2940638392207467911
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 160)) continue;

          
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 420.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(88200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_2940649386995017095
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 70)) continue;

          
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 70)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_2940919866919937415
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 80)) continue;

          
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_2940930862038836615
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 90)) continue;

          
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 90)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_2941482817007576455
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 120)) continue;

          
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 120)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_3063833799189854821
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 0.0 <= mass <= 9.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(40.5 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_3424918508618058399
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 57));
            
          if (not etaWindow1) continue;
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.5 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

          // 150.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(11250.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_3425188988478491295
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 57));
            
          if (not etaWindow1) continue;
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.5 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

          // 200.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(20000.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_3425199983594769055
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 57));
            
          if (not etaWindow1) continue;
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 250.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(31250.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

          // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.5 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_3425470463455201951
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 57));
            
          if (not etaWindow1) continue;
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.5 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

          // 300.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(45000.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_3425477060524968607
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 57));
            
          if (not etaWindow1) continue;
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 330.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(54450.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

          // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.5 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_3425483657594735263
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 57));
            
          if (not etaWindow1) continue;
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 360.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(64800.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

          // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.5 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_3915066037994721069
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1, etaWindow2;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 120)) continue;

          
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_3915066038015692589
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1, etaWindow2;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 120)) continue;

          
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 90)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_4461482972834602413
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000625
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 211));
            
          if (not etaWindow1) continue;
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000625
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 211));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 14.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(98.0 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

          const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_7551966572230413519
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

          
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 1.0 <= mass <= 151982.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.5 * POW10[6]);
  maximum = (long long)(11549264162.0 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_7789845660996549812
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1, etaWindow2;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 80)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -2.697 <= eta <= 2.697
              etaWindow1 = ((-62 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 61));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_9620562409477107818
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 80)) continue;

                        // -2.697 <= eta <= 2.697
              etaWindow1 = ((-62 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 61));
            
          if (not etaWindow1) continue;
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -2.697 <= eta <= 2.697
              etaWindow1 = ((-62 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 61));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_9620565158256341098
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 90)) continue;

                        // -2.697 <= eta <= 2.697
              etaWindow1 = ((-62 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 61));
            
          if (not etaWindow1) continue;
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 90)) continue;

                        // -2.697 <= eta <= 2.697
              etaWindow1 = ((-62 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 61));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
      





bool
MuonMuonCorrelation_12923126501326425857
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 0.00 <= DeltaR <= 1.20
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.441 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_15199048927445899759
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 138));
            
          if (not etaWindow1) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 138));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.961 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_15199048929593776303
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064375
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 184));
            
          if (not etaWindow1) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064375
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 184));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.961 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_15498326754036298551
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 1.60
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(2.561 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_16784489743460462578
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 9)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 9)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 0.00 <= DeltaR <= 1.20
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.441 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_5013507948943010765
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == -1)) continue;
    nobj0++;
      
        const int idx0 = ii;
    bool etaWindow1;bool phiWindow1;
                              // MU3-1: ET >= 7 at BX = -1
      if (not (data->muonIEt.at(idx0) >= 7)) continue;

                        // -1.2016875 <= eta <= 1.2016875
              etaWindow1 = ((-110 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 110));
            
                        // 0.523598775598 <= phi <= 2.61799387799
              phiWindow1 = ((48 <= data->muonIPhiAtVtx.at(idx0)) and (data->muonIPhiAtVtx.at(idx0) <= 239));
            
                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

          if (not etaWindow1) continue;if (not phiWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
            const int idx1 = jj;
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 7)) continue;

                        // -1.2016875 <= eta <= 1.2016875
              etaWindow1 = ((-110 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 110));
            
                        // 3.66519142919 <= phi <= 5.75958653158
              phiWindow1 = ((336 <= data->muonIPhiAtVtx.at(idx1)) and (data->muonIPhiAtVtx.at(idx1) <= 527));
            
                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

          if (not etaWindow1) continue;if (not phiWindow1) continue;
          long long minimum;
  long long maximum;

  
        // 2.618 <= DeltaPhi <= 3.142
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
    minimum = (long long)(2.618 * POW10[3]);
  maximum = (long long)(3.142 * POW10[3]);
  if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
MuonMuonCorrelation_5698493964878099256
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 1.60
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(2.561 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_9513481109949270451
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -1.4083125 <= eta <= 1.4083125
              etaWindow1 = ((-129 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 129));
            
          if (not etaWindow1) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -1.4083125 <= eta <= 1.4083125
              etaWindow1 = ((-129 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 129));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.961 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_9513481109957663155
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 138));
            
          if (not etaWindow1) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 138));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.961 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_9513481247421761971
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064375
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 184));
            
          if (not etaWindow1) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064375
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 184));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.961 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      


bool
QuadJET_15326169532950703320
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET36: ET >= 72 at BX = 0
      if (not (data->jetIEt.at(idx) >= 72)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET36: ET >= 72 at BX = 0
      if (not (data->jetIEt.at(idx) >= 72)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET36: ET >= 72 at BX = 0
      if (not (data->jetIEt.at(idx) >= 72)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET36: ET >= 72 at BX = 0
      if (not (data->jetIEt.at(idx) >= 72)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_1795783097020010207
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx) >= 160)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx) >= 90)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 52));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 52));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_1795850802884464351
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx) >= 160)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET50: ET >= 100 at BX = 0
      if (not (data->jetIEt.at(idx) >= 100)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 52));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx) >= 90)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 52));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_284978008326942669
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_2969440952489109684
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET70: ET >= 140 at BX = 0
      if (not (data->jetIEt.at(idx) >= 140)) continue;

                        // -2.3925 <= eta <= 2.3925
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET55: ET >= 110 at BX = 0
      if (not (data->jetIEt.at(idx) >= 110)) continue;

                        // -2.3925 <= eta <= 2.3925
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -2.3925 <= eta <= 2.3925
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // -2.3925 <= eta <= 2.3925
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_2969443065613019316
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET70: ET >= 140 at BX = 0
      if (not (data->jetIEt.at(idx) >= 140)) continue;

                        // -2.3925 <= eta <= 2.3925
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET55: ET >= 110 at BX = 0
      if (not (data->jetIEt.at(idx) >= 110)) continue;

                        // -2.3925 <= eta <= 2.3925
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -2.3925 <= eta <= 2.3925
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -2.3925 <= eta <= 2.3925
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_6278535506187355856
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET95: ET >= 190 at BX = 0
      if (not (data->jetIEt.at(idx) >= 190)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET75: ET >= 150 at BX = 0
      if (not (data->jetIEt.at(idx) >= 150)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET65: ET >= 130 at BX = 0
      if (not (data->jetIEt.at(idx) >= 130)) continue;

          
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET20: ET >= 40 at BX = 0
      if (not (data->jetIEt.at(idx) >= 40)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadMU_509409160461874775
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(3)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadMU_509409667408098135
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(3)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadMU_509409849236703575
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(3)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_10104274574069012747
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

                        // -1.5225 <= eta <= 1.5225
              etaWindow1 = ((-35 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 34));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_10104274582658947339
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // -1.5225 <= eta <= 1.5225
              etaWindow1 = ((-35 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 34));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_10104274591248881931
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -1.5225 <= eta <= 1.5225
              etaWindow1 = ((-35 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 34));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_10104274591248882187
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -1.5225 <= eta <= 1.5225
              etaWindow1 = ((-35 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 34));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_10104275106644957707
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(idx) >= 60)) continue;

                        // -1.5225 <= eta <= 1.5225
              etaWindow1 = ((-35 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 34));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_14243075932402617139
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_14243075932536834867
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG9: ET >= 18 at BX = 0
      if (not (data->egIEt.at(idx) >= 18)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_14262501724811299635
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_14262501725482388275
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG15: ET >= 30 at BX = 0
      if (not (data->egIEt.at(idx) >= 30)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_14262501742393822003
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG23: ET >= 46 at BX = 0
      if (not (data->egIEt.at(idx) >= 46)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_14262501742796475187
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_14262501759707908915
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG34: ET >= 68 at BX = 0
      if (not (data->egIEt.at(idx) >= 68)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_14262501759976344371
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG36: ET >= 72 at BX = 0
      if (not (data->egIEt.at(idx) >= 72)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_14262501760244779827
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG38: ET >= 76 at BX = 0
      if (not (data->egIEt.at(idx) >= 76)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_14262501776350907187
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG40: ET >= 80 at BX = 0
      if (not (data->egIEt.at(idx) >= 80)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_14262501776619342643
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG42: ET >= 84 at BX = 0
      if (not (data->egIEt.at(idx) >= 84)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_14262501777021995827
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG45: ET >= 90 at BX = 0
      if (not (data->egIEt.at(idx) >= 90)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_145873584
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG50: ET >= 100 at BX = 0
      if (not (data->egIEt.at(idx) >= 100)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_145873712
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG60: ET >= 120 at BX = 0
      if (not (data->egIEt.at(idx) >= 120)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244734408399686654
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG20: ET >= 40 at BX = 0
      if (not (data->egIEt.at(idx) >= 40)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244737706934569982
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG23: ET >= 46 at BX = 0
      if (not (data->egIEt.at(idx) >= 46)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244738805910375166
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244738805910375422
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244741004933630718
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244741004933630974
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244741005469453054
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244743203956886270
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244743203956886526
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244743204492708606
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244875145352219390
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(idx) >= 60)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244875145352219646
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(idx) >= 60)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244875145888041726
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(idx) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244877344375474942
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG32: ET >= 64 at BX = 0
      if (not (data->egIEt.at(idx) >= 64)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244877344911297278
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG32: ET >= 64 at BX = 0
      if (not (data->egIEt.at(idx) >= 64)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_9244879543934552830
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG34: ET >= 68 at BX = 0
      if (not (data->egIEt.at(idx) >= 68)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      
bool
SingleETMHF_306372248967728
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF40: ET >= 80 at BX = 0
      if (not (data->sumIEt.at(ii) >= 80)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_306372248967856
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF50: ET >= 100 at BX = 0
      if (not (data->sumIEt.at(ii) >= 100)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_306372248967984
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF60: ET >= 120 at BX = 0
      if (not (data->sumIEt.at(ii) >= 120)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_306372248968240
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF80: ET >= 160 at BX = 0
      if (not (data->sumIEt.at(ii) >= 160)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_306372248968368
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF90: ET >= 180 at BX = 0
      if (not (data->sumIEt.at(ii) >= 180)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_39215647867820080
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF100: ET >= 200 at BX = 0
      if (not (data->sumIEt.at(ii) >= 200)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_39215647867820208
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF110: ET >= 220 at BX = 0
      if (not (data->sumIEt.at(ii) >= 220)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_39215647867820336
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF120: ET >= 240 at BX = 0
      if (not (data->sumIEt.at(ii) >= 240)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_39215647867820464
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF130: ET >= 260 at BX = 0
      if (not (data->sumIEt.at(ii) >= 260)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_39215647867820592
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF140: ET >= 280 at BX = 0
      if (not (data->sumIEt.at(ii) >= 280)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_39215647867820720
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF150: ET >= 300 at BX = 0
      if (not (data->sumIEt.at(ii) >= 300)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETM_2393532815664
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM120: ET >= 240 at BX = 0
      if (not (data->sumIEt.at(ii) >= 240)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETM_2393532816048
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM150: ET >= 300 at BX = 0
      if (not (data->sumIEt.at(ii) >= 300)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETT_306374079453232
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETT1200: ET >= 2400 at BX = 0
      if (not (data->sumIEt.at(ii) >= 2400)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETT_306374079518768
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETT1600: ET >= 3200 at BX = 0
      if (not (data->sumIEt.at(ii) >= 3200)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETT_306374081517616
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETT2000: ET >= 4000 at BX = 0
      if (not (data->sumIEt.at(ii) >= 4000)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleEXT_10333571492674155900
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_10333571493211026812
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_10333571493479462268
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_1189548080491112364
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_14860848214295529384
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_14860848214295529387
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_15141600570663550655
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_16249626042834147010
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_17417807877912935668
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_17960169865075597331
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_4108951444235007726
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_6102798787913291629
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_6102798788181727085
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_6102799243448260461
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_6106690317781795101
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_6106690317781795102
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_6106690317781795103
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_6106690317781795104
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_6909925150529645277
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_6909925150529645278
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_6909925150529645533
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_6909925150529645534
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_6953200472440552930
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_866206785869629780
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_866206786138065236
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_8736797827952386068
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_9794008929098471889
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_9794008929098471890
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_9794008929098472145
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_9794008929098472146
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_9945386644737729380
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_9945386645006164836
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_9960888781174681113
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_9960888781443116569
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleHTT_19504896816
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT60: ET >= 120 at BX = 0
      if (not (data->sumIEt.at(ii) >= 120)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626710576
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT100: ET >= 200 at BX = 0
      if (not (data->sumIEt.at(ii) >= 200)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626710832
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT120: ET >= 240 at BX = 0
      if (not (data->sumIEt.at(ii) >= 240)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626711344
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT160: ET >= 320 at BX = 0
      if (not (data->sumIEt.at(ii) >= 320)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626726960
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT200: ET >= 400 at BX = 0
      if (not (data->sumIEt.at(ii) >= 400)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626727216
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT220: ET >= 440 at BX = 0
      if (not (data->sumIEt.at(ii) >= 440)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626727472
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT240: ET >= 480 at BX = 0
      if (not (data->sumIEt.at(ii) >= 480)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626727600
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT250: ET >= 500 at BX = 0
      if (not (data->sumIEt.at(ii) >= 500)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626727605
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT255: ET >= 510 at BX = 0
      if (not (data->sumIEt.at(ii) >= 510)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626727728
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT260: ET >= 520 at BX = 0
      if (not (data->sumIEt.at(ii) >= 520)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626727984
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT280: ET >= 560 at BX = 0
      if (not (data->sumIEt.at(ii) >= 560)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626743344
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT300: ET >= 600 at BX = 0
      if (not (data->sumIEt.at(ii) >= 600)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626743600
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT320: ET >= 640 at BX = 0
      if (not (data->sumIEt.at(ii) >= 640)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626743856
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT340: ET >= 680 at BX = 0
      if (not (data->sumIEt.at(ii) >= 680)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626744112
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT360: ET >= 720 at BX = 0
      if (not (data->sumIEt.at(ii) >= 720)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626759728
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT400: ET >= 800 at BX = 0
      if (not (data->sumIEt.at(ii) >= 800)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_2496626760368
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT450: ET >= 900 at BX = 0
      if (not (data->sumIEt.at(ii) >= 900)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      


bool
SingleJET_11235106006895834903
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET8: ET >= 16 at BX = 0
      if (not (data->jetIEt.at(idx) >= 16)) continue;

                        // -2.958 <= eta <= -1.392
              etaWindow1 = ((-68 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -33));
            
                        // 1.392 <= eta <= 2.958
              etaWindow2 = ((32 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 67));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_11401653256114550551
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET10: ET >= 20 at BX = 0
      if (not (data->jetIEt.at(idx) >= 20)) continue;

                        // -2.958 <= eta <= -1.392
              etaWindow1 = ((-68 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -33));
            
                        // 1.392 <= eta <= 2.958
              etaWindow2 = ((32 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 67));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_11401653256131327767
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET12: ET >= 24 at BX = 0
      if (not (data->jetIEt.at(idx) >= 24)) continue;

                        // -2.958 <= eta <= -1.392
              etaWindow1 = ((-68 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -33));
            
                        // 1.392 <= eta <= 2.958
              etaWindow2 = ((32 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 67));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_15945345163595414977
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET20: ET >= 40 at BX = 0
      if (not (data->jetIEt.at(idx) >= 40)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_15945345163599774657
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET20: ET >= 40 at BX = 0
      if (not (data->jetIEt.at(idx) >= 40)) continue;

                        // 3.0015 <= eta <= 5.0
              etaWindow1 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_15945345172520893889
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_15945345172525253569
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // 3.0015 <= eta <= 5.0
              etaWindow1 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_15945345197955153345
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_15945345197959513025
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // 3.0015 <= eta <= 5.0
              etaWindow1 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_15945345223724957121
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(idx) >= 180)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_15945345223729316801
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(idx) >= 180)) continue;

                        // 3.0015 <= eta <= 5.0
              etaWindow1 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_20010310069
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_20010310448
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_20010310832
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(idx) >= 180)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_2022287695184199115
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_2022287695188558795
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

                        // 3.0015 <= eta <= 5.0
              etaWindow1 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_2561319655605
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET115: ET >= 230 at BX = 0
      if (not (data->jetIEt.at(idx) >= 230)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_2561319655728
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_2561319656496
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET180: ET >= 360 at BX = 0
      if (not (data->jetIEt.at(idx) >= 360)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_2561319671856
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET200: ET >= 400 at BX = 0
      if (not (data->jetIEt.at(idx) >= 400)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_3448182530626688965
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx) >= 200)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_3448186928673200069
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_3448191326719711173
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET140: ET >= 280 at BX = 0
      if (not (data->jetIEt.at(idx) >= 280)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_3448195724766222277
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET160: ET >= 320 at BX = 0
      if (not (data->jetIEt.at(idx) >= 320)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_3448200122812733381
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET180: ET >= 360 at BX = 0
      if (not (data->jetIEt.at(idx) >= 360)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_7529292616000999046
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET20: ET >= 40 at BX = 0
      if (not (data->jetIEt.at(idx) >= 40)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_7529294815024254598
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx) >= 60)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_7529294900923600518
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_7529297065587117702
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET43: ET >= 86 at BX = 0
      if (not (data->jetIEt.at(idx) >= 86)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_7529297117126725254
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET46: ET >= 92 at BX = 0
      if (not (data->jetIEt.at(idx) >= 92)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_7529301412094021254
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_7529308009163787910
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(idx) >= 180)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      
bool
SingleMBT0HFM_43640316738250417
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMinBiasHFM0)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // MBT0HFM1: Count >= 1 at BX = 0
      if (not (data->sumIEt.at(ii) >= 1)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleMBT0HFP_43640316738250801
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMinBiasHFP0)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // MBT0HFP1: Count >= 1 at BX = 0
      if (not (data->sumIEt.at(ii) >= 1)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

      


bool
SingleMU_10359208142105920998
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.1038125 <= eta <= 1.1038125
              etaWindow1 = ((-101 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 101));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_10359208142105921382
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.4083125 <= eta <= 1.4083125
              etaWindow1 = ((-129 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 129));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_10359639116570609638
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_10360061329035675622
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU6: ET >= 13 at BX = 0
      if (not (data->muonIEt.at(idx) >= 13)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_10360202066524030950
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx) >= 15)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_10360342804012386278
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU8: ET >= 17 at BX = 0
      if (not (data->muonIEt.at(idx) >= 17)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_10360483541500741606
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU9: ET >= 19 at BX = 0
      if (not (data->muonIEt.at(idx) >= 19)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_1272496
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_14243093768255232179
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // -0.7993125 <= eta <= 0.7993125
              etaWindow1 = ((-73 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 73));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_14769293015645015621
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_14769293018627052229
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_14769293071236239813
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_14769293105595978181
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_14769293122775847365
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU6: ET >= 13 at BX = 0
      if (not (data->muonIEt.at(idx) >= 13)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_14769293135904099909
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx) >= 15)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_14769293139955716549
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx) >= 15)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_14769293157135585733
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU8: ET >= 17 at BX = 0
      if (not (data->muonIEt.at(idx) >= 17)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_16260934496621930532
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -0.7993125 <= eta <= 0.7993125
              etaWindow1 = ((-73 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 73));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_17545683106981072453
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU15: ET >= 31 at BX = 0
      if (not (data->muonIEt.at(idx) >= 31)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_17545683162572296645
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU18: ET >= 37 at BX = 0
      if (not (data->muonIEt.at(idx) >= 37)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_17545685224156598725
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU20: ET >= 41 at BX = 0
      if (not (data->muonIEt.at(idx) >= 41)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_17545685258516337093
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU22: ET >= 45 at BX = 0
      if (not (data->muonIEt.at(idx) >= 45)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_17545685310055944645
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU25: ET >= 51 at BX = 0
      if (not (data->muonIEt.at(idx) >= 51)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_5290897791608380091
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // 1.2451875 <= eta <= 2.45
              etaWindow1 = ((115 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 225));
            
                        // -2.45 <= eta <= -1.2451875
              etaWindow2 = ((-225 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -115));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_6011484727103937211
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // 0.7993125 <= eta <= 1.2451875
              etaWindow1 = ((74 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 114));
            
                        // -1.2451875 <= eta <= -0.7993125
              etaWindow2 = ((-114 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -74));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_6225176159725710459
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // 1.2451875 <= eta <= 2.45
              etaWindow1 = ((115 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 225));
            
                        // -2.45 <= eta <= -1.2451875
              etaWindow2 = ((-225 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -115));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_6225176160372139651
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU22: ET >= 45 at BX = 0
      if (not (data->muonIEt.at(idx) >= 45)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // 1.2451875 <= eta <= 2.45
              etaWindow1 = ((115 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 225));
            
                        // -2.45 <= eta <= -1.2451875
              etaWindow2 = ((-225 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -115));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_6945763095221267579
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // 0.7993125 <= eta <= 1.2451875
              etaWindow1 = ((74 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 114));
            
                        // -1.2451875 <= eta <= -0.7993125
              etaWindow2 = ((-114 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -74));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_6945763095867696771
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU22: ET >= 45 at BX = 0
      if (not (data->muonIEt.at(idx) >= 45)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // 0.7993125 <= eta <= 1.2451875
              etaWindow1 = ((74 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 114));
            
                        // -1.2451875 <= eta <= -0.7993125
              etaWindow2 = ((-114 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -74));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_7069342828816371872
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(idx) >= 25)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

                        // -0.7993125 <= eta <= 0.7993125
              etaWindow1 = ((-73 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 73));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_7181677643621025184
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU18: ET >= 37 at BX = 0
      if (not (data->muonIEt.at(idx) >= 37)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_7270359269352285314
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(idx) >= 25)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

                        // 1.2451875 <= eta <= 2.45
              etaWindow1 = ((115 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 225));
            
                        // -2.45 <= eta <= -1.2451875
              etaWindow2 = ((-225 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -115));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_7990946204847842434
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(idx) >= 25)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

                        // 0.7993125 <= eta <= 1.2451875
              etaWindow1 = ((74 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 114));
            
                        // -1.2451875 <= eta <= -0.7993125
              etaWindow2 = ((-114 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -74));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_9379434261777827232
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU22: ET >= 45 at BX = 0
      if (not (data->muonIEt.at(idx) >= 45)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_9379434265999970464
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU22: ET >= 45 at BX = 0
      if (not (data->muonIEt.at(idx) >= 45)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -0.7993125 <= eta <= 0.7993125
              etaWindow1 = ((-73 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 73));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_9710698557764193463
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU10: ET >= 21 at BX = 0
      if (not (data->muonIEt.at(idx) >= 21)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_9710980032740904119
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(idx) >= 25)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_9711261507717614775
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU14: ET >= 29 at BX = 0
      if (not (data->muonIEt.at(idx) >= 29)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_9711542982694325431
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU16: ET >= 33 at BX = 0
      if (not (data->muonIEt.at(idx) >= 33)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_9711824457671036087
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU18: ET >= 37 at BX = 0
      if (not (data->muonIEt.at(idx) >= 37)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061875 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_10012632024376351534
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU36: ET >= 72 at BX = 0
      if (not (data->tauIEt.at(idx) >= 72)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_12210388642533153582
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU40: ET >= 80 at BX = 0
      if (not (data->tauIEt.at(idx) >= 80)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_14552260448765811502
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU52: ET >= 104 at BX = 0
      if (not (data->tauIEt.at(idx) >= 104)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_16608831008486494024
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU24: ET >= 48 at BX = 0
      if (not (data->tauIEt.at(idx) >= 48)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_16608831042846232392
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU26: ET >= 52 at BX = 0
      if (not (data->tauIEt.at(idx) >= 52)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_16608841934883295048
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU70: ET >= 140 at BX = 0
      if (not (data->tauIEt.at(idx) >= 140)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_3484215725702552004
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU120: ET >= 240 at BX = 0
      if (not (data->tauIEt.at(idx) >= 240)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_3484217924725807556
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU130: ET >= 260 at BX = 0
      if (not (data->tauIEt.at(idx) >= 260)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_9940574430338423598
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU32: ET >= 64 at BX = 0
      if (not (data->tauIEt.at(idx) >= 64)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_9976603227357387566
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU34: ET >= 68 at BX = 0
      if (not (data->tauIEt.at(idx) >= 68)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

                    




bool
TransverseMass_15283531744994796576
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj = 0;
  bool etaWindow1;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
          
                                  // EG32: ET >= 64 at BX = 0
      if (not (data->egIEt.at(ii) >= 64)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEt)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
                          // ETM10: ET >= 20 at BX = 0
      if (not (data->sumIEt.at(jj) >= 20)) continue;
      
          long long minimum;
  long long maximum;
        // 40.0 <= Mt <= 151982.0
      int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->sumIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long cosDeltaPhi = LUT_COS_DPHI_EG_ETM[deltaIPhi];
  long long pt0 = LUT_EG_ET[data->egIEt.at(ii)];
  long long pt1 = LUT_ETM_ET[data->sumIEt.at(jj)];
  long long mt2 = pt0*pt1*(1*POW10[3] - cosDeltaPhi);
    minimum = (long long)(800.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mt2) and (mt2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}
    
                    




bool
TransverseMass_15283531744994797088
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj = 0;
  bool etaWindow1;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
          
                                  // EG32: ET >= 64 at BX = 0
      if (not (data->egIEt.at(ii) >= 64)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEt)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
                          // ETM10: ET >= 20 at BX = 0
      if (not (data->sumIEt.at(jj) >= 20)) continue;
      
          long long minimum;
  long long maximum;
        // 44.0 <= Mt <= 151982.0
      int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->sumIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long cosDeltaPhi = LUT_COS_DPHI_EG_ETM[deltaIPhi];
  long long pt0 = LUT_EG_ET[data->egIEt.at(ii)];
  long long pt1 = LUT_ETM_ET[data->sumIEt.at(jj)];
  long long mt2 = pt0*pt1*(1*POW10[3] - cosDeltaPhi);
    minimum = (long long)(968.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mt2) and (mt2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}
    
                    




bool
TransverseMass_15283531744994797600
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj = 0;
  bool etaWindow1;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
          
                                  // EG32: ET >= 64 at BX = 0
      if (not (data->egIEt.at(ii) >= 64)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEt)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
                          // ETM10: ET >= 20 at BX = 0
      if (not (data->sumIEt.at(jj) >= 20)) continue;
      
          long long minimum;
  long long maximum;
        // 48.0 <= Mt <= 151982.0
      int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->sumIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long cosDeltaPhi = LUT_COS_DPHI_EG_ETM[deltaIPhi];
  long long pt0 = LUT_EG_ET[data->egIEt.at(ii)];
  long long pt1 = LUT_ETM_ET[data->sumIEt.at(jj)];
  long long mt2 = pt0*pt1*(1*POW10[3] - cosDeltaPhi);
    minimum = (long long)(1152.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mt2) and (mt2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}
    
      


bool
TripleEG_10417466438308546634
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG18: ET >= 36 at BX = 0
      if (not (data->egIEt.at(idx) >= 36)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG18: ET >= 36 at BX = 0
      if (not (data->egIEt.at(idx) >= 36)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleEG_10417466496290605130
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG16: ET >= 32 at BX = 0
      if (not (data->egIEt.at(idx) >= 32)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG16: ET >= 32 at BX = 0
      if (not (data->egIEt.at(idx) >= 32)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // EG16: ET >= 32 at BX = 0
      if (not (data->egIEt.at(idx) >= 32)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleEG_12248554337091284411
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG16: ET >= 32 at BX = 0
      if (not (data->egIEt.at(idx) >= 32)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleEG_12248554337191947707
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG16: ET >= 32 at BX = 0
      if (not (data->egIEt.at(idx) >= 32)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG15: ET >= 30 at BX = 0
      if (not (data->egIEt.at(idx) >= 30)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleEG_12248554337275833787
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG18: ET >= 36 at BX = 0
      if (not (data->egIEt.at(idx) >= 36)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG17: ET >= 34 at BX = 0
      if (not (data->egIEt.at(idx) >= 34)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

                        // -2.523 <= eta <= 2.523
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleJET_15795555309634017625
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx) >= 200)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx) >= 160)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET70: ET >= 140 at BX = 0
      if (not (data->jetIEt.at(idx) >= 140)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleJET_15798370060072213465
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET105: ET >= 210 at BX = 0
      if (not (data->jetIEt.at(idx) >= 210)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET85: ET >= 170 at BX = 0
      if (not (data->jetIEt.at(idx) >= 170)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET75: ET >= 150 at BX = 0
      if (not (data->jetIEt.at(idx) >= 150)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleJET_2067252540421294278
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET95: ET >= 190 at BX = 0
      if (not (data->jetIEt.at(idx) >= 190)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET75: ET >= 150 at BX = 0
      if (not (data->jetIEt.at(idx) >= 150)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET65: ET >= 130 at BX = 0
      if (not (data->jetIEt.at(idx) >= 130)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_15692838580664758508
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3p5: ET >= 8 at BX = 0
      if (not (data->muonIEt.at(idx) >= 8)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU2p5: ET >= 6 at BX = 0
      if (not (data->muonIEt.at(idx) >= 6)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_3324682852515662879
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_3324683353497813023
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_3324683533187258399
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_3324685351042537503
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_3324685732743223327
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_3324691511169731615
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_3324691786047638559
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_3324692191841327135
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_3324692466719234079
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_3324692885559266335
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_6936497366859389375
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx) >= 9)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU2p5: ET >= 6 at BX = 0
      if (not (data->muonIEt.at(idx) >= 6)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_9287399899537551596
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3p5: ET >= 8 at BX = 0
      if (not (data->muonIEt.at(idx) >= 8)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU2p5: ET >= 6 at BX = 0
      if (not (data->muonIEt.at(idx) >= 6)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

  

// generate algorithms
bool
L1_AlwaysTrue(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_1189548080491112364(data) or ( not SingleEXT_1189548080491112364(data));
}
bool
L1_BPTX_AND_Ref1_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_10333571492674155900(data);
}
bool
L1_BPTX_AND_Ref3_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_10333571493211026812(data);
}
bool
L1_BPTX_AND_Ref4_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_10333571493479462268(data);
}
bool
L1_BPTX_BeamGas_B1_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_9960888781174681113(data);
}
bool
L1_BPTX_BeamGas_B2_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_9960888781443116569(data);
}
bool
L1_BPTX_BeamGas_Ref1_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_866206785869629780(data);
}
bool
L1_BPTX_BeamGas_Ref2_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_866206786138065236(data);
}
bool
L1_BPTX_NotOR_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_17417807877912935668(data);
}
bool
L1_BPTX_OR_Ref3_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_9945386644737729380(data);
}
bool
L1_BPTX_OR_Ref4_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_9945386645006164836(data);
}
bool
L1_BPTX_RefAND_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_15141600570663550655(data);
}
bool
L1_BptxMinus(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_6102798788181727085(data);
}
bool
L1_BptxOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_6102799243448260461(data);
}
bool
L1_BptxPlus(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_6102798787913291629(data);
}
bool
L1_BptxXOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return (SingleEXT_6102798787913291629(data) and ( not SingleEXT_6102798788181727085(data))) or (SingleEXT_6102798788181727085(data) and ( not SingleEXT_6102798787913291629(data)));
}
bool
L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_5013507948943010765(data);
}
bool
L1_DoubleEG8er2p5_HTT260er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_13299746526184647851(data) and SingleHTT_2496626727728(data);
}
bool
L1_DoubleEG8er2p5_HTT280er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_13299746526184647851(data) and SingleHTT_2496626727984(data);
}
bool
L1_DoubleEG8er2p5_HTT300er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_13299746526184647851(data) and SingleHTT_2496626743344(data);
}
bool
L1_DoubleEG8er2p5_HTT320er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_13299746526184647851(data) and SingleHTT_2496626743600(data);
}
bool
L1_DoubleEG8er2p5_HTT340er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_13299746526184647851(data) and SingleHTT_2496626743856(data);
}
bool
L1_DoubleEG_15_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_7286147403649300475(data);
}
bool
L1_DoubleEG_20_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_7286147931930277883(data);
}
bool
L1_DoubleEG_22_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_7286147940520212475(data);
}
bool
L1_DoubleEG_25_12_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_7286147987764852731(data);
}
bool
L1_DoubleEG_25_14_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_7286148022124591099(data);
}
bool
L1_DoubleEG_27_14_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_7286148030714525691(data);
}
bool
L1_DoubleEG_LooseIso20_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_56167094600470587(data);
}
bool
L1_DoubleEG_LooseIso22_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_56237463344648251(data);
}
bool
L1_DoubleEG_LooseIso22_12_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_56237497704386619(data);
}
bool
L1_DoubleEG_LooseIso25_12_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_56343050820653115(data);
}
bool
L1_DoubleIsoTau32er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_14812192849374407856(data);
}
bool
L1_DoubleIsoTau34er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_977134795165985969(data);
}
bool
L1_DoubleIsoTau36er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_5588820814667115697(data);
}
bool
L1_DoubleJet100er2p3_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_7041035331702023693(data);
}
bool
L1_DoubleJet100er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_7222390773192019523(data);
}
bool
L1_DoubleJet112er2p3_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_7041035331710545453(data);
}
bool
L1_DoubleJet120er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_7222953723145440851(data);
}
bool
L1_DoubleJet150er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_7223798148075572843(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_3424918508618058399(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_3425188988478491295(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_3425199983594769055(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_3425470463455201951(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_3425477060524968607(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_3425483657594735263(data);
}
bool
L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMassOvRm_10967205787862279205(data);
}
bool
L1_DoubleJet40er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_14043965433308858268(data);
}
bool
L1_DoubleJet_100_30_DoubleJet30_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_16307690244847013269(data) and InvariantMass_2940638391876117895(data);
}
bool
L1_DoubleJet_110_35_DoubleJet35_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_16379747838884941845(data) and InvariantMass_2940649386995017095(data);
}
bool
L1_DoubleJet_115_40_DoubleJet40_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_16382562588652064149(data) and InvariantMass_2940919866919937415(data);
}
bool
L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_2561319655605(data) and (InvariantMass_9620562409477107818(data) or InvariantMass_2941482817007576455(data) or InvariantMass_3915066037994721069(data) or InvariantMass_14031536345113029771(data) or InvariantMass_7789845660996549812(data) or InvariantMass_16026464453068411612(data));
}
bool
L1_DoubleJet_120_45_DoubleJet45_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_16451805432922886165(data) and InvariantMass_2940930862038836615(data);
}
bool
L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_2561319655728(data) and (InvariantMass_9620565158256341098(data) or InvariantMass_2941482817007576455(data) or InvariantMass_3915066038015692589(data) or InvariantMass_14031539093892099211(data) or InvariantMass_12401534428203007157(data) or InvariantMass_16026640374949827292(data));
}
bool
L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_2940638392207467911(data) and DoubleMU_14585778268989477695(data);
}
bool
L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMassOvRm_10944688471548334117(data);
}
bool
L1_DoubleJet_80_30_Mass_Min420_Mu8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_2940638392207467911(data) and SingleMU_14769293157135585733(data);
}
bool
L1_DoubleJet_90_30_DoubleJet30_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_4162612533456677351(data) and InvariantMass_2940638391876117895(data);
}
bool
L1_DoubleLooseIsoEG22er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_2355036583129339571(data);
}
bool
L1_DoubleLooseIsoEG24er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_2931778810409473715(data);
}
bool
L1_DoubleMu0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14585777620730815295(data);
}
bool
L1_DoubleMu0_Mass_Min1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_7551966572230413519(data);
}
bool
L1_DoubleMu0_OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14585778097856672575(data);
}
bool
L1_DoubleMu0_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14585778268989477695(data);
}
bool
L1_DoubleMu0_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_3139255731352238604(data);
}
bool
L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_15498326754036298551(data) and CaloMuonCorrelation_8802031396140112104(data);
}
bool
L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_9513481109949270451(data);
}
bool
L1_DoubleMu0er1p5_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_3947425258490794706(data);
}
bool
L1_DoubleMu0er1p5_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14617142003772573591(data);
}
bool
L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_9513481109957663155(data);
}
bool
L1_DoubleMu0er1p5_SQ_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_15199048927445899759(data);
}
bool
L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_9513481247421761971(data);
}
bool
L1_DoubleMu0er2p0_SQ_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_15199048929593776303(data);
}
bool
L1_DoubleMu10_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_12394181525097886766(data);
}
bool
L1_DoubleMu18er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_16323903523977050720(data);
}
bool
L1_DoubleMu3_OS_DoubleEG7p5Upsilon(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_4461482972834602413(data) and InvariantMass_16981538589298500419(data);
}
bool
L1_DoubleMu3_SQ_ETMHF50_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14585786515326686015(data) and SingleETMHF_306372248967856(data) and SingleHTT_19504896816(data);
}
bool
L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14585786515326686015(data) and SingleETMHF_306372248967856(data) and SingleJET_7529301412094021254(data);
}
bool
L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14585786515326686015(data) and SingleETMHF_306372248967856(data) and (SingleJET_7529301412094021254(data) or DoubleJET_14043965433308858268(data));
}
bool
L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14585786515326686015(data) and SingleETMHF_306372248967984(data) and SingleJET_7529301412094021254(data);
}
bool
L1_DoubleMu3_SQ_HTT220er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14585786515326686015(data) and SingleHTT_2496626727216(data);
}
bool
L1_DoubleMu3_SQ_HTT240er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14585786515326686015(data) and SingleHTT_2496626727472(data);
}
bool
L1_DoubleMu3_SQ_HTT260er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14585786515326686015(data) and SingleHTT_2496626727728(data);
}
bool
L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_5698493964878099256(data) and CaloMuonCorrelation_8802031396140113640(data);
}
bool
L1_DoubleMu4_SQ_EG9er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14585789264105755455(data) and SingleEG_14243075932536834867(data);
}
bool
L1_DoubleMu4_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_3229327723899648524(data);
}
bool
L1_DoubleMu4_SQ_OS_dR_Max1p2(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_16784489743460462578(data);
}
bool
L1_DoubleMu4p5_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_2011765979326275391(data);
}
bool
L1_DoubleMu4p5_SQ_OS_dR_Max1p2(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_12923126501326425857(data);
}
bool
L1_DoubleMu4p5er2p0_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_13627348644483379947(data);
}
bool
L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_2342552854377181621(data);
}
bool
L1_DoubleMu5Upsilon_OS_DoubleEG3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_13689376201502793133(data) and InvariantMass_2443380592745462540(data);
}
bool
L1_DoubleMu5_SQ_EG9er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14585792012884824895(data) and SingleEG_14243075932536834867(data);
}
bool
L1_DoubleMu9_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_14585803008001102655(data);
}
bool
L1_DoubleMu_12_5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_16961157256621881348(data);
}
bool
L1_DoubleMu_15_5_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_16961159554147985412(data);
}
bool
L1_DoubleMu_15_7(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_16961163303935834116(data);
}
bool
L1_DoubleMu_15_7_Mass_Min1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_12340992865530516744(data);
}
bool
L1_DoubleMu_15_7_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_16961163952194496516(data);
}
bool
L1_DoubleTau70er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_17539608616528615651(data);
}
bool
L1_ETM120(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETM_2393532815664(data);
}
bool
L1_ETM150(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETM_2393532816048(data);
}
bool
L1_ETMHF100(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_39215647867820080(data);
}
bool
L1_ETMHF100_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_39215647867820080(data) and SingleHTT_19504896816(data);
}
bool
L1_ETMHF110(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_39215647867820208(data);
}
bool
L1_ETMHF110_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_39215647867820208(data) and SingleHTT_19504896816(data);
}
bool
L1_ETMHF110_HTT60er_NotSecondBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_39215647867820208(data) and SingleHTT_19504896816(data) and ((SingleEXT_6909925150529645534(data)) or ( not SingleEXT_9794008929098472145(data)) or ( not SingleEXT_1189548080491112364(data)) or ( not SingleEXT_9794008929098471889(data)) or ( not SingleEXT_9794008929098471890(data)));
}
bool
L1_ETMHF120(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_39215647867820336(data);
}
bool
L1_ETMHF120_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_39215647867820336(data) and SingleHTT_19504896816(data);
}
bool
L1_ETMHF120_NotSecondBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_39215647867820336(data) and ((SingleEXT_6909925150529645534(data)) or ( not SingleEXT_9794008929098472145(data)) or ( not SingleEXT_1189548080491112364(data)) or ( not SingleEXT_9794008929098471889(data)) or ( not SingleEXT_9794008929098471890(data)));
}
bool
L1_ETMHF130(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_39215647867820464(data);
}
bool
L1_ETMHF130_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_39215647867820464(data) and SingleHTT_19504896816(data);
}
bool
L1_ETMHF140(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_39215647867820592(data);
}
bool
L1_ETMHF150(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_39215647867820720(data);
}
bool
L1_ETMHF90_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_306372248968368(data) and SingleHTT_19504896816(data);
}
bool
L1_ETT1200(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETT_306374079453232(data);
}
bool
L1_ETT1600(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETT_306374079518768(data);
}
bool
L1_ETT2000(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETT_306374081517616(data);
}
bool
L1_FirstBunchAfterTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_9794008929098472146(data) and SingleEXT_9794008929098472145(data) and ( not SingleEXT_6102799243448260461(data)) and ( not SingleEXT_6909925150529645277(data)) and ( not SingleEXT_6909925150529645278(data));
}
bool
L1_FirstBunchBeforeTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return ( not SingleEXT_6909925150529645534(data)) and ( not SingleEXT_6909925150529645533(data)) and ( not SingleEXT_6102799243448260461(data)) and SingleEXT_9794008929098471889(data) and SingleEXT_9794008929098471890(data);
}
bool
L1_FirstBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return ( not SingleEXT_6909925150529645534(data)) and ( not SingleEXT_6909925150529645533(data)) and SingleEXT_1189548080491112364(data) and SingleEXT_9794008929098471889(data) and SingleEXT_9794008929098471890(data);
}
bool
L1_FirstCollisionInOrbit(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_4108951444235007726(data);
}
bool
L1_FirstCollisionInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_16249626042834147010(data);
}
bool
L1_HCAL_LaserMon_Trig(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_14860848214295529384(data);
}
bool
L1_HCAL_LaserMon_Veto(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_14860848214295529387(data);
}
bool
L1_HTT120er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_2496626710832(data);
}
bool
L1_HTT160er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_2496626711344(data);
}
bool
L1_HTT200er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_2496626726960(data);
}
bool
L1_HTT255er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_2496626727605(data);
}
bool
L1_HTT280er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_2496626727984(data);
}
bool
L1_HTT280er_QuadJet_70_55_40_35_er2p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_2496626727984(data) and QuadJET_2969440952489109684(data);
}
bool
L1_HTT320er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_2496626743600(data);
}
bool
L1_HTT320er_QuadJet_70_55_40_40_er2p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_2496626743600(data) and QuadJET_2969443065613019316(data);
}
bool
L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_2496626743600(data) and QuadJET_1795783097020010207(data);
}
bool
L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_2496626743600(data) and QuadJET_1795850802884464351(data);
}
bool
L1_HTT360er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_2496626744112(data);
}
bool
L1_HTT400er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_2496626759728(data);
}
bool
L1_HTT450er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_2496626760368(data);
}
bool
L1_IsoEG32er2p5_Mt40(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TransverseMass_15283531744994796576(data);
}
bool
L1_IsoEG32er2p5_Mt44(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TransverseMass_15283531744994797088(data);
}
bool
L1_IsoEG32er2p5_Mt48(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TransverseMass_15283531744994797600(data);
}
bool
L1_IsoTau40er2p1_ETMHF100(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_12210388642533153582(data) and SingleETMHF_39215647867820080(data);
}
bool
L1_IsoTau40er2p1_ETMHF110(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_12210388642533153582(data) and SingleETMHF_39215647867820208(data);
}
bool
L1_IsoTau40er2p1_ETMHF120(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_12210388642533153582(data) and SingleETMHF_39215647867820336(data);
}
bool
L1_IsoTau40er2p1_ETMHF90(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_12210388642533153582(data) and SingleETMHF_306372248968368(data);
}
bool
L1_IsolatedBunch(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return ( not SingleEXT_6909925150529645534(data)) and ( not SingleEXT_6909925150529645533(data)) and SingleEXT_1189548080491112364(data) and ( not SingleEXT_6909925150529645277(data)) and ( not SingleEXT_6909925150529645278(data));
}
bool
L1_LastBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_9794008929098472146(data) and SingleEXT_9794008929098472145(data) and SingleEXT_1189548080491112364(data) and ( not SingleEXT_6909925150529645277(data)) and ( not SingleEXT_6909925150529645278(data));
}
bool
L1_LastCollisionInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_6953200472440552930(data);
}
bool
L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_980305370107485378(data);
}
bool
L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_18379087122140179561(data);
}
bool
L1_LooseIsoEG24er2p1_HTT100er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244738805910375422(data) and SingleHTT_2496626710576(data);
}
bool
L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_980305438827093186(data);
}
bool
L1_LooseIsoEG26er2p1_HTT100er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244741004933630974(data) and SingleHTT_2496626710576(data);
}
bool
L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_9154320878437213226(data);
}
bool
L1_LooseIsoEG28er2p1_HTT100er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244743203956886526(data) and SingleHTT_2496626710576(data);
}
bool
L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_9154320878437278762(data);
}
bool
L1_LooseIsoEG30er2p1_HTT100er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244875145352219646(data) and SingleHTT_2496626710576(data);
}
bool
L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_9154320878441210922(data);
}
bool
L1_MinimumBiasHF0_AND_BptxAND(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return (SingleMBT0HFP_43640316738250801(data) and SingleMBT0HFM_43640316738250417(data)) and SingleEXT_1189548080491112364(data);
}
bool
L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_12420459208926927530(data) and CaloCaloCorrelation_3813196582576312175(data);
}
bool
L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_9314160007219977660(data) and CaloCaloCorrelation_1626121879521236767(data);
}
bool
L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_14690273421121723050(data) and CaloCaloCorrelation_3813196582576378703(data);
}
bool
L1_Mu18er2p1_Tau24er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_7181677643621025184(data) and SingleTAU_16608831008486494024(data);
}
bool
L1_Mu18er2p1_Tau26er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_7181677643621025184(data) and SingleTAU_16608831042846232392(data);
}
bool
L1_Mu20_EG10er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_17545685224156598725(data) and SingleEG_14262501724811299635(data);
}
bool
L1_Mu22er2p1_IsoTau32er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_9379434261777827232(data) and SingleTAU_9940574430338423598(data);
}
bool
L1_Mu22er2p1_IsoTau34er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_9379434261777827232(data) and SingleTAU_9976603227357387566(data);
}
bool
L1_Mu22er2p1_IsoTau36er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_9379434261777827232(data) and SingleTAU_10012632024376351534(data);
}
bool
L1_Mu22er2p1_IsoTau40er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_9379434261777827232(data) and SingleTAU_12210388642533153582(data);
}
bool
L1_Mu22er2p1_Tau70er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_9379434261777827232(data) and SingleTAU_16608841934883295048(data);
}
bool
L1_Mu3_Jet120er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_2018052583404954969(data);
}
bool
L1_Mu3_Jet120er2p5_dR_Max0p8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_2018052583404955481(data);
}
bool
L1_Mu3_Jet16er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_17980860017930427617(data);
}
bool
L1_Mu3_Jet30er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293071236239813(data) and SingleJET_7529294815024254598(data);
}
bool
L1_Mu3_Jet35er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_15675017008716733697(data);
}
bool
L1_Mu3_Jet60er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_4145801962648263985(data);
}
bool
L1_Mu3_Jet80er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_4145801962648264017(data);
}
bool
L1_Mu3er1p5_Jet100er2p5_ETMHF40(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_10359639116570609638(data) and SingleJET_3448182530626688965(data) and SingleETMHF_306372248967728(data);
}
bool
L1_Mu3er1p5_Jet100er2p5_ETMHF50(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_10359639116570609638(data) and SingleJET_3448182530626688965(data) and SingleETMHF_306372248967856(data);
}
bool
L1_Mu5_EG23er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293105595978181(data) and SingleEG_14262501742393822003(data);
}
bool
L1_Mu5_LooseIsoEG20er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293105595978181(data) and SingleEG_9244734408399686654(data);
}
bool
L1_Mu6_DoubleEG10er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293122775847365(data) and DoubleEG_7286147382174463995(data);
}
bool
L1_Mu6_DoubleEG12er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293122775847365(data) and DoubleEG_7286147425124136955(data);
}
bool
L1_Mu6_DoubleEG15er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293122775847365(data) and DoubleEG_7286147489548646395(data);
}
bool
L1_Mu6_DoubleEG17er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293122775847365(data) and DoubleEG_7286147532498319355(data);
}
bool
L1_Mu6_HTT240er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293122775847365(data) and SingleHTT_2496626727472(data);
}
bool
L1_Mu6_HTT250er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293122775847365(data) and SingleHTT_2496626727600(data);
}
bool
L1_Mu7_EG23er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293139955716549(data) and SingleEG_14262501742393822003(data);
}
bool
L1_Mu7_LooseIsoEG20er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293139955716549(data) and SingleEG_9244734408399686654(data);
}
bool
L1_Mu7_LooseIsoEG23er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293139955716549(data) and SingleEG_9244737706934569982(data);
}
bool
L1_NotBptxOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return not SingleEXT_6102799243448260461(data);
}
bool
L1_QuadJet36er2p5_IsoTau52er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadJET_15326169532950703320(data) and SingleTAU_14552260448765811502(data);
}
bool
L1_QuadJet60er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadJET_284978008326942669(data);
}
bool
L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadJET_6278535506187355856(data) and DoubleJET_17548339888472803228(data) and (SingleJET_15945345163595414977(data) or SingleJET_15945345163599774657(data));
}
bool
L1_QuadMu0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadMU_509409160461874775(data);
}
bool
L1_QuadMu0_OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadMU_509409667408098135(data);
}
bool
L1_QuadMu0_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadMU_509409849236703575(data);
}
bool
L1_SecondBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return ( not SingleEXT_6909925150529645534(data)) and SingleEXT_9794008929098472145(data) and SingleEXT_1189548080491112364(data) and SingleEXT_9794008929098471889(data) and SingleEXT_9794008929098471890(data);
}
bool
L1_SecondLastBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_9794008929098472146(data) and SingleEXT_9794008929098472145(data) and SingleEXT_1189548080491112364(data) and SingleEXT_9794008929098471889(data) and ( not SingleEXT_6909925150529645278(data));
}
bool
L1_SingleEG10er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_14262501724811299635(data);
}
bool
L1_SingleEG15er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_14262501725482388275(data);
}
bool
L1_SingleEG26er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_14262501742796475187(data);
}
bool
L1_SingleEG34er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_14262501759707908915(data);
}
bool
L1_SingleEG36er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_14262501759976344371(data);
}
bool
L1_SingleEG38er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_14262501760244779827(data);
}
bool
L1_SingleEG40er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_14262501776350907187(data);
}
bool
L1_SingleEG42er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_14262501776619342643(data);
}
bool
L1_SingleEG45er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_14262501777021995827(data);
}
bool
L1_SingleEG50(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_145873584(data);
}
bool
L1_SingleEG60(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_145873712(data);
}
bool
L1_SingleEG8er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_14243075932402617139(data);
}
bool
L1_SingleIsoEG24er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_10104274574069012747(data);
}
bool
L1_SingleIsoEG24er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244738805910375166(data);
}
bool
L1_SingleIsoEG26er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_10104274582658947339(data);
}
bool
L1_SingleIsoEG26er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244741004933630718(data);
}
bool
L1_SingleIsoEG26er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244741005469453054(data);
}
bool
L1_SingleIsoEG28er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_10104274591248881931(data);
}
bool
L1_SingleIsoEG28er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244743203956886270(data);
}
bool
L1_SingleIsoEG28er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244743204492708606(data);
}
bool
L1_SingleIsoEG30er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244875145352219390(data);
}
bool
L1_SingleIsoEG30er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244875145888041726(data);
}
bool
L1_SingleIsoEG32er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244877344375474942(data);
}
bool
L1_SingleIsoEG32er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244877344911297278(data);
}
bool
L1_SingleIsoEG34er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_9244879543934552830(data);
}
bool
L1_SingleJet10erHE(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_11401653256114550551(data);
}
bool
L1_SingleJet120(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_2561319655728(data);
}
bool
L1_SingleJet120_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_2022287695184199115(data) or SingleJET_2022287695188558795(data);
}
bool
L1_SingleJet120er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_3448186928673200069(data);
}
bool
L1_SingleJet12erHE(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_11401653256131327767(data);
}
bool
L1_SingleJet140er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_3448191326719711173(data);
}
bool
L1_SingleJet140er2p5_ETMHF80(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_3448191326719711173(data) and SingleETMHF_306372248968240(data);
}
bool
L1_SingleJet140er2p5_ETMHF90(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_3448191326719711173(data) and SingleETMHF_306372248968368(data);
}
bool
L1_SingleJet160er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_3448195724766222277(data);
}
bool
L1_SingleJet180(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_2561319656496(data);
}
bool
L1_SingleJet180er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_3448200122812733381(data);
}
bool
L1_SingleJet200(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_2561319671856(data);
}
bool
L1_SingleJet20er2p5_NotBptxOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_7529292616000999046(data) and ( not SingleEXT_6102799243448260461(data));
}
bool
L1_SingleJet20er2p5_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_7529292616000999046(data) and ( not SingleEXT_6909925150529645533(data)) and ( not SingleEXT_6102799243448260461(data)) and ( not SingleEXT_6909925150529645277(data));
}
bool
L1_SingleJet35(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_20010310069(data);
}
bool
L1_SingleJet35_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_15945345172520893889(data) or SingleJET_15945345172525253569(data);
}
bool
L1_SingleJet35er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_7529294900923600518(data);
}
bool
L1_SingleJet43er2p5_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_7529297065587117702(data) and ( not SingleEXT_6909925150529645533(data)) and ( not SingleEXT_6102799243448260461(data)) and ( not SingleEXT_6909925150529645277(data));
}
bool
L1_SingleJet46er2p5_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_7529297117126725254(data) and ( not SingleEXT_6909925150529645533(data)) and ( not SingleEXT_6102799243448260461(data)) and ( not SingleEXT_6909925150529645277(data));
}
bool
L1_SingleJet60(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_20010310448(data);
}
bool
L1_SingleJet60_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_15945345197955153345(data) or SingleJET_15945345197959513025(data);
}
bool
L1_SingleJet60er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_7529301412094021254(data);
}
bool
L1_SingleJet8erHE(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_11235106006895834903(data);
}
bool
L1_SingleJet90(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_20010310832(data);
}
bool
L1_SingleJet90_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_15945345223724957121(data) or SingleJET_15945345223729316801(data);
}
bool
L1_SingleJet90er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_7529308009163787910(data);
}
bool
L1_SingleLooseIsoEG28er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_10104274591248882187(data);
}
bool
L1_SingleLooseIsoEG30er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_10104275106644957707(data);
}
bool
L1_SingleMu0_BMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_16260934496621930532(data);
}
bool
L1_SingleMu0_DQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293015645015621(data);
}
bool
L1_SingleMu0_EMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_6225176159725710459(data);
}
bool
L1_SingleMu0_OMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_6945763095221267579(data);
}
bool
L1_SingleMu10er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_9710698557764193463(data);
}
bool
L1_SingleMu12_DQ_BMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_7069342828816371872(data);
}
bool
L1_SingleMu12_DQ_EMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_7270359269352285314(data);
}
bool
L1_SingleMu12_DQ_OMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_7990946204847842434(data);
}
bool
L1_SingleMu12er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_9710980032740904119(data);
}
bool
L1_SingleMu14er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_9711261507717614775(data);
}
bool
L1_SingleMu15_DQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_17545683106981072453(data);
}
bool
L1_SingleMu16er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_9711542982694325431(data);
}
bool
L1_SingleMu18(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_17545683162572296645(data);
}
bool
L1_SingleMu18er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_9711824457671036087(data);
}
bool
L1_SingleMu20(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_17545685224156598725(data);
}
bool
L1_SingleMu22(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_17545685258516337093(data);
}
bool
L1_SingleMu22_BMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_9379434265999970464(data);
}
bool
L1_SingleMu22_EMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_6225176160372139651(data);
}
bool
L1_SingleMu22_OMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_6945763095867696771(data);
}
bool
L1_SingleMu25(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_17545685310055944645(data);
}
bool
L1_SingleMu3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293071236239813(data);
}
bool
L1_SingleMu5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293105595978181(data);
}
bool
L1_SingleMu6er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_10360061329035675622(data);
}
bool
L1_SingleMu7(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293139955716549(data);
}
bool
L1_SingleMu7_DQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293135904099909(data);
}
bool
L1_SingleMu7er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_10360202066524030950(data);
}
bool
L1_SingleMu8er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_10360342804012386278(data);
}
bool
L1_SingleMu9er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_10360483541500741606(data);
}
bool
L1_SingleMuCosmics(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_1272496(data);
}
bool
L1_SingleMuCosmics_BMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14243093768255232179(data);
}
bool
L1_SingleMuCosmics_EMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_5290897791608380091(data);
}
bool
L1_SingleMuCosmics_OMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_6011484727103937211(data);
}
bool
L1_SingleMuOpen(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293018627052229(data);
}
bool
L1_SingleMuOpen_NotBptxOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_14769293018627052229(data) and ( not SingleEXT_6102799243448260461(data));
}
bool
L1_SingleMuOpen_er1p1_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_10359208142105920998(data) and ( not SingleEXT_6909925150529645533(data)) and ( not SingleEXT_6102799243448260461(data)) and ( not SingleEXT_6909925150529645277(data));
}
bool
L1_SingleMuOpen_er1p4_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_10359208142105921382(data) and ( not SingleEXT_6909925150529645533(data)) and ( not SingleEXT_6102799243448260461(data)) and ( not SingleEXT_6909925150529645277(data));
}
bool
L1_SingleTau120er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_3484215725702552004(data);
}
bool
L1_SingleTau130er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_3484217924725807556(data);
}
bool
L1_TOTEM_1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_6106690317781795101(data);
}
bool
L1_TOTEM_2(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_6106690317781795102(data);
}
bool
L1_TOTEM_3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_6106690317781795103(data);
}
bool
L1_TOTEM_4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_6106690317781795104(data);
}
bool
L1_TripleEG16er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_10417466496290605130(data);
}
bool
L1_TripleEG_16_12_8_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_12248554337091284411(data);
}
bool
L1_TripleEG_16_15_8_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_12248554337191947707(data);
}
bool
L1_TripleEG_18_17_8_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_12248554337275833787(data);
}
bool
L1_TripleEG_18_18_12_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_10417466438308546634(data);
}
bool
L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleJET_15795555309634017625(data) and DoubleJET_209751802956826525(data);
}
bool
L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleJET_15798370060072213465(data) and DoubleJET_254798794346809245(data);
}
bool
L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleJET_2067252540421294278(data) and DoubleJET_17548339888472803228(data);
}
bool
L1_TripleMu0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_3324682852515662879(data);
}
bool
L1_TripleMu0_OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_3324683353497813023(data);
}
bool
L1_TripleMu0_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_3324683533187258399(data);
}
bool
L1_TripleMu3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_3324691511169731615(data);
}
bool
L1_TripleMu3_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_3324692191841327135(data);
}
bool
L1_TripleMu_5SQ_3SQ_0OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_3324685732743223327(data);
}
bool
L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_3324685732743223327(data) and InvariantMass_3063833799189854821(data);
}
bool
L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_3324685351042537503(data) and InvariantMass_3063833799189854821(data);
}
bool
L1_TripleMu_5_3_3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_3324691786047638559(data);
}
bool
L1_TripleMu_5_3_3_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_3324692466719234079(data);
}
bool
L1_TripleMu_5_3p5_2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_9287399899537551596(data);
}
bool
L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_9287399899537551596(data) and InvariantMass_15191958030943548804(data);
}
bool
L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_15692838580664758508(data) and InvariantMass_15192153509407276420(data);
}
bool
L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_6936497366859389375(data) and InvariantMass_15191958030943548804(data);
}
bool
L1_TripleMu_5_5_3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_3324692885559266335(data);
}
bool
L1_UnpairedBunchBptxMinus(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_8736797827952386068(data);
}
bool
L1_UnpairedBunchBptxPlus(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_17960169865075597331(data);
}
bool
L1_ZeroBias(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_1189548080491112364(data);
}
bool
L1_ZeroBias_copy(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_1189548080491112364(data);
}


std::string getNameFromId(const int index)
{
  static const std::pair<int, std::string> id2name[] = {
          std::make_pair(458, "L1_AlwaysTrue"),          std::make_pair(486, "L1_BPTX_AND_Ref1_VME"),          std::make_pair(487, "L1_BPTX_AND_Ref3_VME"),          std::make_pair(488, "L1_BPTX_AND_Ref4_VME"),          std::make_pair(491, "L1_BPTX_BeamGas_B1_VME"),          std::make_pair(492, "L1_BPTX_BeamGas_B2_VME"),          std::make_pair(489, "L1_BPTX_BeamGas_Ref1_VME"),          std::make_pair(490, "L1_BPTX_BeamGas_Ref2_VME"),          std::make_pair(482, "L1_BPTX_NotOR_VME"),          std::make_pair(483, "L1_BPTX_OR_Ref3_VME"),          std::make_pair(484, "L1_BPTX_OR_Ref4_VME"),          std::make_pair(485, "L1_BPTX_RefAND_VME"),          std::make_pair(467, "L1_BptxMinus"),          std::make_pair(464, "L1_BptxOR"),          std::make_pair(466, "L1_BptxPlus"),          std::make_pair(465, "L1_BptxXOR"),          std::make_pair(494, "L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142"),          std::make_pair(247, "L1_DoubleEG8er2p5_HTT260er"),          std::make_pair(248, "L1_DoubleEG8er2p5_HTT280er"),          std::make_pair(249, "L1_DoubleEG8er2p5_HTT300er"),          std::make_pair(250, "L1_DoubleEG8er2p5_HTT320er"),          std::make_pair(251, "L1_DoubleEG8er2p5_HTT340er"),          std::make_pair(205, "L1_DoubleEG_15_10_er2p5"),          std::make_pair(206, "L1_DoubleEG_20_10_er2p5"),          std::make_pair(207, "L1_DoubleEG_22_10_er2p5"),          std::make_pair(208, "L1_DoubleEG_25_12_er2p5"),          std::make_pair(209, "L1_DoubleEG_25_14_er2p5"),          std::make_pair(210, "L1_DoubleEG_27_14_er2p5"),          std::make_pair(212, "L1_DoubleEG_LooseIso20_10_er2p5"),          std::make_pair(213, "L1_DoubleEG_LooseIso22_10_er2p5"),          std::make_pair(214, "L1_DoubleEG_LooseIso22_12_er2p5"),          std::make_pair(215, "L1_DoubleEG_LooseIso25_12_er2p5"),          std::make_pair(275, "L1_DoubleIsoTau32er2p1"),          std::make_pair(276, "L1_DoubleIsoTau34er2p1"),          std::make_pair(277, "L1_DoubleIsoTau36er2p1"),          std::make_pair(345, "L1_DoubleJet100er2p3_dEta_Max1p6"),          std::make_pair(341, "L1_DoubleJet100er2p5"),          std::make_pair(346, "L1_DoubleJet112er2p3_dEta_Max1p6"),          std::make_pair(342, "L1_DoubleJet120er2p5"),          std::make_pair(343, "L1_DoubleJet150er2p5"),          std::make_pair(348, "L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5"),          std::make_pair(349, "L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5"),          std::make_pair(350, "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5"),          std::make_pair(351, "L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5"),          std::make_pair(352, "L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5"),          std::make_pair(353, "L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5"),          std::make_pair(363, "L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp"),          std::make_pair(340, "L1_DoubleJet40er2p5"),          std::make_pair(356, "L1_DoubleJet_100_30_DoubleJet30_Mass_Min620"),          std::make_pair(357, "L1_DoubleJet_110_35_DoubleJet35_Mass_Min620"),          std::make_pair(358, "L1_DoubleJet_115_40_DoubleJet40_Mass_Min620"),          std::make_pair(360, "L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28"),          std::make_pair(359, "L1_DoubleJet_120_45_DoubleJet45_Mass_Min620"),          std::make_pair(361, "L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28"),          std::make_pair(366, "L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ"),          std::make_pair(364, "L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp"),          std::make_pair(365, "L1_DoubleJet_80_30_Mass_Min420_Mu8"),          std::make_pair(355, "L1_DoubleJet_90_30_DoubleJet30_Mass_Min620"),          std::make_pair(217, "L1_DoubleLooseIsoEG22er2p1"),          std::make_pair(218, "L1_DoubleLooseIsoEG24er2p1"),          std::make_pair(40, "L1_DoubleMu0"),          std::make_pair(43, "L1_DoubleMu0_Mass_Min1"),          std::make_pair(39, "L1_DoubleMu0_OQ"),          std::make_pair(41, "L1_DoubleMu0_SQ"),          std::make_pair(42, "L1_DoubleMu0_SQ_OS"),          std::make_pair(142, "L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8"),          std::make_pair(59, "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4"),          std::make_pair(55, "L1_DoubleMu0er1p5_SQ"),          std::make_pair(56, "L1_DoubleMu0er1p5_SQ_OS"),          std::make_pair(58, "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"),          std::make_pair(57, "L1_DoubleMu0er1p5_SQ_dR_Max1p4"),          std::make_pair(54, "L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4"),          std::make_pair(53, "L1_DoubleMu0er2p0_SQ_dR_Max1p4"),          std::make_pair(45, "L1_DoubleMu10_SQ"),          std::make_pair(51, "L1_DoubleMu18er2p1"),          std::make_pair(112, "L1_DoubleMu3_OS_DoubleEG7p5Upsilon"),          std::make_pair(145, "L1_DoubleMu3_SQ_ETMHF50_HTT60er"),          std::make_pair(147, "L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5"),          std::make_pair(146, "L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5"),          std::make_pair(148, "L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5"),          std::make_pair(150, "L1_DoubleMu3_SQ_HTT220er"),          std::make_pair(151, "L1_DoubleMu3_SQ_HTT240er"),          std::make_pair(152, "L1_DoubleMu3_SQ_HTT260er"),          std::make_pair(143, "L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8"),          std::make_pair(109, "L1_DoubleMu4_SQ_EG9er2p5"),          std::make_pair(60, "L1_DoubleMu4_SQ_OS"),          std::make_pair(61, "L1_DoubleMu4_SQ_OS_dR_Max1p2"),          std::make_pair(62, "L1_DoubleMu4p5_SQ_OS"),          std::make_pair(63, "L1_DoubleMu4p5_SQ_OS_dR_Max1p2"),          std::make_pair(64, "L1_DoubleMu4p5er2p0_SQ_OS"),          std::make_pair(65, "L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18"),          std::make_pair(113, "L1_DoubleMu5Upsilon_OS_DoubleEG3"),          std::make_pair(110, "L1_DoubleMu5_SQ_EG9er2p5"),          std::make_pair(44, "L1_DoubleMu9_SQ"),          std::make_pair(46, "L1_DoubleMu_12_5"),          std::make_pair(47, "L1_DoubleMu_15_5_SQ"),          std::make_pair(48, "L1_DoubleMu_15_7"),          std::make_pair(50, "L1_DoubleMu_15_7_Mass_Min1"),          std::make_pair(49, "L1_DoubleMu_15_7_SQ"),          std::make_pair(273, "L1_DoubleTau70er2p1"),          std::make_pair(416, "L1_ETM120"),          std::make_pair(417, "L1_ETM150"),          std::make_pair(421, "L1_ETMHF100"),          std::make_pair(429, "L1_ETMHF100_HTT60er"),          std::make_pair(422, "L1_ETMHF110"),          std::make_pair(430, "L1_ETMHF110_HTT60er"),          std::make_pair(444, "L1_ETMHF110_HTT60er_NotSecondBunchInTrain"),          std::make_pair(423, "L1_ETMHF120"),          std::make_pair(431, "L1_ETMHF120_HTT60er"),          std::make_pair(443, "L1_ETMHF120_NotSecondBunchInTrain"),          std::make_pair(424, "L1_ETMHF130"),          std::make_pair(432, "L1_ETMHF130_HTT60er"),          std::make_pair(425, "L1_ETMHF140"),          std::make_pair(426, "L1_ETMHF150"),          std::make_pair(428, "L1_ETMHF90_HTT60er"),          std::make_pair(410, "L1_ETT1200"),          std::make_pair(411, "L1_ETT1600"),          std::make_pair(412, "L1_ETT2000"),          std::make_pair(477, "L1_FirstBunchAfterTrain"),          std::make_pair(472, "L1_FirstBunchBeforeTrain"),          std::make_pair(473, "L1_FirstBunchInTrain"),          std::make_pair(480, "L1_FirstCollisionInOrbit"),          std::make_pair(479, "L1_FirstCollisionInTrain"),          std::make_pair(500, "L1_HCAL_LaserMon_Trig"),          std::make_pair(501, "L1_HCAL_LaserMon_Veto"),          std::make_pair(398, "L1_HTT120er"),          std::make_pair(399, "L1_HTT160er"),          std::make_pair(400, "L1_HTT200er"),          std::make_pair(401, "L1_HTT255er"),          std::make_pair(402, "L1_HTT280er"),          std::make_pair(384, "L1_HTT280er_QuadJet_70_55_40_35_er2p4"),          std::make_pair(403, "L1_HTT320er"),          std::make_pair(385, "L1_HTT320er_QuadJet_70_55_40_40_er2p4"),          std::make_pair(386, "L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3"),          std::make_pair(387, "L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3"),          std::make_pair(404, "L1_HTT360er"),          std::make_pair(405, "L1_HTT400er"),          std::make_pair(406, "L1_HTT450er"),          std::make_pair(197, "L1_IsoEG32er2p5_Mt40"),          std::make_pair(198, "L1_IsoEG32er2p5_Mt44"),          std::make_pair(199, "L1_IsoEG32er2p5_Mt48"),          std::make_pair(294, "L1_IsoTau40er2p1_ETMHF100"),          std::make_pair(295, "L1_IsoTau40er2p1_ETMHF110"),          std::make_pair(296, "L1_IsoTau40er2p1_ETMHF120"),          std::make_pair(293, "L1_IsoTau40er2p1_ETMHF90"),          std::make_pair(471, "L1_IsolatedBunch"),          std::make_pair(476, "L1_LastBunchInTrain"),          std::make_pair(478, "L1_LastCollisionInTrain"),          std::make_pair(257, "L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3"),          std::make_pair(259, "L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3"),          std::make_pair(238, "L1_LooseIsoEG24er2p1_HTT100er"),          std::make_pair(258, "L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3"),          std::make_pair(239, "L1_LooseIsoEG26er2p1_HTT100er"),          std::make_pair(234, "L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3"),          std::make_pair(240, "L1_LooseIsoEG28er2p1_HTT100er"),          std::make_pair(235, "L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3"),          std::make_pair(241, "L1_LooseIsoEG30er2p1_HTT100er"),          std::make_pair(236, "L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3"),          std::make_pair(461, "L1_MinimumBiasHF0_AND_BptxAND"),          std::make_pair(134, "L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6"),          std::make_pair(136, "L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6"),          std::make_pair(135, "L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6"),          std::make_pair(281, "L1_Mu18er2p1_Tau24er2p1"),          std::make_pair(282, "L1_Mu18er2p1_Tau26er2p1"),          std::make_pair(98, "L1_Mu20_EG10er2p5"),          std::make_pair(284, "L1_Mu22er2p1_IsoTau32er2p1"),          std::make_pair(285, "L1_Mu22er2p1_IsoTau34er2p1"),          std::make_pair(286, "L1_Mu22er2p1_IsoTau36er2p1"),          std::make_pair(287, "L1_Mu22er2p1_IsoTau40er2p1"),          std::make_pair(289, "L1_Mu22er2p1_Tau70er2p1"),          std::make_pair(126, "L1_Mu3_Jet120er2p5_dR_Max0p4"),          std::make_pair(125, "L1_Mu3_Jet120er2p5_dR_Max0p8"),          std::make_pair(121, "L1_Mu3_Jet16er2p5_dR_Max0p4"),          std::make_pair(119, "L1_Mu3_Jet30er2p5"),          std::make_pair(122, "L1_Mu3_Jet35er2p5_dR_Max0p4"),          std::make_pair(123, "L1_Mu3_Jet60er2p5_dR_Max0p4"),          std::make_pair(124, "L1_Mu3_Jet80er2p5_dR_Max0p4"),          std::make_pair(128, "L1_Mu3er1p5_Jet100er2p5_ETMHF40"),          std::make_pair(129, "L1_Mu3er1p5_Jet100er2p5_ETMHF50"),          std::make_pair(96, "L1_Mu5_EG23er2p5"),          std::make_pair(100, "L1_Mu5_LooseIsoEG20er2p5"),          std::make_pair(104, "L1_Mu6_DoubleEG10er2p5"),          std::make_pair(105, "L1_Mu6_DoubleEG12er2p5"),          std::make_pair(106, "L1_Mu6_DoubleEG15er2p5"),          std::make_pair(107, "L1_Mu6_DoubleEG17er2p5"),          std::make_pair(131, "L1_Mu6_HTT240er"),          std::make_pair(132, "L1_Mu6_HTT250er"),          std::make_pair(97, "L1_Mu7_EG23er2p5"),          std::make_pair(101, "L1_Mu7_LooseIsoEG20er2p5"),          std::make_pair(102, "L1_Mu7_LooseIsoEG23er2p5"),          std::make_pair(463, "L1_NotBptxOR"),          std::make_pair(298, "L1_QuadJet36er2p5_IsoTau52er2p1"),          std::make_pair(382, "L1_QuadJet60er2p5"),          std::make_pair(376, "L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0"),          std::make_pair(89, "L1_QuadMu0"),          std::make_pair(88, "L1_QuadMu0_OQ"),          std::make_pair(90, "L1_QuadMu0_SQ"),          std::make_pair(474, "L1_SecondBunchInTrain"),          std::make_pair(475, "L1_SecondLastBunchInTrain"),          std::make_pair(164, "L1_SingleEG10er2p5"),          std::make_pair(165, "L1_SingleEG15er2p5"),          std::make_pair(166, "L1_SingleEG26er2p5"),          std::make_pair(167, "L1_SingleEG34er2p5"),          std::make_pair(168, "L1_SingleEG36er2p5"),          std::make_pair(169, "L1_SingleEG38er2p5"),          std::make_pair(170, "L1_SingleEG40er2p5"),          std::make_pair(171, "L1_SingleEG42er2p5"),          std::make_pair(172, "L1_SingleEG45er2p5"),          std::make_pair(173, "L1_SingleEG50"),          std::make_pair(174, "L1_SingleEG60"),          std::make_pair(163, "L1_SingleEG8er2p5"),          std::make_pair(184, "L1_SingleIsoEG24er1p5"),          std::make_pair(183, "L1_SingleIsoEG24er2p1"),          std::make_pair(187, "L1_SingleIsoEG26er1p5"),          std::make_pair(186, "L1_SingleIsoEG26er2p1"),          std::make_pair(185, "L1_SingleIsoEG26er2p5"),          std::make_pair(190, "L1_SingleIsoEG28er1p5"),          std::make_pair(189, "L1_SingleIsoEG28er2p1"),          std::make_pair(188, "L1_SingleIsoEG28er2p5"),          std::make_pair(192, "L1_SingleIsoEG30er2p1"),          std::make_pair(191, "L1_SingleIsoEG30er2p5"),          std::make_pair(194, "L1_SingleIsoEG32er2p1"),          std::make_pair(193, "L1_SingleIsoEG32er2p5"),          std::make_pair(195, "L1_SingleIsoEG34er2p5"),          std::make_pair(330, "L1_SingleJet10erHE"),          std::make_pair(312, "L1_SingleJet120"),          std::make_pair(327, "L1_SingleJet120_FWD3p0"),          std::make_pair(319, "L1_SingleJet120er2p5"),          std::make_pair(331, "L1_SingleJet12erHE"),          std::make_pair(320, "L1_SingleJet140er2p5"),          std::make_pair(333, "L1_SingleJet140er2p5_ETMHF80"),          std::make_pair(334, "L1_SingleJet140er2p5_ETMHF90"),          std::make_pair(321, "L1_SingleJet160er2p5"),          std::make_pair(313, "L1_SingleJet180"),          std::make_pair(322, "L1_SingleJet180er2p5"),          std::make_pair(314, "L1_SingleJet200"),          std::make_pair(450, "L1_SingleJet20er2p5_NotBptxOR"),          std::make_pair(451, "L1_SingleJet20er2p5_NotBptxOR_3BX"),          std::make_pair(309, "L1_SingleJet35"),          std::make_pair(324, "L1_SingleJet35_FWD3p0"),          std::make_pair(316, "L1_SingleJet35er2p5"),          std::make_pair(452, "L1_SingleJet43er2p5_NotBptxOR_3BX"),          std::make_pair(454, "L1_SingleJet46er2p5_NotBptxOR_3BX"),          std::make_pair(310, "L1_SingleJet60"),          std::make_pair(325, "L1_SingleJet60_FWD3p0"),          std::make_pair(317, "L1_SingleJet60er2p5"),          std::make_pair(329, "L1_SingleJet8erHE"),          std::make_pair(311, "L1_SingleJet90"),          std::make_pair(326, "L1_SingleJet90_FWD3p0"),          std::make_pair(318, "L1_SingleJet90er2p5"),          std::make_pair(180, "L1_SingleLooseIsoEG28er1p5"),          std::make_pair(181, "L1_SingleLooseIsoEG30er1p5"),          std::make_pair(6, "L1_SingleMu0_BMTF"),          std::make_pair(5, "L1_SingleMu0_DQ"),          std::make_pair(8, "L1_SingleMu0_EMTF"),          std::make_pair(7, "L1_SingleMu0_OMTF"),          std::make_pair(29, "L1_SingleMu10er1p5"),          std::make_pair(13, "L1_SingleMu12_DQ_BMTF"),          std::make_pair(15, "L1_SingleMu12_DQ_EMTF"),          std::make_pair(14, "L1_SingleMu12_DQ_OMTF"),          std::make_pair(30, "L1_SingleMu12er1p5"),          std::make_pair(31, "L1_SingleMu14er1p5"),          std::make_pair(16, "L1_SingleMu15_DQ"),          std::make_pair(32, "L1_SingleMu16er1p5"),          std::make_pair(17, "L1_SingleMu18"),          std::make_pair(33, "L1_SingleMu18er1p5"),          std::make_pair(18, "L1_SingleMu20"),          std::make_pair(19, "L1_SingleMu22"),          std::make_pair(20, "L1_SingleMu22_BMTF"),          std::make_pair(22, "L1_SingleMu22_EMTF"),          std::make_pair(21, "L1_SingleMu22_OMTF"),          std::make_pair(23, "L1_SingleMu25"),          std::make_pair(9, "L1_SingleMu3"),          std::make_pair(10, "L1_SingleMu5"),          std::make_pair(25, "L1_SingleMu6er1p5"),          std::make_pair(12, "L1_SingleMu7"),          std::make_pair(11, "L1_SingleMu7_DQ"),          std::make_pair(26, "L1_SingleMu7er1p5"),          std::make_pair(27, "L1_SingleMu8er1p5"),          std::make_pair(28, "L1_SingleMu9er1p5"),          std::make_pair(0, "L1_SingleMuCosmics"),          std::make_pair(1, "L1_SingleMuCosmics_BMTF"),          std::make_pair(3, "L1_SingleMuCosmics_EMTF"),          std::make_pair(2, "L1_SingleMuCosmics_OMTF"),          std::make_pair(4, "L1_SingleMuOpen"),          std::make_pair(446, "L1_SingleMuOpen_NotBptxOR"),          std::make_pair(448, "L1_SingleMuOpen_er1p1_NotBptxOR_3BX"),          std::make_pair(447, "L1_SingleMuOpen_er1p4_NotBptxOR_3BX"),          std::make_pair(270, "L1_SingleTau120er2p1"),          std::make_pair(271, "L1_SingleTau130er2p1"),          std::make_pair(503, "L1_TOTEM_1"),          std::make_pair(504, "L1_TOTEM_2"),          std::make_pair(505, "L1_TOTEM_3"),          std::make_pair(506, "L1_TOTEM_4"),          std::make_pair(228, "L1_TripleEG16er2p5"),          std::make_pair(224, "L1_TripleEG_16_12_8_er2p5"),          std::make_pair(225, "L1_TripleEG_16_15_8_er2p5"),          std::make_pair(226, "L1_TripleEG_18_17_8_er2p5"),          std::make_pair(227, "L1_TripleEG_18_18_12_er2p5"),          std::make_pair(373, "L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5"),          std::make_pair(374, "L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5"),          std::make_pair(372, "L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5"),          std::make_pair(72, "L1_TripleMu0"),          std::make_pair(71, "L1_TripleMu0_OQ"),          std::make_pair(73, "L1_TripleMu0_SQ"),          std::make_pair(74, "L1_TripleMu3"),          std::make_pair(75, "L1_TripleMu3_SQ"),          std::make_pair(76, "L1_TripleMu_5SQ_3SQ_0OQ"),          std::make_pair(85, "L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9"),          std::make_pair(86, "L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"),          std::make_pair(78, "L1_TripleMu_5_3_3"),          std::make_pair(79, "L1_TripleMu_5_3_3_SQ"),          std::make_pair(77, "L1_TripleMu_5_3p5_2p5"),          std::make_pair(83, "L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17"),          std::make_pair(82, "L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17"),          std::make_pair(84, "L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17"),          std::make_pair(80, "L1_TripleMu_5_5_3"),          std::make_pair(469, "L1_UnpairedBunchBptxMinus"),          std::make_pair(468, "L1_UnpairedBunchBptxPlus"),          std::make_pair(459, "L1_ZeroBias"),          std::make_pair(460, "L1_ZeroBias_copy")      };

  static const std::map<int, std::string> Id2Name(id2name, id2name + sizeof(id2name) / sizeof(id2name[0]));
  const std::map<int, std::string>::const_iterator rc = Id2Name.find(index);
  std::string name;
  if (rc != Id2Name.end()) name = rc->second;
  return name;
}


int getIdFromName(const std::string& name)
{
  static const std::pair<std::string, int> name2id[] = {
          std::make_pair("L1_AlwaysTrue", 458),          std::make_pair("L1_BPTX_AND_Ref1_VME", 486),          std::make_pair("L1_BPTX_AND_Ref3_VME", 487),          std::make_pair("L1_BPTX_AND_Ref4_VME", 488),          std::make_pair("L1_BPTX_BeamGas_B1_VME", 491),          std::make_pair("L1_BPTX_BeamGas_B2_VME", 492),          std::make_pair("L1_BPTX_BeamGas_Ref1_VME", 489),          std::make_pair("L1_BPTX_BeamGas_Ref2_VME", 490),          std::make_pair("L1_BPTX_NotOR_VME", 482),          std::make_pair("L1_BPTX_OR_Ref3_VME", 483),          std::make_pair("L1_BPTX_OR_Ref4_VME", 484),          std::make_pair("L1_BPTX_RefAND_VME", 485),          std::make_pair("L1_BptxMinus", 467),          std::make_pair("L1_BptxOR", 464),          std::make_pair("L1_BptxPlus", 466),          std::make_pair("L1_BptxXOR", 465),          std::make_pair("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", 494),          std::make_pair("L1_DoubleEG8er2p5_HTT260er", 247),          std::make_pair("L1_DoubleEG8er2p5_HTT280er", 248),          std::make_pair("L1_DoubleEG8er2p5_HTT300er", 249),          std::make_pair("L1_DoubleEG8er2p5_HTT320er", 250),          std::make_pair("L1_DoubleEG8er2p5_HTT340er", 251),          std::make_pair("L1_DoubleEG_15_10_er2p5", 205),          std::make_pair("L1_DoubleEG_20_10_er2p5", 206),          std::make_pair("L1_DoubleEG_22_10_er2p5", 207),          std::make_pair("L1_DoubleEG_25_12_er2p5", 208),          std::make_pair("L1_DoubleEG_25_14_er2p5", 209),          std::make_pair("L1_DoubleEG_27_14_er2p5", 210),          std::make_pair("L1_DoubleEG_LooseIso20_10_er2p5", 212),          std::make_pair("L1_DoubleEG_LooseIso22_10_er2p5", 213),          std::make_pair("L1_DoubleEG_LooseIso22_12_er2p5", 214),          std::make_pair("L1_DoubleEG_LooseIso25_12_er2p5", 215),          std::make_pair("L1_DoubleIsoTau32er2p1", 275),          std::make_pair("L1_DoubleIsoTau34er2p1", 276),          std::make_pair("L1_DoubleIsoTau36er2p1", 277),          std::make_pair("L1_DoubleJet100er2p3_dEta_Max1p6", 345),          std::make_pair("L1_DoubleJet100er2p5", 341),          std::make_pair("L1_DoubleJet112er2p3_dEta_Max1p6", 346),          std::make_pair("L1_DoubleJet120er2p5", 342),          std::make_pair("L1_DoubleJet150er2p5", 343),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5", 348),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5", 349),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", 350),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", 351),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", 352),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", 353),          std::make_pair("L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp", 363),          std::make_pair("L1_DoubleJet40er2p5", 340),          std::make_pair("L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", 356),          std::make_pair("L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", 357),          std::make_pair("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", 358),          std::make_pair("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28", 360),          std::make_pair("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620", 359),          std::make_pair("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28", 361),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ", 366),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp", 364),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_Mu8", 365),          std::make_pair("L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", 355),          std::make_pair("L1_DoubleLooseIsoEG22er2p1", 217),          std::make_pair("L1_DoubleLooseIsoEG24er2p1", 218),          std::make_pair("L1_DoubleMu0", 40),          std::make_pair("L1_DoubleMu0_Mass_Min1", 43),          std::make_pair("L1_DoubleMu0_OQ", 39),          std::make_pair("L1_DoubleMu0_SQ", 41),          std::make_pair("L1_DoubleMu0_SQ_OS", 42),          std::make_pair("L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8", 142),          std::make_pair("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", 59),          std::make_pair("L1_DoubleMu0er1p5_SQ", 55),          std::make_pair("L1_DoubleMu0er1p5_SQ_OS", 56),          std::make_pair("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", 58),          std::make_pair("L1_DoubleMu0er1p5_SQ_dR_Max1p4", 57),          std::make_pair("L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", 54),          std::make_pair("L1_DoubleMu0er2p0_SQ_dR_Max1p4", 53),          std::make_pair("L1_DoubleMu10_SQ", 45),          std::make_pair("L1_DoubleMu18er2p1", 51),          std::make_pair("L1_DoubleMu3_OS_DoubleEG7p5Upsilon", 112),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_HTT60er", 145),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5", 147),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5", 146),          std::make_pair("L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5", 148),          std::make_pair("L1_DoubleMu3_SQ_HTT220er", 150),          std::make_pair("L1_DoubleMu3_SQ_HTT240er", 151),          std::make_pair("L1_DoubleMu3_SQ_HTT260er", 152),          std::make_pair("L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8", 143),          std::make_pair("L1_DoubleMu4_SQ_EG9er2p5", 109),          std::make_pair("L1_DoubleMu4_SQ_OS", 60),          std::make_pair("L1_DoubleMu4_SQ_OS_dR_Max1p2", 61),          std::make_pair("L1_DoubleMu4p5_SQ_OS", 62),          std::make_pair("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", 63),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS", 64),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", 65),          std::make_pair("L1_DoubleMu5Upsilon_OS_DoubleEG3", 113),          std::make_pair("L1_DoubleMu5_SQ_EG9er2p5", 110),          std::make_pair("L1_DoubleMu9_SQ", 44),          std::make_pair("L1_DoubleMu_12_5", 46),          std::make_pair("L1_DoubleMu_15_5_SQ", 47),          std::make_pair("L1_DoubleMu_15_7", 48),          std::make_pair("L1_DoubleMu_15_7_Mass_Min1", 50),          std::make_pair("L1_DoubleMu_15_7_SQ", 49),          std::make_pair("L1_DoubleTau70er2p1", 273),          std::make_pair("L1_ETM120", 416),          std::make_pair("L1_ETM150", 417),          std::make_pair("L1_ETMHF100", 421),          std::make_pair("L1_ETMHF100_HTT60er", 429),          std::make_pair("L1_ETMHF110", 422),          std::make_pair("L1_ETMHF110_HTT60er", 430),          std::make_pair("L1_ETMHF110_HTT60er_NotSecondBunchInTrain", 444),          std::make_pair("L1_ETMHF120", 423),          std::make_pair("L1_ETMHF120_HTT60er", 431),          std::make_pair("L1_ETMHF120_NotSecondBunchInTrain", 443),          std::make_pair("L1_ETMHF130", 424),          std::make_pair("L1_ETMHF130_HTT60er", 432),          std::make_pair("L1_ETMHF140", 425),          std::make_pair("L1_ETMHF150", 426),          std::make_pair("L1_ETMHF90_HTT60er", 428),          std::make_pair("L1_ETT1200", 410),          std::make_pair("L1_ETT1600", 411),          std::make_pair("L1_ETT2000", 412),          std::make_pair("L1_FirstBunchAfterTrain", 477),          std::make_pair("L1_FirstBunchBeforeTrain", 472),          std::make_pair("L1_FirstBunchInTrain", 473),          std::make_pair("L1_FirstCollisionInOrbit", 480),          std::make_pair("L1_FirstCollisionInTrain", 479),          std::make_pair("L1_HCAL_LaserMon_Trig", 500),          std::make_pair("L1_HCAL_LaserMon_Veto", 501),          std::make_pair("L1_HTT120er", 398),          std::make_pair("L1_HTT160er", 399),          std::make_pair("L1_HTT200er", 400),          std::make_pair("L1_HTT255er", 401),          std::make_pair("L1_HTT280er", 402),          std::make_pair("L1_HTT280er_QuadJet_70_55_40_35_er2p4", 384),          std::make_pair("L1_HTT320er", 403),          std::make_pair("L1_HTT320er_QuadJet_70_55_40_40_er2p4", 385),          std::make_pair("L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", 386),          std::make_pair("L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", 387),          std::make_pair("L1_HTT360er", 404),          std::make_pair("L1_HTT400er", 405),          std::make_pair("L1_HTT450er", 406),          std::make_pair("L1_IsoEG32er2p5_Mt40", 197),          std::make_pair("L1_IsoEG32er2p5_Mt44", 198),          std::make_pair("L1_IsoEG32er2p5_Mt48", 199),          std::make_pair("L1_IsoTau40er2p1_ETMHF100", 294),          std::make_pair("L1_IsoTau40er2p1_ETMHF110", 295),          std::make_pair("L1_IsoTau40er2p1_ETMHF120", 296),          std::make_pair("L1_IsoTau40er2p1_ETMHF90", 293),          std::make_pair("L1_IsolatedBunch", 471),          std::make_pair("L1_LastBunchInTrain", 476),          std::make_pair("L1_LastCollisionInTrain", 478),          std::make_pair("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", 257),          std::make_pair("L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3", 259),          std::make_pair("L1_LooseIsoEG24er2p1_HTT100er", 238),          std::make_pair("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", 258),          std::make_pair("L1_LooseIsoEG26er2p1_HTT100er", 239),          std::make_pair("L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3", 234),          std::make_pair("L1_LooseIsoEG28er2p1_HTT100er", 240),          std::make_pair("L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3", 235),          std::make_pair("L1_LooseIsoEG30er2p1_HTT100er", 241),          std::make_pair("L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3", 236),          std::make_pair("L1_MinimumBiasHF0_AND_BptxAND", 461),          std::make_pair("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", 134),          std::make_pair("L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6", 136),          std::make_pair("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", 135),          std::make_pair("L1_Mu18er2p1_Tau24er2p1", 281),          std::make_pair("L1_Mu18er2p1_Tau26er2p1", 282),          std::make_pair("L1_Mu20_EG10er2p5", 98),          std::make_pair("L1_Mu22er2p1_IsoTau32er2p1", 284),          std::make_pair("L1_Mu22er2p1_IsoTau34er2p1", 285),          std::make_pair("L1_Mu22er2p1_IsoTau36er2p1", 286),          std::make_pair("L1_Mu22er2p1_IsoTau40er2p1", 287),          std::make_pair("L1_Mu22er2p1_Tau70er2p1", 289),          std::make_pair("L1_Mu3_Jet120er2p5_dR_Max0p4", 126),          std::make_pair("L1_Mu3_Jet120er2p5_dR_Max0p8", 125),          std::make_pair("L1_Mu3_Jet16er2p5_dR_Max0p4", 121),          std::make_pair("L1_Mu3_Jet30er2p5", 119),          std::make_pair("L1_Mu3_Jet35er2p5_dR_Max0p4", 122),          std::make_pair("L1_Mu3_Jet60er2p5_dR_Max0p4", 123),          std::make_pair("L1_Mu3_Jet80er2p5_dR_Max0p4", 124),          std::make_pair("L1_Mu3er1p5_Jet100er2p5_ETMHF40", 128),          std::make_pair("L1_Mu3er1p5_Jet100er2p5_ETMHF50", 129),          std::make_pair("L1_Mu5_EG23er2p5", 96),          std::make_pair("L1_Mu5_LooseIsoEG20er2p5", 100),          std::make_pair("L1_Mu6_DoubleEG10er2p5", 104),          std::make_pair("L1_Mu6_DoubleEG12er2p5", 105),          std::make_pair("L1_Mu6_DoubleEG15er2p5", 106),          std::make_pair("L1_Mu6_DoubleEG17er2p5", 107),          std::make_pair("L1_Mu6_HTT240er", 131),          std::make_pair("L1_Mu6_HTT250er", 132),          std::make_pair("L1_Mu7_EG23er2p5", 97),          std::make_pair("L1_Mu7_LooseIsoEG20er2p5", 101),          std::make_pair("L1_Mu7_LooseIsoEG23er2p5", 102),          std::make_pair("L1_NotBptxOR", 463),          std::make_pair("L1_QuadJet36er2p5_IsoTau52er2p1", 298),          std::make_pair("L1_QuadJet60er2p5", 382),          std::make_pair("L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0", 376),          std::make_pair("L1_QuadMu0", 89),          std::make_pair("L1_QuadMu0_OQ", 88),          std::make_pair("L1_QuadMu0_SQ", 90),          std::make_pair("L1_SecondBunchInTrain", 474),          std::make_pair("L1_SecondLastBunchInTrain", 475),          std::make_pair("L1_SingleEG10er2p5", 164),          std::make_pair("L1_SingleEG15er2p5", 165),          std::make_pair("L1_SingleEG26er2p5", 166),          std::make_pair("L1_SingleEG34er2p5", 167),          std::make_pair("L1_SingleEG36er2p5", 168),          std::make_pair("L1_SingleEG38er2p5", 169),          std::make_pair("L1_SingleEG40er2p5", 170),          std::make_pair("L1_SingleEG42er2p5", 171),          std::make_pair("L1_SingleEG45er2p5", 172),          std::make_pair("L1_SingleEG50", 173),          std::make_pair("L1_SingleEG60", 174),          std::make_pair("L1_SingleEG8er2p5", 163),          std::make_pair("L1_SingleIsoEG24er1p5", 184),          std::make_pair("L1_SingleIsoEG24er2p1", 183),          std::make_pair("L1_SingleIsoEG26er1p5", 187),          std::make_pair("L1_SingleIsoEG26er2p1", 186),          std::make_pair("L1_SingleIsoEG26er2p5", 185),          std::make_pair("L1_SingleIsoEG28er1p5", 190),          std::make_pair("L1_SingleIsoEG28er2p1", 189),          std::make_pair("L1_SingleIsoEG28er2p5", 188),          std::make_pair("L1_SingleIsoEG30er2p1", 192),          std::make_pair("L1_SingleIsoEG30er2p5", 191),          std::make_pair("L1_SingleIsoEG32er2p1", 194),          std::make_pair("L1_SingleIsoEG32er2p5", 193),          std::make_pair("L1_SingleIsoEG34er2p5", 195),          std::make_pair("L1_SingleJet10erHE", 330),          std::make_pair("L1_SingleJet120", 312),          std::make_pair("L1_SingleJet120_FWD3p0", 327),          std::make_pair("L1_SingleJet120er2p5", 319),          std::make_pair("L1_SingleJet12erHE", 331),          std::make_pair("L1_SingleJet140er2p5", 320),          std::make_pair("L1_SingleJet140er2p5_ETMHF80", 333),          std::make_pair("L1_SingleJet140er2p5_ETMHF90", 334),          std::make_pair("L1_SingleJet160er2p5", 321),          std::make_pair("L1_SingleJet180", 313),          std::make_pair("L1_SingleJet180er2p5", 322),          std::make_pair("L1_SingleJet200", 314),          std::make_pair("L1_SingleJet20er2p5_NotBptxOR", 450),          std::make_pair("L1_SingleJet20er2p5_NotBptxOR_3BX", 451),          std::make_pair("L1_SingleJet35", 309),          std::make_pair("L1_SingleJet35_FWD3p0", 324),          std::make_pair("L1_SingleJet35er2p5", 316),          std::make_pair("L1_SingleJet43er2p5_NotBptxOR_3BX", 452),          std::make_pair("L1_SingleJet46er2p5_NotBptxOR_3BX", 454),          std::make_pair("L1_SingleJet60", 310),          std::make_pair("L1_SingleJet60_FWD3p0", 325),          std::make_pair("L1_SingleJet60er2p5", 317),          std::make_pair("L1_SingleJet8erHE", 329),          std::make_pair("L1_SingleJet90", 311),          std::make_pair("L1_SingleJet90_FWD3p0", 326),          std::make_pair("L1_SingleJet90er2p5", 318),          std::make_pair("L1_SingleLooseIsoEG28er1p5", 180),          std::make_pair("L1_SingleLooseIsoEG30er1p5", 181),          std::make_pair("L1_SingleMu0_BMTF", 6),          std::make_pair("L1_SingleMu0_DQ", 5),          std::make_pair("L1_SingleMu0_EMTF", 8),          std::make_pair("L1_SingleMu0_OMTF", 7),          std::make_pair("L1_SingleMu10er1p5", 29),          std::make_pair("L1_SingleMu12_DQ_BMTF", 13),          std::make_pair("L1_SingleMu12_DQ_EMTF", 15),          std::make_pair("L1_SingleMu12_DQ_OMTF", 14),          std::make_pair("L1_SingleMu12er1p5", 30),          std::make_pair("L1_SingleMu14er1p5", 31),          std::make_pair("L1_SingleMu15_DQ", 16),          std::make_pair("L1_SingleMu16er1p5", 32),          std::make_pair("L1_SingleMu18", 17),          std::make_pair("L1_SingleMu18er1p5", 33),          std::make_pair("L1_SingleMu20", 18),          std::make_pair("L1_SingleMu22", 19),          std::make_pair("L1_SingleMu22_BMTF", 20),          std::make_pair("L1_SingleMu22_EMTF", 22),          std::make_pair("L1_SingleMu22_OMTF", 21),          std::make_pair("L1_SingleMu25", 23),          std::make_pair("L1_SingleMu3", 9),          std::make_pair("L1_SingleMu5", 10),          std::make_pair("L1_SingleMu6er1p5", 25),          std::make_pair("L1_SingleMu7", 12),          std::make_pair("L1_SingleMu7_DQ", 11),          std::make_pair("L1_SingleMu7er1p5", 26),          std::make_pair("L1_SingleMu8er1p5", 27),          std::make_pair("L1_SingleMu9er1p5", 28),          std::make_pair("L1_SingleMuCosmics", 0),          std::make_pair("L1_SingleMuCosmics_BMTF", 1),          std::make_pair("L1_SingleMuCosmics_EMTF", 3),          std::make_pair("L1_SingleMuCosmics_OMTF", 2),          std::make_pair("L1_SingleMuOpen", 4),          std::make_pair("L1_SingleMuOpen_NotBptxOR", 446),          std::make_pair("L1_SingleMuOpen_er1p1_NotBptxOR_3BX", 448),          std::make_pair("L1_SingleMuOpen_er1p4_NotBptxOR_3BX", 447),          std::make_pair("L1_SingleTau120er2p1", 270),          std::make_pair("L1_SingleTau130er2p1", 271),          std::make_pair("L1_TOTEM_1", 503),          std::make_pair("L1_TOTEM_2", 504),          std::make_pair("L1_TOTEM_3", 505),          std::make_pair("L1_TOTEM_4", 506),          std::make_pair("L1_TripleEG16er2p5", 228),          std::make_pair("L1_TripleEG_16_12_8_er2p5", 224),          std::make_pair("L1_TripleEG_16_15_8_er2p5", 225),          std::make_pair("L1_TripleEG_18_17_8_er2p5", 226),          std::make_pair("L1_TripleEG_18_18_12_er2p5", 227),          std::make_pair("L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5", 373),          std::make_pair("L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5", 374),          std::make_pair("L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5", 372),          std::make_pair("L1_TripleMu0", 72),          std::make_pair("L1_TripleMu0_OQ", 71),          std::make_pair("L1_TripleMu0_SQ", 73),          std::make_pair("L1_TripleMu3", 74),          std::make_pair("L1_TripleMu3_SQ", 75),          std::make_pair("L1_TripleMu_5SQ_3SQ_0OQ", 76),          std::make_pair("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", 85),          std::make_pair("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", 86),          std::make_pair("L1_TripleMu_5_3_3", 78),          std::make_pair("L1_TripleMu_5_3_3_SQ", 79),          std::make_pair("L1_TripleMu_5_3p5_2p5", 77),          std::make_pair("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", 83),          std::make_pair("L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", 82),          std::make_pair("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", 84),          std::make_pair("L1_TripleMu_5_5_3", 80),          std::make_pair("L1_UnpairedBunchBptxMinus", 469),          std::make_pair("L1_UnpairedBunchBptxPlus", 468),          std::make_pair("L1_ZeroBias", 459),          std::make_pair("L1_ZeroBias_copy", 460)      };

  static const std::map<std::string, int> Name2Id(name2id, name2id + sizeof(name2id) / sizeof(name2id[0]));
  const std::map<std::string, int>::const_iterator rc = Name2Id.find(name);
  int id = -1;
  if (rc != Name2Id.end()) id = rc->second;
  return id;
}


AlgorithmFunction getFuncFromId(const int index)
{
  static const std::pair<int, AlgorithmFunction> id2func[] = {
          std::make_pair(458, &L1_AlwaysTrue),          std::make_pair(486, &L1_BPTX_AND_Ref1_VME),          std::make_pair(487, &L1_BPTX_AND_Ref3_VME),          std::make_pair(488, &L1_BPTX_AND_Ref4_VME),          std::make_pair(491, &L1_BPTX_BeamGas_B1_VME),          std::make_pair(492, &L1_BPTX_BeamGas_B2_VME),          std::make_pair(489, &L1_BPTX_BeamGas_Ref1_VME),          std::make_pair(490, &L1_BPTX_BeamGas_Ref2_VME),          std::make_pair(482, &L1_BPTX_NotOR_VME),          std::make_pair(483, &L1_BPTX_OR_Ref3_VME),          std::make_pair(484, &L1_BPTX_OR_Ref4_VME),          std::make_pair(485, &L1_BPTX_RefAND_VME),          std::make_pair(467, &L1_BptxMinus),          std::make_pair(464, &L1_BptxOR),          std::make_pair(466, &L1_BptxPlus),          std::make_pair(465, &L1_BptxXOR),          std::make_pair(494, &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142),          std::make_pair(247, &L1_DoubleEG8er2p5_HTT260er),          std::make_pair(248, &L1_DoubleEG8er2p5_HTT280er),          std::make_pair(249, &L1_DoubleEG8er2p5_HTT300er),          std::make_pair(250, &L1_DoubleEG8er2p5_HTT320er),          std::make_pair(251, &L1_DoubleEG8er2p5_HTT340er),          std::make_pair(205, &L1_DoubleEG_15_10_er2p5),          std::make_pair(206, &L1_DoubleEG_20_10_er2p5),          std::make_pair(207, &L1_DoubleEG_22_10_er2p5),          std::make_pair(208, &L1_DoubleEG_25_12_er2p5),          std::make_pair(209, &L1_DoubleEG_25_14_er2p5),          std::make_pair(210, &L1_DoubleEG_27_14_er2p5),          std::make_pair(212, &L1_DoubleEG_LooseIso20_10_er2p5),          std::make_pair(213, &L1_DoubleEG_LooseIso22_10_er2p5),          std::make_pair(214, &L1_DoubleEG_LooseIso22_12_er2p5),          std::make_pair(215, &L1_DoubleEG_LooseIso25_12_er2p5),          std::make_pair(275, &L1_DoubleIsoTau32er2p1),          std::make_pair(276, &L1_DoubleIsoTau34er2p1),          std::make_pair(277, &L1_DoubleIsoTau36er2p1),          std::make_pair(345, &L1_DoubleJet100er2p3_dEta_Max1p6),          std::make_pair(341, &L1_DoubleJet100er2p5),          std::make_pair(346, &L1_DoubleJet112er2p3_dEta_Max1p6),          std::make_pair(342, &L1_DoubleJet120er2p5),          std::make_pair(343, &L1_DoubleJet150er2p5),          std::make_pair(348, &L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5),          std::make_pair(349, &L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5),          std::make_pair(350, &L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5),          std::make_pair(351, &L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5),          std::make_pair(352, &L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5),          std::make_pair(353, &L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5),          std::make_pair(363, &L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp),          std::make_pair(340, &L1_DoubleJet40er2p5),          std::make_pair(356, &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620),          std::make_pair(357, &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620),          std::make_pair(358, &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620),          std::make_pair(360, &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28),          std::make_pair(359, &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620),          std::make_pair(361, &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28),          std::make_pair(366, &L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ),          std::make_pair(364, &L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp),          std::make_pair(365, &L1_DoubleJet_80_30_Mass_Min420_Mu8),          std::make_pair(355, &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620),          std::make_pair(217, &L1_DoubleLooseIsoEG22er2p1),          std::make_pair(218, &L1_DoubleLooseIsoEG24er2p1),          std::make_pair(40, &L1_DoubleMu0),          std::make_pair(43, &L1_DoubleMu0_Mass_Min1),          std::make_pair(39, &L1_DoubleMu0_OQ),          std::make_pair(41, &L1_DoubleMu0_SQ),          std::make_pair(42, &L1_DoubleMu0_SQ_OS),          std::make_pair(142, &L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8),          std::make_pair(59, &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4),          std::make_pair(55, &L1_DoubleMu0er1p5_SQ),          std::make_pair(56, &L1_DoubleMu0er1p5_SQ_OS),          std::make_pair(58, &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4),          std::make_pair(57, &L1_DoubleMu0er1p5_SQ_dR_Max1p4),          std::make_pair(54, &L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4),          std::make_pair(53, &L1_DoubleMu0er2p0_SQ_dR_Max1p4),          std::make_pair(45, &L1_DoubleMu10_SQ),          std::make_pair(51, &L1_DoubleMu18er2p1),          std::make_pair(112, &L1_DoubleMu3_OS_DoubleEG7p5Upsilon),          std::make_pair(145, &L1_DoubleMu3_SQ_ETMHF50_HTT60er),          std::make_pair(147, &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5),          std::make_pair(146, &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5),          std::make_pair(148, &L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5),          std::make_pair(150, &L1_DoubleMu3_SQ_HTT220er),          std::make_pair(151, &L1_DoubleMu3_SQ_HTT240er),          std::make_pair(152, &L1_DoubleMu3_SQ_HTT260er),          std::make_pair(143, &L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8),          std::make_pair(109, &L1_DoubleMu4_SQ_EG9er2p5),          std::make_pair(60, &L1_DoubleMu4_SQ_OS),          std::make_pair(61, &L1_DoubleMu4_SQ_OS_dR_Max1p2),          std::make_pair(62, &L1_DoubleMu4p5_SQ_OS),          std::make_pair(63, &L1_DoubleMu4p5_SQ_OS_dR_Max1p2),          std::make_pair(64, &L1_DoubleMu4p5er2p0_SQ_OS),          std::make_pair(65, &L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18),          std::make_pair(113, &L1_DoubleMu5Upsilon_OS_DoubleEG3),          std::make_pair(110, &L1_DoubleMu5_SQ_EG9er2p5),          std::make_pair(44, &L1_DoubleMu9_SQ),          std::make_pair(46, &L1_DoubleMu_12_5),          std::make_pair(47, &L1_DoubleMu_15_5_SQ),          std::make_pair(48, &L1_DoubleMu_15_7),          std::make_pair(50, &L1_DoubleMu_15_7_Mass_Min1),          std::make_pair(49, &L1_DoubleMu_15_7_SQ),          std::make_pair(273, &L1_DoubleTau70er2p1),          std::make_pair(416, &L1_ETM120),          std::make_pair(417, &L1_ETM150),          std::make_pair(421, &L1_ETMHF100),          std::make_pair(429, &L1_ETMHF100_HTT60er),          std::make_pair(422, &L1_ETMHF110),          std::make_pair(430, &L1_ETMHF110_HTT60er),          std::make_pair(444, &L1_ETMHF110_HTT60er_NotSecondBunchInTrain),          std::make_pair(423, &L1_ETMHF120),          std::make_pair(431, &L1_ETMHF120_HTT60er),          std::make_pair(443, &L1_ETMHF120_NotSecondBunchInTrain),          std::make_pair(424, &L1_ETMHF130),          std::make_pair(432, &L1_ETMHF130_HTT60er),          std::make_pair(425, &L1_ETMHF140),          std::make_pair(426, &L1_ETMHF150),          std::make_pair(428, &L1_ETMHF90_HTT60er),          std::make_pair(410, &L1_ETT1200),          std::make_pair(411, &L1_ETT1600),          std::make_pair(412, &L1_ETT2000),          std::make_pair(477, &L1_FirstBunchAfterTrain),          std::make_pair(472, &L1_FirstBunchBeforeTrain),          std::make_pair(473, &L1_FirstBunchInTrain),          std::make_pair(480, &L1_FirstCollisionInOrbit),          std::make_pair(479, &L1_FirstCollisionInTrain),          std::make_pair(500, &L1_HCAL_LaserMon_Trig),          std::make_pair(501, &L1_HCAL_LaserMon_Veto),          std::make_pair(398, &L1_HTT120er),          std::make_pair(399, &L1_HTT160er),          std::make_pair(400, &L1_HTT200er),          std::make_pair(401, &L1_HTT255er),          std::make_pair(402, &L1_HTT280er),          std::make_pair(384, &L1_HTT280er_QuadJet_70_55_40_35_er2p4),          std::make_pair(403, &L1_HTT320er),          std::make_pair(385, &L1_HTT320er_QuadJet_70_55_40_40_er2p4),          std::make_pair(386, &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3),          std::make_pair(387, &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3),          std::make_pair(404, &L1_HTT360er),          std::make_pair(405, &L1_HTT400er),          std::make_pair(406, &L1_HTT450er),          std::make_pair(197, &L1_IsoEG32er2p5_Mt40),          std::make_pair(198, &L1_IsoEG32er2p5_Mt44),          std::make_pair(199, &L1_IsoEG32er2p5_Mt48),          std::make_pair(294, &L1_IsoTau40er2p1_ETMHF100),          std::make_pair(295, &L1_IsoTau40er2p1_ETMHF110),          std::make_pair(296, &L1_IsoTau40er2p1_ETMHF120),          std::make_pair(293, &L1_IsoTau40er2p1_ETMHF90),          std::make_pair(471, &L1_IsolatedBunch),          std::make_pair(476, &L1_LastBunchInTrain),          std::make_pair(478, &L1_LastCollisionInTrain),          std::make_pair(257, &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3),          std::make_pair(259, &L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3),          std::make_pair(238, &L1_LooseIsoEG24er2p1_HTT100er),          std::make_pair(258, &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3),          std::make_pair(239, &L1_LooseIsoEG26er2p1_HTT100er),          std::make_pair(234, &L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair(240, &L1_LooseIsoEG28er2p1_HTT100er),          std::make_pair(235, &L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair(241, &L1_LooseIsoEG30er2p1_HTT100er),          std::make_pair(236, &L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair(461, &L1_MinimumBiasHF0_AND_BptxAND),          std::make_pair(134, &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6),          std::make_pair(136, &L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6),          std::make_pair(135, &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6),          std::make_pair(281, &L1_Mu18er2p1_Tau24er2p1),          std::make_pair(282, &L1_Mu18er2p1_Tau26er2p1),          std::make_pair(98, &L1_Mu20_EG10er2p5),          std::make_pair(284, &L1_Mu22er2p1_IsoTau32er2p1),          std::make_pair(285, &L1_Mu22er2p1_IsoTau34er2p1),          std::make_pair(286, &L1_Mu22er2p1_IsoTau36er2p1),          std::make_pair(287, &L1_Mu22er2p1_IsoTau40er2p1),          std::make_pair(289, &L1_Mu22er2p1_Tau70er2p1),          std::make_pair(126, &L1_Mu3_Jet120er2p5_dR_Max0p4),          std::make_pair(125, &L1_Mu3_Jet120er2p5_dR_Max0p8),          std::make_pair(121, &L1_Mu3_Jet16er2p5_dR_Max0p4),          std::make_pair(119, &L1_Mu3_Jet30er2p5),          std::make_pair(122, &L1_Mu3_Jet35er2p5_dR_Max0p4),          std::make_pair(123, &L1_Mu3_Jet60er2p5_dR_Max0p4),          std::make_pair(124, &L1_Mu3_Jet80er2p5_dR_Max0p4),          std::make_pair(128, &L1_Mu3er1p5_Jet100er2p5_ETMHF40),          std::make_pair(129, &L1_Mu3er1p5_Jet100er2p5_ETMHF50),          std::make_pair(96, &L1_Mu5_EG23er2p5),          std::make_pair(100, &L1_Mu5_LooseIsoEG20er2p5),          std::make_pair(104, &L1_Mu6_DoubleEG10er2p5),          std::make_pair(105, &L1_Mu6_DoubleEG12er2p5),          std::make_pair(106, &L1_Mu6_DoubleEG15er2p5),          std::make_pair(107, &L1_Mu6_DoubleEG17er2p5),          std::make_pair(131, &L1_Mu6_HTT240er),          std::make_pair(132, &L1_Mu6_HTT250er),          std::make_pair(97, &L1_Mu7_EG23er2p5),          std::make_pair(101, &L1_Mu7_LooseIsoEG20er2p5),          std::make_pair(102, &L1_Mu7_LooseIsoEG23er2p5),          std::make_pair(463, &L1_NotBptxOR),          std::make_pair(298, &L1_QuadJet36er2p5_IsoTau52er2p1),          std::make_pair(382, &L1_QuadJet60er2p5),          std::make_pair(376, &L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0),          std::make_pair(89, &L1_QuadMu0),          std::make_pair(88, &L1_QuadMu0_OQ),          std::make_pair(90, &L1_QuadMu0_SQ),          std::make_pair(474, &L1_SecondBunchInTrain),          std::make_pair(475, &L1_SecondLastBunchInTrain),          std::make_pair(164, &L1_SingleEG10er2p5),          std::make_pair(165, &L1_SingleEG15er2p5),          std::make_pair(166, &L1_SingleEG26er2p5),          std::make_pair(167, &L1_SingleEG34er2p5),          std::make_pair(168, &L1_SingleEG36er2p5),          std::make_pair(169, &L1_SingleEG38er2p5),          std::make_pair(170, &L1_SingleEG40er2p5),          std::make_pair(171, &L1_SingleEG42er2p5),          std::make_pair(172, &L1_SingleEG45er2p5),          std::make_pair(173, &L1_SingleEG50),          std::make_pair(174, &L1_SingleEG60),          std::make_pair(163, &L1_SingleEG8er2p5),          std::make_pair(184, &L1_SingleIsoEG24er1p5),          std::make_pair(183, &L1_SingleIsoEG24er2p1),          std::make_pair(187, &L1_SingleIsoEG26er1p5),          std::make_pair(186, &L1_SingleIsoEG26er2p1),          std::make_pair(185, &L1_SingleIsoEG26er2p5),          std::make_pair(190, &L1_SingleIsoEG28er1p5),          std::make_pair(189, &L1_SingleIsoEG28er2p1),          std::make_pair(188, &L1_SingleIsoEG28er2p5),          std::make_pair(192, &L1_SingleIsoEG30er2p1),          std::make_pair(191, &L1_SingleIsoEG30er2p5),          std::make_pair(194, &L1_SingleIsoEG32er2p1),          std::make_pair(193, &L1_SingleIsoEG32er2p5),          std::make_pair(195, &L1_SingleIsoEG34er2p5),          std::make_pair(330, &L1_SingleJet10erHE),          std::make_pair(312, &L1_SingleJet120),          std::make_pair(327, &L1_SingleJet120_FWD3p0),          std::make_pair(319, &L1_SingleJet120er2p5),          std::make_pair(331, &L1_SingleJet12erHE),          std::make_pair(320, &L1_SingleJet140er2p5),          std::make_pair(333, &L1_SingleJet140er2p5_ETMHF80),          std::make_pair(334, &L1_SingleJet140er2p5_ETMHF90),          std::make_pair(321, &L1_SingleJet160er2p5),          std::make_pair(313, &L1_SingleJet180),          std::make_pair(322, &L1_SingleJet180er2p5),          std::make_pair(314, &L1_SingleJet200),          std::make_pair(450, &L1_SingleJet20er2p5_NotBptxOR),          std::make_pair(451, &L1_SingleJet20er2p5_NotBptxOR_3BX),          std::make_pair(309, &L1_SingleJet35),          std::make_pair(324, &L1_SingleJet35_FWD3p0),          std::make_pair(316, &L1_SingleJet35er2p5),          std::make_pair(452, &L1_SingleJet43er2p5_NotBptxOR_3BX),          std::make_pair(454, &L1_SingleJet46er2p5_NotBptxOR_3BX),          std::make_pair(310, &L1_SingleJet60),          std::make_pair(325, &L1_SingleJet60_FWD3p0),          std::make_pair(317, &L1_SingleJet60er2p5),          std::make_pair(329, &L1_SingleJet8erHE),          std::make_pair(311, &L1_SingleJet90),          std::make_pair(326, &L1_SingleJet90_FWD3p0),          std::make_pair(318, &L1_SingleJet90er2p5),          std::make_pair(180, &L1_SingleLooseIsoEG28er1p5),          std::make_pair(181, &L1_SingleLooseIsoEG30er1p5),          std::make_pair(6, &L1_SingleMu0_BMTF),          std::make_pair(5, &L1_SingleMu0_DQ),          std::make_pair(8, &L1_SingleMu0_EMTF),          std::make_pair(7, &L1_SingleMu0_OMTF),          std::make_pair(29, &L1_SingleMu10er1p5),          std::make_pair(13, &L1_SingleMu12_DQ_BMTF),          std::make_pair(15, &L1_SingleMu12_DQ_EMTF),          std::make_pair(14, &L1_SingleMu12_DQ_OMTF),          std::make_pair(30, &L1_SingleMu12er1p5),          std::make_pair(31, &L1_SingleMu14er1p5),          std::make_pair(16, &L1_SingleMu15_DQ),          std::make_pair(32, &L1_SingleMu16er1p5),          std::make_pair(17, &L1_SingleMu18),          std::make_pair(33, &L1_SingleMu18er1p5),          std::make_pair(18, &L1_SingleMu20),          std::make_pair(19, &L1_SingleMu22),          std::make_pair(20, &L1_SingleMu22_BMTF),          std::make_pair(22, &L1_SingleMu22_EMTF),          std::make_pair(21, &L1_SingleMu22_OMTF),          std::make_pair(23, &L1_SingleMu25),          std::make_pair(9, &L1_SingleMu3),          std::make_pair(10, &L1_SingleMu5),          std::make_pair(25, &L1_SingleMu6er1p5),          std::make_pair(12, &L1_SingleMu7),          std::make_pair(11, &L1_SingleMu7_DQ),          std::make_pair(26, &L1_SingleMu7er1p5),          std::make_pair(27, &L1_SingleMu8er1p5),          std::make_pair(28, &L1_SingleMu9er1p5),          std::make_pair(0, &L1_SingleMuCosmics),          std::make_pair(1, &L1_SingleMuCosmics_BMTF),          std::make_pair(3, &L1_SingleMuCosmics_EMTF),          std::make_pair(2, &L1_SingleMuCosmics_OMTF),          std::make_pair(4, &L1_SingleMuOpen),          std::make_pair(446, &L1_SingleMuOpen_NotBptxOR),          std::make_pair(448, &L1_SingleMuOpen_er1p1_NotBptxOR_3BX),          std::make_pair(447, &L1_SingleMuOpen_er1p4_NotBptxOR_3BX),          std::make_pair(270, &L1_SingleTau120er2p1),          std::make_pair(271, &L1_SingleTau130er2p1),          std::make_pair(503, &L1_TOTEM_1),          std::make_pair(504, &L1_TOTEM_2),          std::make_pair(505, &L1_TOTEM_3),          std::make_pair(506, &L1_TOTEM_4),          std::make_pair(228, &L1_TripleEG16er2p5),          std::make_pair(224, &L1_TripleEG_16_12_8_er2p5),          std::make_pair(225, &L1_TripleEG_16_15_8_er2p5),          std::make_pair(226, &L1_TripleEG_18_17_8_er2p5),          std::make_pair(227, &L1_TripleEG_18_18_12_er2p5),          std::make_pair(373, &L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5),          std::make_pair(374, &L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5),          std::make_pair(372, &L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5),          std::make_pair(72, &L1_TripleMu0),          std::make_pair(71, &L1_TripleMu0_OQ),          std::make_pair(73, &L1_TripleMu0_SQ),          std::make_pair(74, &L1_TripleMu3),          std::make_pair(75, &L1_TripleMu3_SQ),          std::make_pair(76, &L1_TripleMu_5SQ_3SQ_0OQ),          std::make_pair(85, &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9),          std::make_pair(86, &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9),          std::make_pair(78, &L1_TripleMu_5_3_3),          std::make_pair(79, &L1_TripleMu_5_3_3_SQ),          std::make_pair(77, &L1_TripleMu_5_3p5_2p5),          std::make_pair(83, &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17),          std::make_pair(82, &L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17),          std::make_pair(84, &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17),          std::make_pair(80, &L1_TripleMu_5_5_3),          std::make_pair(469, &L1_UnpairedBunchBptxMinus),          std::make_pair(468, &L1_UnpairedBunchBptxPlus),          std::make_pair(459, &L1_ZeroBias),          std::make_pair(460, &L1_ZeroBias_copy)      };

  static const std::map<int, AlgorithmFunction> Id2Func(id2func, id2func + sizeof(id2func) / sizeof(id2func[0]));
  const std::map<int, AlgorithmFunction>::const_iterator rc = Id2Func.find(index);
  AlgorithmFunction fp = 0;
  if (rc != Id2Func.end()) fp = rc->second;
  return fp;
}


AlgorithmFunction getFuncFromName(const std::string& name)
{
  static const std::pair<std::string, AlgorithmFunction> name2func[] = {
          std::make_pair("L1_AlwaysTrue", &L1_AlwaysTrue),          std::make_pair("L1_BPTX_AND_Ref1_VME", &L1_BPTX_AND_Ref1_VME),          std::make_pair("L1_BPTX_AND_Ref3_VME", &L1_BPTX_AND_Ref3_VME),          std::make_pair("L1_BPTX_AND_Ref4_VME", &L1_BPTX_AND_Ref4_VME),          std::make_pair("L1_BPTX_BeamGas_B1_VME", &L1_BPTX_BeamGas_B1_VME),          std::make_pair("L1_BPTX_BeamGas_B2_VME", &L1_BPTX_BeamGas_B2_VME),          std::make_pair("L1_BPTX_BeamGas_Ref1_VME", &L1_BPTX_BeamGas_Ref1_VME),          std::make_pair("L1_BPTX_BeamGas_Ref2_VME", &L1_BPTX_BeamGas_Ref2_VME),          std::make_pair("L1_BPTX_NotOR_VME", &L1_BPTX_NotOR_VME),          std::make_pair("L1_BPTX_OR_Ref3_VME", &L1_BPTX_OR_Ref3_VME),          std::make_pair("L1_BPTX_OR_Ref4_VME", &L1_BPTX_OR_Ref4_VME),          std::make_pair("L1_BPTX_RefAND_VME", &L1_BPTX_RefAND_VME),          std::make_pair("L1_BptxMinus", &L1_BptxMinus),          std::make_pair("L1_BptxOR", &L1_BptxOR),          std::make_pair("L1_BptxPlus", &L1_BptxPlus),          std::make_pair("L1_BptxXOR", &L1_BptxXOR),          std::make_pair("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142),          std::make_pair("L1_DoubleEG8er2p5_HTT260er", &L1_DoubleEG8er2p5_HTT260er),          std::make_pair("L1_DoubleEG8er2p5_HTT280er", &L1_DoubleEG8er2p5_HTT280er),          std::make_pair("L1_DoubleEG8er2p5_HTT300er", &L1_DoubleEG8er2p5_HTT300er),          std::make_pair("L1_DoubleEG8er2p5_HTT320er", &L1_DoubleEG8er2p5_HTT320er),          std::make_pair("L1_DoubleEG8er2p5_HTT340er", &L1_DoubleEG8er2p5_HTT340er),          std::make_pair("L1_DoubleEG_15_10_er2p5", &L1_DoubleEG_15_10_er2p5),          std::make_pair("L1_DoubleEG_20_10_er2p5", &L1_DoubleEG_20_10_er2p5),          std::make_pair("L1_DoubleEG_22_10_er2p5", &L1_DoubleEG_22_10_er2p5),          std::make_pair("L1_DoubleEG_25_12_er2p5", &L1_DoubleEG_25_12_er2p5),          std::make_pair("L1_DoubleEG_25_14_er2p5", &L1_DoubleEG_25_14_er2p5),          std::make_pair("L1_DoubleEG_27_14_er2p5", &L1_DoubleEG_27_14_er2p5),          std::make_pair("L1_DoubleEG_LooseIso20_10_er2p5", &L1_DoubleEG_LooseIso20_10_er2p5),          std::make_pair("L1_DoubleEG_LooseIso22_10_er2p5", &L1_DoubleEG_LooseIso22_10_er2p5),          std::make_pair("L1_DoubleEG_LooseIso22_12_er2p5", &L1_DoubleEG_LooseIso22_12_er2p5),          std::make_pair("L1_DoubleEG_LooseIso25_12_er2p5", &L1_DoubleEG_LooseIso25_12_er2p5),          std::make_pair("L1_DoubleIsoTau32er2p1", &L1_DoubleIsoTau32er2p1),          std::make_pair("L1_DoubleIsoTau34er2p1", &L1_DoubleIsoTau34er2p1),          std::make_pair("L1_DoubleIsoTau36er2p1", &L1_DoubleIsoTau36er2p1),          std::make_pair("L1_DoubleJet100er2p3_dEta_Max1p6", &L1_DoubleJet100er2p3_dEta_Max1p6),          std::make_pair("L1_DoubleJet100er2p5", &L1_DoubleJet100er2p5),          std::make_pair("L1_DoubleJet112er2p3_dEta_Max1p6", &L1_DoubleJet112er2p3_dEta_Max1p6),          std::make_pair("L1_DoubleJet120er2p5", &L1_DoubleJet120er2p5),          std::make_pair("L1_DoubleJet150er2p5", &L1_DoubleJet150er2p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5),          std::make_pair("L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp", &L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp),          std::make_pair("L1_DoubleJet40er2p5", &L1_DoubleJet40er2p5),          std::make_pair("L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620),          std::make_pair("L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620),          std::make_pair("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620),          std::make_pair("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28),          std::make_pair("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620),          std::make_pair("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ", &L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp", &L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_Mu8", &L1_DoubleJet_80_30_Mass_Min420_Mu8),          std::make_pair("L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620),          std::make_pair("L1_DoubleLooseIsoEG22er2p1", &L1_DoubleLooseIsoEG22er2p1),          std::make_pair("L1_DoubleLooseIsoEG24er2p1", &L1_DoubleLooseIsoEG24er2p1),          std::make_pair("L1_DoubleMu0", &L1_DoubleMu0),          std::make_pair("L1_DoubleMu0_Mass_Min1", &L1_DoubleMu0_Mass_Min1),          std::make_pair("L1_DoubleMu0_OQ", &L1_DoubleMu0_OQ),          std::make_pair("L1_DoubleMu0_SQ", &L1_DoubleMu0_SQ),          std::make_pair("L1_DoubleMu0_SQ_OS", &L1_DoubleMu0_SQ_OS),          std::make_pair("L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8),          std::make_pair("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4),          std::make_pair("L1_DoubleMu0er1p5_SQ", &L1_DoubleMu0er1p5_SQ),          std::make_pair("L1_DoubleMu0er1p5_SQ_OS", &L1_DoubleMu0er1p5_SQ_OS),          std::make_pair("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4),          std::make_pair("L1_DoubleMu0er1p5_SQ_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_dR_Max1p4),          std::make_pair("L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4),          std::make_pair("L1_DoubleMu0er2p0_SQ_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_dR_Max1p4),          std::make_pair("L1_DoubleMu10_SQ", &L1_DoubleMu10_SQ),          std::make_pair("L1_DoubleMu18er2p1", &L1_DoubleMu18er2p1),          std::make_pair("L1_DoubleMu3_OS_DoubleEG7p5Upsilon", &L1_DoubleMu3_OS_DoubleEG7p5Upsilon),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_HTT60er", &L1_DoubleMu3_SQ_ETMHF50_HTT60er),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5),          std::make_pair("L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5),          std::make_pair("L1_DoubleMu3_SQ_HTT220er", &L1_DoubleMu3_SQ_HTT220er),          std::make_pair("L1_DoubleMu3_SQ_HTT240er", &L1_DoubleMu3_SQ_HTT240er),          std::make_pair("L1_DoubleMu3_SQ_HTT260er", &L1_DoubleMu3_SQ_HTT260er),          std::make_pair("L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8),          std::make_pair("L1_DoubleMu4_SQ_EG9er2p5", &L1_DoubleMu4_SQ_EG9er2p5),          std::make_pair("L1_DoubleMu4_SQ_OS", &L1_DoubleMu4_SQ_OS),          std::make_pair("L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2),          std::make_pair("L1_DoubleMu4p5_SQ_OS", &L1_DoubleMu4p5_SQ_OS),          std::make_pair("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS", &L1_DoubleMu4p5er2p0_SQ_OS),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", &L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18),          std::make_pair("L1_DoubleMu5Upsilon_OS_DoubleEG3", &L1_DoubleMu5Upsilon_OS_DoubleEG3),          std::make_pair("L1_DoubleMu5_SQ_EG9er2p5", &L1_DoubleMu5_SQ_EG9er2p5),          std::make_pair("L1_DoubleMu9_SQ", &L1_DoubleMu9_SQ),          std::make_pair("L1_DoubleMu_12_5", &L1_DoubleMu_12_5),          std::make_pair("L1_DoubleMu_15_5_SQ", &L1_DoubleMu_15_5_SQ),          std::make_pair("L1_DoubleMu_15_7", &L1_DoubleMu_15_7),          std::make_pair("L1_DoubleMu_15_7_Mass_Min1", &L1_DoubleMu_15_7_Mass_Min1),          std::make_pair("L1_DoubleMu_15_7_SQ", &L1_DoubleMu_15_7_SQ),          std::make_pair("L1_DoubleTau70er2p1", &L1_DoubleTau70er2p1),          std::make_pair("L1_ETM120", &L1_ETM120),          std::make_pair("L1_ETM150", &L1_ETM150),          std::make_pair("L1_ETMHF100", &L1_ETMHF100),          std::make_pair("L1_ETMHF100_HTT60er", &L1_ETMHF100_HTT60er),          std::make_pair("L1_ETMHF110", &L1_ETMHF110),          std::make_pair("L1_ETMHF110_HTT60er", &L1_ETMHF110_HTT60er),          std::make_pair("L1_ETMHF110_HTT60er_NotSecondBunchInTrain", &L1_ETMHF110_HTT60er_NotSecondBunchInTrain),          std::make_pair("L1_ETMHF120", &L1_ETMHF120),          std::make_pair("L1_ETMHF120_HTT60er", &L1_ETMHF120_HTT60er),          std::make_pair("L1_ETMHF120_NotSecondBunchInTrain", &L1_ETMHF120_NotSecondBunchInTrain),          std::make_pair("L1_ETMHF130", &L1_ETMHF130),          std::make_pair("L1_ETMHF130_HTT60er", &L1_ETMHF130_HTT60er),          std::make_pair("L1_ETMHF140", &L1_ETMHF140),          std::make_pair("L1_ETMHF150", &L1_ETMHF150),          std::make_pair("L1_ETMHF90_HTT60er", &L1_ETMHF90_HTT60er),          std::make_pair("L1_ETT1200", &L1_ETT1200),          std::make_pair("L1_ETT1600", &L1_ETT1600),          std::make_pair("L1_ETT2000", &L1_ETT2000),          std::make_pair("L1_FirstBunchAfterTrain", &L1_FirstBunchAfterTrain),          std::make_pair("L1_FirstBunchBeforeTrain", &L1_FirstBunchBeforeTrain),          std::make_pair("L1_FirstBunchInTrain", &L1_FirstBunchInTrain),          std::make_pair("L1_FirstCollisionInOrbit", &L1_FirstCollisionInOrbit),          std::make_pair("L1_FirstCollisionInTrain", &L1_FirstCollisionInTrain),          std::make_pair("L1_HCAL_LaserMon_Trig", &L1_HCAL_LaserMon_Trig),          std::make_pair("L1_HCAL_LaserMon_Veto", &L1_HCAL_LaserMon_Veto),          std::make_pair("L1_HTT120er", &L1_HTT120er),          std::make_pair("L1_HTT160er", &L1_HTT160er),          std::make_pair("L1_HTT200er", &L1_HTT200er),          std::make_pair("L1_HTT255er", &L1_HTT255er),          std::make_pair("L1_HTT280er", &L1_HTT280er),          std::make_pair("L1_HTT280er_QuadJet_70_55_40_35_er2p4", &L1_HTT280er_QuadJet_70_55_40_35_er2p4),          std::make_pair("L1_HTT320er", &L1_HTT320er),          std::make_pair("L1_HTT320er_QuadJet_70_55_40_40_er2p4", &L1_HTT320er_QuadJet_70_55_40_40_er2p4),          std::make_pair("L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3),          std::make_pair("L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3),          std::make_pair("L1_HTT360er", &L1_HTT360er),          std::make_pair("L1_HTT400er", &L1_HTT400er),          std::make_pair("L1_HTT450er", &L1_HTT450er),          std::make_pair("L1_IsoEG32er2p5_Mt40", &L1_IsoEG32er2p5_Mt40),          std::make_pair("L1_IsoEG32er2p5_Mt44", &L1_IsoEG32er2p5_Mt44),          std::make_pair("L1_IsoEG32er2p5_Mt48", &L1_IsoEG32er2p5_Mt48),          std::make_pair("L1_IsoTau40er2p1_ETMHF100", &L1_IsoTau40er2p1_ETMHF100),          std::make_pair("L1_IsoTau40er2p1_ETMHF110", &L1_IsoTau40er2p1_ETMHF110),          std::make_pair("L1_IsoTau40er2p1_ETMHF120", &L1_IsoTau40er2p1_ETMHF120),          std::make_pair("L1_IsoTau40er2p1_ETMHF90", &L1_IsoTau40er2p1_ETMHF90),          std::make_pair("L1_IsolatedBunch", &L1_IsolatedBunch),          std::make_pair("L1_LastBunchInTrain", &L1_LastBunchInTrain),          std::make_pair("L1_LastCollisionInTrain", &L1_LastCollisionInTrain),          std::make_pair("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3),          std::make_pair("L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3),          std::make_pair("L1_LooseIsoEG24er2p1_HTT100er", &L1_LooseIsoEG24er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3),          std::make_pair("L1_LooseIsoEG26er2p1_HTT100er", &L1_LooseIsoEG26er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair("L1_LooseIsoEG28er2p1_HTT100er", &L1_LooseIsoEG28er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair("L1_LooseIsoEG30er2p1_HTT100er", &L1_LooseIsoEG30er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair("L1_MinimumBiasHF0_AND_BptxAND", &L1_MinimumBiasHF0_AND_BptxAND),          std::make_pair("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6),          std::make_pair("L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6),          std::make_pair("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6),          std::make_pair("L1_Mu18er2p1_Tau24er2p1", &L1_Mu18er2p1_Tau24er2p1),          std::make_pair("L1_Mu18er2p1_Tau26er2p1", &L1_Mu18er2p1_Tau26er2p1),          std::make_pair("L1_Mu20_EG10er2p5", &L1_Mu20_EG10er2p5),          std::make_pair("L1_Mu22er2p1_IsoTau32er2p1", &L1_Mu22er2p1_IsoTau32er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau34er2p1", &L1_Mu22er2p1_IsoTau34er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau36er2p1", &L1_Mu22er2p1_IsoTau36er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau40er2p1", &L1_Mu22er2p1_IsoTau40er2p1),          std::make_pair("L1_Mu22er2p1_Tau70er2p1", &L1_Mu22er2p1_Tau70er2p1),          std::make_pair("L1_Mu3_Jet120er2p5_dR_Max0p4", &L1_Mu3_Jet120er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet120er2p5_dR_Max0p8", &L1_Mu3_Jet120er2p5_dR_Max0p8),          std::make_pair("L1_Mu3_Jet16er2p5_dR_Max0p4", &L1_Mu3_Jet16er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet30er2p5", &L1_Mu3_Jet30er2p5),          std::make_pair("L1_Mu3_Jet35er2p5_dR_Max0p4", &L1_Mu3_Jet35er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet60er2p5_dR_Max0p4", &L1_Mu3_Jet60er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet80er2p5_dR_Max0p4", &L1_Mu3_Jet80er2p5_dR_Max0p4),          std::make_pair("L1_Mu3er1p5_Jet100er2p5_ETMHF40", &L1_Mu3er1p5_Jet100er2p5_ETMHF40),          std::make_pair("L1_Mu3er1p5_Jet100er2p5_ETMHF50", &L1_Mu3er1p5_Jet100er2p5_ETMHF50),          std::make_pair("L1_Mu5_EG23er2p5", &L1_Mu5_EG23er2p5),          std::make_pair("L1_Mu5_LooseIsoEG20er2p5", &L1_Mu5_LooseIsoEG20er2p5),          std::make_pair("L1_Mu6_DoubleEG10er2p5", &L1_Mu6_DoubleEG10er2p5),          std::make_pair("L1_Mu6_DoubleEG12er2p5", &L1_Mu6_DoubleEG12er2p5),          std::make_pair("L1_Mu6_DoubleEG15er2p5", &L1_Mu6_DoubleEG15er2p5),          std::make_pair("L1_Mu6_DoubleEG17er2p5", &L1_Mu6_DoubleEG17er2p5),          std::make_pair("L1_Mu6_HTT240er", &L1_Mu6_HTT240er),          std::make_pair("L1_Mu6_HTT250er", &L1_Mu6_HTT250er),          std::make_pair("L1_Mu7_EG23er2p5", &L1_Mu7_EG23er2p5),          std::make_pair("L1_Mu7_LooseIsoEG20er2p5", &L1_Mu7_LooseIsoEG20er2p5),          std::make_pair("L1_Mu7_LooseIsoEG23er2p5", &L1_Mu7_LooseIsoEG23er2p5),          std::make_pair("L1_NotBptxOR", &L1_NotBptxOR),          std::make_pair("L1_QuadJet36er2p5_IsoTau52er2p1", &L1_QuadJet36er2p5_IsoTau52er2p1),          std::make_pair("L1_QuadJet60er2p5", &L1_QuadJet60er2p5),          std::make_pair("L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0", &L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0),          std::make_pair("L1_QuadMu0", &L1_QuadMu0),          std::make_pair("L1_QuadMu0_OQ", &L1_QuadMu0_OQ),          std::make_pair("L1_QuadMu0_SQ", &L1_QuadMu0_SQ),          std::make_pair("L1_SecondBunchInTrain", &L1_SecondBunchInTrain),          std::make_pair("L1_SecondLastBunchInTrain", &L1_SecondLastBunchInTrain),          std::make_pair("L1_SingleEG10er2p5", &L1_SingleEG10er2p5),          std::make_pair("L1_SingleEG15er2p5", &L1_SingleEG15er2p5),          std::make_pair("L1_SingleEG26er2p5", &L1_SingleEG26er2p5),          std::make_pair("L1_SingleEG34er2p5", &L1_SingleEG34er2p5),          std::make_pair("L1_SingleEG36er2p5", &L1_SingleEG36er2p5),          std::make_pair("L1_SingleEG38er2p5", &L1_SingleEG38er2p5),          std::make_pair("L1_SingleEG40er2p5", &L1_SingleEG40er2p5),          std::make_pair("L1_SingleEG42er2p5", &L1_SingleEG42er2p5),          std::make_pair("L1_SingleEG45er2p5", &L1_SingleEG45er2p5),          std::make_pair("L1_SingleEG50", &L1_SingleEG50),          std::make_pair("L1_SingleEG60", &L1_SingleEG60),          std::make_pair("L1_SingleEG8er2p5", &L1_SingleEG8er2p5),          std::make_pair("L1_SingleIsoEG24er1p5", &L1_SingleIsoEG24er1p5),          std::make_pair("L1_SingleIsoEG24er2p1", &L1_SingleIsoEG24er2p1),          std::make_pair("L1_SingleIsoEG26er1p5", &L1_SingleIsoEG26er1p5),          std::make_pair("L1_SingleIsoEG26er2p1", &L1_SingleIsoEG26er2p1),          std::make_pair("L1_SingleIsoEG26er2p5", &L1_SingleIsoEG26er2p5),          std::make_pair("L1_SingleIsoEG28er1p5", &L1_SingleIsoEG28er1p5),          std::make_pair("L1_SingleIsoEG28er2p1", &L1_SingleIsoEG28er2p1),          std::make_pair("L1_SingleIsoEG28er2p5", &L1_SingleIsoEG28er2p5),          std::make_pair("L1_SingleIsoEG30er2p1", &L1_SingleIsoEG30er2p1),          std::make_pair("L1_SingleIsoEG30er2p5", &L1_SingleIsoEG30er2p5),          std::make_pair("L1_SingleIsoEG32er2p1", &L1_SingleIsoEG32er2p1),          std::make_pair("L1_SingleIsoEG32er2p5", &L1_SingleIsoEG32er2p5),          std::make_pair("L1_SingleIsoEG34er2p5", &L1_SingleIsoEG34er2p5),          std::make_pair("L1_SingleJet10erHE", &L1_SingleJet10erHE),          std::make_pair("L1_SingleJet120", &L1_SingleJet120),          std::make_pair("L1_SingleJet120_FWD3p0", &L1_SingleJet120_FWD3p0),          std::make_pair("L1_SingleJet120er2p5", &L1_SingleJet120er2p5),          std::make_pair("L1_SingleJet12erHE", &L1_SingleJet12erHE),          std::make_pair("L1_SingleJet140er2p5", &L1_SingleJet140er2p5),          std::make_pair("L1_SingleJet140er2p5_ETMHF80", &L1_SingleJet140er2p5_ETMHF80),          std::make_pair("L1_SingleJet140er2p5_ETMHF90", &L1_SingleJet140er2p5_ETMHF90),          std::make_pair("L1_SingleJet160er2p5", &L1_SingleJet160er2p5),          std::make_pair("L1_SingleJet180", &L1_SingleJet180),          std::make_pair("L1_SingleJet180er2p5", &L1_SingleJet180er2p5),          std::make_pair("L1_SingleJet200", &L1_SingleJet200),          std::make_pair("L1_SingleJet20er2p5_NotBptxOR", &L1_SingleJet20er2p5_NotBptxOR),          std::make_pair("L1_SingleJet20er2p5_NotBptxOR_3BX", &L1_SingleJet20er2p5_NotBptxOR_3BX),          std::make_pair("L1_SingleJet35", &L1_SingleJet35),          std::make_pair("L1_SingleJet35_FWD3p0", &L1_SingleJet35_FWD3p0),          std::make_pair("L1_SingleJet35er2p5", &L1_SingleJet35er2p5),          std::make_pair("L1_SingleJet43er2p5_NotBptxOR_3BX", &L1_SingleJet43er2p5_NotBptxOR_3BX),          std::make_pair("L1_SingleJet46er2p5_NotBptxOR_3BX", &L1_SingleJet46er2p5_NotBptxOR_3BX),          std::make_pair("L1_SingleJet60", &L1_SingleJet60),          std::make_pair("L1_SingleJet60_FWD3p0", &L1_SingleJet60_FWD3p0),          std::make_pair("L1_SingleJet60er2p5", &L1_SingleJet60er2p5),          std::make_pair("L1_SingleJet8erHE", &L1_SingleJet8erHE),          std::make_pair("L1_SingleJet90", &L1_SingleJet90),          std::make_pair("L1_SingleJet90_FWD3p0", &L1_SingleJet90_FWD3p0),          std::make_pair("L1_SingleJet90er2p5", &L1_SingleJet90er2p5),          std::make_pair("L1_SingleLooseIsoEG28er1p5", &L1_SingleLooseIsoEG28er1p5),          std::make_pair("L1_SingleLooseIsoEG30er1p5", &L1_SingleLooseIsoEG30er1p5),          std::make_pair("L1_SingleMu0_BMTF", &L1_SingleMu0_BMTF),          std::make_pair("L1_SingleMu0_DQ", &L1_SingleMu0_DQ),          std::make_pair("L1_SingleMu0_EMTF", &L1_SingleMu0_EMTF),          std::make_pair("L1_SingleMu0_OMTF", &L1_SingleMu0_OMTF),          std::make_pair("L1_SingleMu10er1p5", &L1_SingleMu10er1p5),          std::make_pair("L1_SingleMu12_DQ_BMTF", &L1_SingleMu12_DQ_BMTF),          std::make_pair("L1_SingleMu12_DQ_EMTF", &L1_SingleMu12_DQ_EMTF),          std::make_pair("L1_SingleMu12_DQ_OMTF", &L1_SingleMu12_DQ_OMTF),          std::make_pair("L1_SingleMu12er1p5", &L1_SingleMu12er1p5),          std::make_pair("L1_SingleMu14er1p5", &L1_SingleMu14er1p5),          std::make_pair("L1_SingleMu15_DQ", &L1_SingleMu15_DQ),          std::make_pair("L1_SingleMu16er1p5", &L1_SingleMu16er1p5),          std::make_pair("L1_SingleMu18", &L1_SingleMu18),          std::make_pair("L1_SingleMu18er1p5", &L1_SingleMu18er1p5),          std::make_pair("L1_SingleMu20", &L1_SingleMu20),          std::make_pair("L1_SingleMu22", &L1_SingleMu22),          std::make_pair("L1_SingleMu22_BMTF", &L1_SingleMu22_BMTF),          std::make_pair("L1_SingleMu22_EMTF", &L1_SingleMu22_EMTF),          std::make_pair("L1_SingleMu22_OMTF", &L1_SingleMu22_OMTF),          std::make_pair("L1_SingleMu25", &L1_SingleMu25),          std::make_pair("L1_SingleMu3", &L1_SingleMu3),          std::make_pair("L1_SingleMu5", &L1_SingleMu5),          std::make_pair("L1_SingleMu6er1p5", &L1_SingleMu6er1p5),          std::make_pair("L1_SingleMu7", &L1_SingleMu7),          std::make_pair("L1_SingleMu7_DQ", &L1_SingleMu7_DQ),          std::make_pair("L1_SingleMu7er1p5", &L1_SingleMu7er1p5),          std::make_pair("L1_SingleMu8er1p5", &L1_SingleMu8er1p5),          std::make_pair("L1_SingleMu9er1p5", &L1_SingleMu9er1p5),          std::make_pair("L1_SingleMuCosmics", &L1_SingleMuCosmics),          std::make_pair("L1_SingleMuCosmics_BMTF", &L1_SingleMuCosmics_BMTF),          std::make_pair("L1_SingleMuCosmics_EMTF", &L1_SingleMuCosmics_EMTF),          std::make_pair("L1_SingleMuCosmics_OMTF", &L1_SingleMuCosmics_OMTF),          std::make_pair("L1_SingleMuOpen", &L1_SingleMuOpen),          std::make_pair("L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR),          std::make_pair("L1_SingleMuOpen_er1p1_NotBptxOR_3BX", &L1_SingleMuOpen_er1p1_NotBptxOR_3BX),          std::make_pair("L1_SingleMuOpen_er1p4_NotBptxOR_3BX", &L1_SingleMuOpen_er1p4_NotBptxOR_3BX),          std::make_pair("L1_SingleTau120er2p1", &L1_SingleTau120er2p1),          std::make_pair("L1_SingleTau130er2p1", &L1_SingleTau130er2p1),          std::make_pair("L1_TOTEM_1", &L1_TOTEM_1),          std::make_pair("L1_TOTEM_2", &L1_TOTEM_2),          std::make_pair("L1_TOTEM_3", &L1_TOTEM_3),          std::make_pair("L1_TOTEM_4", &L1_TOTEM_4),          std::make_pair("L1_TripleEG16er2p5", &L1_TripleEG16er2p5),          std::make_pair("L1_TripleEG_16_12_8_er2p5", &L1_TripleEG_16_12_8_er2p5),          std::make_pair("L1_TripleEG_16_15_8_er2p5", &L1_TripleEG_16_15_8_er2p5),          std::make_pair("L1_TripleEG_18_17_8_er2p5", &L1_TripleEG_18_17_8_er2p5),          std::make_pair("L1_TripleEG_18_18_12_er2p5", &L1_TripleEG_18_18_12_er2p5),          std::make_pair("L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5", &L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5),          std::make_pair("L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5", &L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5),          std::make_pair("L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5", &L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5),          std::make_pair("L1_TripleMu0", &L1_TripleMu0),          std::make_pair("L1_TripleMu0_OQ", &L1_TripleMu0_OQ),          std::make_pair("L1_TripleMu0_SQ", &L1_TripleMu0_SQ),          std::make_pair("L1_TripleMu3", &L1_TripleMu3),          std::make_pair("L1_TripleMu3_SQ", &L1_TripleMu3_SQ),          std::make_pair("L1_TripleMu_5SQ_3SQ_0OQ", &L1_TripleMu_5SQ_3SQ_0OQ),          std::make_pair("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9),          std::make_pair("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9),          std::make_pair("L1_TripleMu_5_3_3", &L1_TripleMu_5_3_3),          std::make_pair("L1_TripleMu_5_3_3_SQ", &L1_TripleMu_5_3_3_SQ),          std::make_pair("L1_TripleMu_5_3p5_2p5", &L1_TripleMu_5_3p5_2p5),          std::make_pair("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17),          std::make_pair("L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17),          std::make_pair("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17),          std::make_pair("L1_TripleMu_5_5_3", &L1_TripleMu_5_5_3),          std::make_pair("L1_UnpairedBunchBptxMinus", &L1_UnpairedBunchBptxMinus),          std::make_pair("L1_UnpairedBunchBptxPlus", &L1_UnpairedBunchBptxPlus),          std::make_pair("L1_ZeroBias", &L1_ZeroBias),          std::make_pair("L1_ZeroBias_copy", &L1_ZeroBias_copy)      };

  static const std::map<std::string, AlgorithmFunction> Name2Func(name2func, name2func + sizeof(name2func) / sizeof(name2func[0]));
  const std::map<std::string, AlgorithmFunction>::const_iterator rc = Name2Func.find(name);
  AlgorithmFunction fp = 0;
  if (rc != Name2Func.end()) fp = rc->second;
  if (fp == 0)
  {
    std::stringstream ss;
    ss << "fat> algorithm '" << name << "' is not defined in L1Menu_Collisions2018_v1_0_0\n";
    throw std::runtime_error(ss.str());
  }
  return fp;
}


bool addFuncFromName(std::map<std::string, std::function<bool()>> &L1SeedFun,
                     L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade,
                     L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  static const std::pair<std::string, AlgorithmFunction> name2func[] = {
          std::make_pair("L1_AlwaysTrue", &L1_AlwaysTrue),          std::make_pair("L1_BPTX_AND_Ref1_VME", &L1_BPTX_AND_Ref1_VME),          std::make_pair("L1_BPTX_AND_Ref3_VME", &L1_BPTX_AND_Ref3_VME),          std::make_pair("L1_BPTX_AND_Ref4_VME", &L1_BPTX_AND_Ref4_VME),          std::make_pair("L1_BPTX_BeamGas_B1_VME", &L1_BPTX_BeamGas_B1_VME),          std::make_pair("L1_BPTX_BeamGas_B2_VME", &L1_BPTX_BeamGas_B2_VME),          std::make_pair("L1_BPTX_BeamGas_Ref1_VME", &L1_BPTX_BeamGas_Ref1_VME),          std::make_pair("L1_BPTX_BeamGas_Ref2_VME", &L1_BPTX_BeamGas_Ref2_VME),          std::make_pair("L1_BPTX_NotOR_VME", &L1_BPTX_NotOR_VME),          std::make_pair("L1_BPTX_OR_Ref3_VME", &L1_BPTX_OR_Ref3_VME),          std::make_pair("L1_BPTX_OR_Ref4_VME", &L1_BPTX_OR_Ref4_VME),          std::make_pair("L1_BPTX_RefAND_VME", &L1_BPTX_RefAND_VME),          std::make_pair("L1_BptxMinus", &L1_BptxMinus),          std::make_pair("L1_BptxOR", &L1_BptxOR),          std::make_pair("L1_BptxPlus", &L1_BptxPlus),          std::make_pair("L1_BptxXOR", &L1_BptxXOR),          std::make_pair("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142),          std::make_pair("L1_DoubleEG8er2p5_HTT260er", &L1_DoubleEG8er2p5_HTT260er),          std::make_pair("L1_DoubleEG8er2p5_HTT280er", &L1_DoubleEG8er2p5_HTT280er),          std::make_pair("L1_DoubleEG8er2p5_HTT300er", &L1_DoubleEG8er2p5_HTT300er),          std::make_pair("L1_DoubleEG8er2p5_HTT320er", &L1_DoubleEG8er2p5_HTT320er),          std::make_pair("L1_DoubleEG8er2p5_HTT340er", &L1_DoubleEG8er2p5_HTT340er),          std::make_pair("L1_DoubleEG_15_10_er2p5", &L1_DoubleEG_15_10_er2p5),          std::make_pair("L1_DoubleEG_20_10_er2p5", &L1_DoubleEG_20_10_er2p5),          std::make_pair("L1_DoubleEG_22_10_er2p5", &L1_DoubleEG_22_10_er2p5),          std::make_pair("L1_DoubleEG_25_12_er2p5", &L1_DoubleEG_25_12_er2p5),          std::make_pair("L1_DoubleEG_25_14_er2p5", &L1_DoubleEG_25_14_er2p5),          std::make_pair("L1_DoubleEG_27_14_er2p5", &L1_DoubleEG_27_14_er2p5),          std::make_pair("L1_DoubleEG_LooseIso20_10_er2p5", &L1_DoubleEG_LooseIso20_10_er2p5),          std::make_pair("L1_DoubleEG_LooseIso22_10_er2p5", &L1_DoubleEG_LooseIso22_10_er2p5),          std::make_pair("L1_DoubleEG_LooseIso22_12_er2p5", &L1_DoubleEG_LooseIso22_12_er2p5),          std::make_pair("L1_DoubleEG_LooseIso25_12_er2p5", &L1_DoubleEG_LooseIso25_12_er2p5),          std::make_pair("L1_DoubleIsoTau32er2p1", &L1_DoubleIsoTau32er2p1),          std::make_pair("L1_DoubleIsoTau34er2p1", &L1_DoubleIsoTau34er2p1),          std::make_pair("L1_DoubleIsoTau36er2p1", &L1_DoubleIsoTau36er2p1),          std::make_pair("L1_DoubleJet100er2p3_dEta_Max1p6", &L1_DoubleJet100er2p3_dEta_Max1p6),          std::make_pair("L1_DoubleJet100er2p5", &L1_DoubleJet100er2p5),          std::make_pair("L1_DoubleJet112er2p3_dEta_Max1p6", &L1_DoubleJet112er2p3_dEta_Max1p6),          std::make_pair("L1_DoubleJet120er2p5", &L1_DoubleJet120er2p5),          std::make_pair("L1_DoubleJet150er2p5", &L1_DoubleJet150er2p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5),          std::make_pair("L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp", &L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp),          std::make_pair("L1_DoubleJet40er2p5", &L1_DoubleJet40er2p5),          std::make_pair("L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620),          std::make_pair("L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620),          std::make_pair("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620),          std::make_pair("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28),          std::make_pair("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620),          std::make_pair("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ", &L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp", &L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_Mu8", &L1_DoubleJet_80_30_Mass_Min420_Mu8),          std::make_pair("L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620),          std::make_pair("L1_DoubleLooseIsoEG22er2p1", &L1_DoubleLooseIsoEG22er2p1),          std::make_pair("L1_DoubleLooseIsoEG24er2p1", &L1_DoubleLooseIsoEG24er2p1),          std::make_pair("L1_DoubleMu0", &L1_DoubleMu0),          std::make_pair("L1_DoubleMu0_Mass_Min1", &L1_DoubleMu0_Mass_Min1),          std::make_pair("L1_DoubleMu0_OQ", &L1_DoubleMu0_OQ),          std::make_pair("L1_DoubleMu0_SQ", &L1_DoubleMu0_SQ),          std::make_pair("L1_DoubleMu0_SQ_OS", &L1_DoubleMu0_SQ_OS),          std::make_pair("L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8),          std::make_pair("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4),          std::make_pair("L1_DoubleMu0er1p5_SQ", &L1_DoubleMu0er1p5_SQ),          std::make_pair("L1_DoubleMu0er1p5_SQ_OS", &L1_DoubleMu0er1p5_SQ_OS),          std::make_pair("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4),          std::make_pair("L1_DoubleMu0er1p5_SQ_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_dR_Max1p4),          std::make_pair("L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4),          std::make_pair("L1_DoubleMu0er2p0_SQ_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_dR_Max1p4),          std::make_pair("L1_DoubleMu10_SQ", &L1_DoubleMu10_SQ),          std::make_pair("L1_DoubleMu18er2p1", &L1_DoubleMu18er2p1),          std::make_pair("L1_DoubleMu3_OS_DoubleEG7p5Upsilon", &L1_DoubleMu3_OS_DoubleEG7p5Upsilon),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_HTT60er", &L1_DoubleMu3_SQ_ETMHF50_HTT60er),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5),          std::make_pair("L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5),          std::make_pair("L1_DoubleMu3_SQ_HTT220er", &L1_DoubleMu3_SQ_HTT220er),          std::make_pair("L1_DoubleMu3_SQ_HTT240er", &L1_DoubleMu3_SQ_HTT240er),          std::make_pair("L1_DoubleMu3_SQ_HTT260er", &L1_DoubleMu3_SQ_HTT260er),          std::make_pair("L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8),          std::make_pair("L1_DoubleMu4_SQ_EG9er2p5", &L1_DoubleMu4_SQ_EG9er2p5),          std::make_pair("L1_DoubleMu4_SQ_OS", &L1_DoubleMu4_SQ_OS),          std::make_pair("L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2),          std::make_pair("L1_DoubleMu4p5_SQ_OS", &L1_DoubleMu4p5_SQ_OS),          std::make_pair("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS", &L1_DoubleMu4p5er2p0_SQ_OS),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", &L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18),          std::make_pair("L1_DoubleMu5Upsilon_OS_DoubleEG3", &L1_DoubleMu5Upsilon_OS_DoubleEG3),          std::make_pair("L1_DoubleMu5_SQ_EG9er2p5", &L1_DoubleMu5_SQ_EG9er2p5),          std::make_pair("L1_DoubleMu9_SQ", &L1_DoubleMu9_SQ),          std::make_pair("L1_DoubleMu_12_5", &L1_DoubleMu_12_5),          std::make_pair("L1_DoubleMu_15_5_SQ", &L1_DoubleMu_15_5_SQ),          std::make_pair("L1_DoubleMu_15_7", &L1_DoubleMu_15_7),          std::make_pair("L1_DoubleMu_15_7_Mass_Min1", &L1_DoubleMu_15_7_Mass_Min1),          std::make_pair("L1_DoubleMu_15_7_SQ", &L1_DoubleMu_15_7_SQ),          std::make_pair("L1_DoubleTau70er2p1", &L1_DoubleTau70er2p1),          std::make_pair("L1_ETM120", &L1_ETM120),          std::make_pair("L1_ETM150", &L1_ETM150),          std::make_pair("L1_ETMHF100", &L1_ETMHF100),          std::make_pair("L1_ETMHF100_HTT60er", &L1_ETMHF100_HTT60er),          std::make_pair("L1_ETMHF110", &L1_ETMHF110),          std::make_pair("L1_ETMHF110_HTT60er", &L1_ETMHF110_HTT60er),          std::make_pair("L1_ETMHF110_HTT60er_NotSecondBunchInTrain", &L1_ETMHF110_HTT60er_NotSecondBunchInTrain),          std::make_pair("L1_ETMHF120", &L1_ETMHF120),          std::make_pair("L1_ETMHF120_HTT60er", &L1_ETMHF120_HTT60er),          std::make_pair("L1_ETMHF120_NotSecondBunchInTrain", &L1_ETMHF120_NotSecondBunchInTrain),          std::make_pair("L1_ETMHF130", &L1_ETMHF130),          std::make_pair("L1_ETMHF130_HTT60er", &L1_ETMHF130_HTT60er),          std::make_pair("L1_ETMHF140", &L1_ETMHF140),          std::make_pair("L1_ETMHF150", &L1_ETMHF150),          std::make_pair("L1_ETMHF90_HTT60er", &L1_ETMHF90_HTT60er),          std::make_pair("L1_ETT1200", &L1_ETT1200),          std::make_pair("L1_ETT1600", &L1_ETT1600),          std::make_pair("L1_ETT2000", &L1_ETT2000),          std::make_pair("L1_FirstBunchAfterTrain", &L1_FirstBunchAfterTrain),          std::make_pair("L1_FirstBunchBeforeTrain", &L1_FirstBunchBeforeTrain),          std::make_pair("L1_FirstBunchInTrain", &L1_FirstBunchInTrain),          std::make_pair("L1_FirstCollisionInOrbit", &L1_FirstCollisionInOrbit),          std::make_pair("L1_FirstCollisionInTrain", &L1_FirstCollisionInTrain),          std::make_pair("L1_HCAL_LaserMon_Trig", &L1_HCAL_LaserMon_Trig),          std::make_pair("L1_HCAL_LaserMon_Veto", &L1_HCAL_LaserMon_Veto),          std::make_pair("L1_HTT120er", &L1_HTT120er),          std::make_pair("L1_HTT160er", &L1_HTT160er),          std::make_pair("L1_HTT200er", &L1_HTT200er),          std::make_pair("L1_HTT255er", &L1_HTT255er),          std::make_pair("L1_HTT280er", &L1_HTT280er),          std::make_pair("L1_HTT280er_QuadJet_70_55_40_35_er2p4", &L1_HTT280er_QuadJet_70_55_40_35_er2p4),          std::make_pair("L1_HTT320er", &L1_HTT320er),          std::make_pair("L1_HTT320er_QuadJet_70_55_40_40_er2p4", &L1_HTT320er_QuadJet_70_55_40_40_er2p4),          std::make_pair("L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3),          std::make_pair("L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3),          std::make_pair("L1_HTT360er", &L1_HTT360er),          std::make_pair("L1_HTT400er", &L1_HTT400er),          std::make_pair("L1_HTT450er", &L1_HTT450er),          std::make_pair("L1_IsoEG32er2p5_Mt40", &L1_IsoEG32er2p5_Mt40),          std::make_pair("L1_IsoEG32er2p5_Mt44", &L1_IsoEG32er2p5_Mt44),          std::make_pair("L1_IsoEG32er2p5_Mt48", &L1_IsoEG32er2p5_Mt48),          std::make_pair("L1_IsoTau40er2p1_ETMHF100", &L1_IsoTau40er2p1_ETMHF100),          std::make_pair("L1_IsoTau40er2p1_ETMHF110", &L1_IsoTau40er2p1_ETMHF110),          std::make_pair("L1_IsoTau40er2p1_ETMHF120", &L1_IsoTau40er2p1_ETMHF120),          std::make_pair("L1_IsoTau40er2p1_ETMHF90", &L1_IsoTau40er2p1_ETMHF90),          std::make_pair("L1_IsolatedBunch", &L1_IsolatedBunch),          std::make_pair("L1_LastBunchInTrain", &L1_LastBunchInTrain),          std::make_pair("L1_LastCollisionInTrain", &L1_LastCollisionInTrain),          std::make_pair("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3),          std::make_pair("L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3),          std::make_pair("L1_LooseIsoEG24er2p1_HTT100er", &L1_LooseIsoEG24er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3),          std::make_pair("L1_LooseIsoEG26er2p1_HTT100er", &L1_LooseIsoEG26er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair("L1_LooseIsoEG28er2p1_HTT100er", &L1_LooseIsoEG28er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair("L1_LooseIsoEG30er2p1_HTT100er", &L1_LooseIsoEG30er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair("L1_MinimumBiasHF0_AND_BptxAND", &L1_MinimumBiasHF0_AND_BptxAND),          std::make_pair("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6),          std::make_pair("L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6),          std::make_pair("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6),          std::make_pair("L1_Mu18er2p1_Tau24er2p1", &L1_Mu18er2p1_Tau24er2p1),          std::make_pair("L1_Mu18er2p1_Tau26er2p1", &L1_Mu18er2p1_Tau26er2p1),          std::make_pair("L1_Mu20_EG10er2p5", &L1_Mu20_EG10er2p5),          std::make_pair("L1_Mu22er2p1_IsoTau32er2p1", &L1_Mu22er2p1_IsoTau32er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau34er2p1", &L1_Mu22er2p1_IsoTau34er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau36er2p1", &L1_Mu22er2p1_IsoTau36er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau40er2p1", &L1_Mu22er2p1_IsoTau40er2p1),          std::make_pair("L1_Mu22er2p1_Tau70er2p1", &L1_Mu22er2p1_Tau70er2p1),          std::make_pair("L1_Mu3_Jet120er2p5_dR_Max0p4", &L1_Mu3_Jet120er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet120er2p5_dR_Max0p8", &L1_Mu3_Jet120er2p5_dR_Max0p8),          std::make_pair("L1_Mu3_Jet16er2p5_dR_Max0p4", &L1_Mu3_Jet16er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet30er2p5", &L1_Mu3_Jet30er2p5),          std::make_pair("L1_Mu3_Jet35er2p5_dR_Max0p4", &L1_Mu3_Jet35er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet60er2p5_dR_Max0p4", &L1_Mu3_Jet60er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet80er2p5_dR_Max0p4", &L1_Mu3_Jet80er2p5_dR_Max0p4),          std::make_pair("L1_Mu3er1p5_Jet100er2p5_ETMHF40", &L1_Mu3er1p5_Jet100er2p5_ETMHF40),          std::make_pair("L1_Mu3er1p5_Jet100er2p5_ETMHF50", &L1_Mu3er1p5_Jet100er2p5_ETMHF50),          std::make_pair("L1_Mu5_EG23er2p5", &L1_Mu5_EG23er2p5),          std::make_pair("L1_Mu5_LooseIsoEG20er2p5", &L1_Mu5_LooseIsoEG20er2p5),          std::make_pair("L1_Mu6_DoubleEG10er2p5", &L1_Mu6_DoubleEG10er2p5),          std::make_pair("L1_Mu6_DoubleEG12er2p5", &L1_Mu6_DoubleEG12er2p5),          std::make_pair("L1_Mu6_DoubleEG15er2p5", &L1_Mu6_DoubleEG15er2p5),          std::make_pair("L1_Mu6_DoubleEG17er2p5", &L1_Mu6_DoubleEG17er2p5),          std::make_pair("L1_Mu6_HTT240er", &L1_Mu6_HTT240er),          std::make_pair("L1_Mu6_HTT250er", &L1_Mu6_HTT250er),          std::make_pair("L1_Mu7_EG23er2p5", &L1_Mu7_EG23er2p5),          std::make_pair("L1_Mu7_LooseIsoEG20er2p5", &L1_Mu7_LooseIsoEG20er2p5),          std::make_pair("L1_Mu7_LooseIsoEG23er2p5", &L1_Mu7_LooseIsoEG23er2p5),          std::make_pair("L1_NotBptxOR", &L1_NotBptxOR),          std::make_pair("L1_QuadJet36er2p5_IsoTau52er2p1", &L1_QuadJet36er2p5_IsoTau52er2p1),          std::make_pair("L1_QuadJet60er2p5", &L1_QuadJet60er2p5),          std::make_pair("L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0", &L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0),          std::make_pair("L1_QuadMu0", &L1_QuadMu0),          std::make_pair("L1_QuadMu0_OQ", &L1_QuadMu0_OQ),          std::make_pair("L1_QuadMu0_SQ", &L1_QuadMu0_SQ),          std::make_pair("L1_SecondBunchInTrain", &L1_SecondBunchInTrain),          std::make_pair("L1_SecondLastBunchInTrain", &L1_SecondLastBunchInTrain),          std::make_pair("L1_SingleEG10er2p5", &L1_SingleEG10er2p5),          std::make_pair("L1_SingleEG15er2p5", &L1_SingleEG15er2p5),          std::make_pair("L1_SingleEG26er2p5", &L1_SingleEG26er2p5),          std::make_pair("L1_SingleEG34er2p5", &L1_SingleEG34er2p5),          std::make_pair("L1_SingleEG36er2p5", &L1_SingleEG36er2p5),          std::make_pair("L1_SingleEG38er2p5", &L1_SingleEG38er2p5),          std::make_pair("L1_SingleEG40er2p5", &L1_SingleEG40er2p5),          std::make_pair("L1_SingleEG42er2p5", &L1_SingleEG42er2p5),          std::make_pair("L1_SingleEG45er2p5", &L1_SingleEG45er2p5),          std::make_pair("L1_SingleEG50", &L1_SingleEG50),          std::make_pair("L1_SingleEG60", &L1_SingleEG60),          std::make_pair("L1_SingleEG8er2p5", &L1_SingleEG8er2p5),          std::make_pair("L1_SingleIsoEG24er1p5", &L1_SingleIsoEG24er1p5),          std::make_pair("L1_SingleIsoEG24er2p1", &L1_SingleIsoEG24er2p1),          std::make_pair("L1_SingleIsoEG26er1p5", &L1_SingleIsoEG26er1p5),          std::make_pair("L1_SingleIsoEG26er2p1", &L1_SingleIsoEG26er2p1),          std::make_pair("L1_SingleIsoEG26er2p5", &L1_SingleIsoEG26er2p5),          std::make_pair("L1_SingleIsoEG28er1p5", &L1_SingleIsoEG28er1p5),          std::make_pair("L1_SingleIsoEG28er2p1", &L1_SingleIsoEG28er2p1),          std::make_pair("L1_SingleIsoEG28er2p5", &L1_SingleIsoEG28er2p5),          std::make_pair("L1_SingleIsoEG30er2p1", &L1_SingleIsoEG30er2p1),          std::make_pair("L1_SingleIsoEG30er2p5", &L1_SingleIsoEG30er2p5),          std::make_pair("L1_SingleIsoEG32er2p1", &L1_SingleIsoEG32er2p1),          std::make_pair("L1_SingleIsoEG32er2p5", &L1_SingleIsoEG32er2p5),          std::make_pair("L1_SingleIsoEG34er2p5", &L1_SingleIsoEG34er2p5),          std::make_pair("L1_SingleJet10erHE", &L1_SingleJet10erHE),          std::make_pair("L1_SingleJet120", &L1_SingleJet120),          std::make_pair("L1_SingleJet120_FWD3p0", &L1_SingleJet120_FWD3p0),          std::make_pair("L1_SingleJet120er2p5", &L1_SingleJet120er2p5),          std::make_pair("L1_SingleJet12erHE", &L1_SingleJet12erHE),          std::make_pair("L1_SingleJet140er2p5", &L1_SingleJet140er2p5),          std::make_pair("L1_SingleJet140er2p5_ETMHF80", &L1_SingleJet140er2p5_ETMHF80),          std::make_pair("L1_SingleJet140er2p5_ETMHF90", &L1_SingleJet140er2p5_ETMHF90),          std::make_pair("L1_SingleJet160er2p5", &L1_SingleJet160er2p5),          std::make_pair("L1_SingleJet180", &L1_SingleJet180),          std::make_pair("L1_SingleJet180er2p5", &L1_SingleJet180er2p5),          std::make_pair("L1_SingleJet200", &L1_SingleJet200),          std::make_pair("L1_SingleJet20er2p5_NotBptxOR", &L1_SingleJet20er2p5_NotBptxOR),          std::make_pair("L1_SingleJet20er2p5_NotBptxOR_3BX", &L1_SingleJet20er2p5_NotBptxOR_3BX),          std::make_pair("L1_SingleJet35", &L1_SingleJet35),          std::make_pair("L1_SingleJet35_FWD3p0", &L1_SingleJet35_FWD3p0),          std::make_pair("L1_SingleJet35er2p5", &L1_SingleJet35er2p5),          std::make_pair("L1_SingleJet43er2p5_NotBptxOR_3BX", &L1_SingleJet43er2p5_NotBptxOR_3BX),          std::make_pair("L1_SingleJet46er2p5_NotBptxOR_3BX", &L1_SingleJet46er2p5_NotBptxOR_3BX),          std::make_pair("L1_SingleJet60", &L1_SingleJet60),          std::make_pair("L1_SingleJet60_FWD3p0", &L1_SingleJet60_FWD3p0),          std::make_pair("L1_SingleJet60er2p5", &L1_SingleJet60er2p5),          std::make_pair("L1_SingleJet8erHE", &L1_SingleJet8erHE),          std::make_pair("L1_SingleJet90", &L1_SingleJet90),          std::make_pair("L1_SingleJet90_FWD3p0", &L1_SingleJet90_FWD3p0),          std::make_pair("L1_SingleJet90er2p5", &L1_SingleJet90er2p5),          std::make_pair("L1_SingleLooseIsoEG28er1p5", &L1_SingleLooseIsoEG28er1p5),          std::make_pair("L1_SingleLooseIsoEG30er1p5", &L1_SingleLooseIsoEG30er1p5),          std::make_pair("L1_SingleMu0_BMTF", &L1_SingleMu0_BMTF),          std::make_pair("L1_SingleMu0_DQ", &L1_SingleMu0_DQ),          std::make_pair("L1_SingleMu0_EMTF", &L1_SingleMu0_EMTF),          std::make_pair("L1_SingleMu0_OMTF", &L1_SingleMu0_OMTF),          std::make_pair("L1_SingleMu10er1p5", &L1_SingleMu10er1p5),          std::make_pair("L1_SingleMu12_DQ_BMTF", &L1_SingleMu12_DQ_BMTF),          std::make_pair("L1_SingleMu12_DQ_EMTF", &L1_SingleMu12_DQ_EMTF),          std::make_pair("L1_SingleMu12_DQ_OMTF", &L1_SingleMu12_DQ_OMTF),          std::make_pair("L1_SingleMu12er1p5", &L1_SingleMu12er1p5),          std::make_pair("L1_SingleMu14er1p5", &L1_SingleMu14er1p5),          std::make_pair("L1_SingleMu15_DQ", &L1_SingleMu15_DQ),          std::make_pair("L1_SingleMu16er1p5", &L1_SingleMu16er1p5),          std::make_pair("L1_SingleMu18", &L1_SingleMu18),          std::make_pair("L1_SingleMu18er1p5", &L1_SingleMu18er1p5),          std::make_pair("L1_SingleMu20", &L1_SingleMu20),          std::make_pair("L1_SingleMu22", &L1_SingleMu22),          std::make_pair("L1_SingleMu22_BMTF", &L1_SingleMu22_BMTF),          std::make_pair("L1_SingleMu22_EMTF", &L1_SingleMu22_EMTF),          std::make_pair("L1_SingleMu22_OMTF", &L1_SingleMu22_OMTF),          std::make_pair("L1_SingleMu25", &L1_SingleMu25),          std::make_pair("L1_SingleMu3", &L1_SingleMu3),          std::make_pair("L1_SingleMu5", &L1_SingleMu5),          std::make_pair("L1_SingleMu6er1p5", &L1_SingleMu6er1p5),          std::make_pair("L1_SingleMu7", &L1_SingleMu7),          std::make_pair("L1_SingleMu7_DQ", &L1_SingleMu7_DQ),          std::make_pair("L1_SingleMu7er1p5", &L1_SingleMu7er1p5),          std::make_pair("L1_SingleMu8er1p5", &L1_SingleMu8er1p5),          std::make_pair("L1_SingleMu9er1p5", &L1_SingleMu9er1p5),          std::make_pair("L1_SingleMuCosmics", &L1_SingleMuCosmics),          std::make_pair("L1_SingleMuCosmics_BMTF", &L1_SingleMuCosmics_BMTF),          std::make_pair("L1_SingleMuCosmics_EMTF", &L1_SingleMuCosmics_EMTF),          std::make_pair("L1_SingleMuCosmics_OMTF", &L1_SingleMuCosmics_OMTF),          std::make_pair("L1_SingleMuOpen", &L1_SingleMuOpen),          std::make_pair("L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR),          std::make_pair("L1_SingleMuOpen_er1p1_NotBptxOR_3BX", &L1_SingleMuOpen_er1p1_NotBptxOR_3BX),          std::make_pair("L1_SingleMuOpen_er1p4_NotBptxOR_3BX", &L1_SingleMuOpen_er1p4_NotBptxOR_3BX),          std::make_pair("L1_SingleTau120er2p1", &L1_SingleTau120er2p1),          std::make_pair("L1_SingleTau130er2p1", &L1_SingleTau130er2p1),          std::make_pair("L1_TOTEM_1", &L1_TOTEM_1),          std::make_pair("L1_TOTEM_2", &L1_TOTEM_2),          std::make_pair("L1_TOTEM_3", &L1_TOTEM_3),          std::make_pair("L1_TOTEM_4", &L1_TOTEM_4),          std::make_pair("L1_TripleEG16er2p5", &L1_TripleEG16er2p5),          std::make_pair("L1_TripleEG_16_12_8_er2p5", &L1_TripleEG_16_12_8_er2p5),          std::make_pair("L1_TripleEG_16_15_8_er2p5", &L1_TripleEG_16_15_8_er2p5),          std::make_pair("L1_TripleEG_18_17_8_er2p5", &L1_TripleEG_18_17_8_er2p5),          std::make_pair("L1_TripleEG_18_18_12_er2p5", &L1_TripleEG_18_18_12_er2p5),          std::make_pair("L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5", &L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5),          std::make_pair("L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5", &L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5),          std::make_pair("L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5", &L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5),          std::make_pair("L1_TripleMu0", &L1_TripleMu0),          std::make_pair("L1_TripleMu0_OQ", &L1_TripleMu0_OQ),          std::make_pair("L1_TripleMu0_SQ", &L1_TripleMu0_SQ),          std::make_pair("L1_TripleMu3", &L1_TripleMu3),          std::make_pair("L1_TripleMu3_SQ", &L1_TripleMu3_SQ),          std::make_pair("L1_TripleMu_5SQ_3SQ_0OQ", &L1_TripleMu_5SQ_3SQ_0OQ),          std::make_pair("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9),          std::make_pair("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9),          std::make_pair("L1_TripleMu_5_3_3", &L1_TripleMu_5_3_3),          std::make_pair("L1_TripleMu_5_3_3_SQ", &L1_TripleMu_5_3_3_SQ),          std::make_pair("L1_TripleMu_5_3p5_2p5", &L1_TripleMu_5_3p5_2p5),          std::make_pair("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17),          std::make_pair("L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17),          std::make_pair("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17),          std::make_pair("L1_TripleMu_5_5_3", &L1_TripleMu_5_5_3),          std::make_pair("L1_UnpairedBunchBptxMinus", &L1_UnpairedBunchBptxMinus),          std::make_pair("L1_UnpairedBunchBptxPlus", &L1_UnpairedBunchBptxPlus),          std::make_pair("L1_ZeroBias", &L1_ZeroBias),          std::make_pair("L1_ZeroBias_copy", &L1_ZeroBias_copy)      };

  for (auto pair : name2func)
  {
    L1SeedFun[pair.first] = std::bind(pair.second, upgrade, calo_tower);
  }

  return true;
}
// eof
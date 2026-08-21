
#ifndef MYSTRINGS_H
#define MYSTRINGS_H

#include <string>

namespace mystrings {
  string addSubscript(string base, string subscript);
  string addSuperscript(string base, string superscript);
  string formatHadronName(string hadron);
  string formatHadronDaughters(string hadron);
  string invMassOf(string system);
  string massOf(string particle);
  string getPtString(string subscript);
  string getZString(string subscript);
  string getRatioString(string num, string den);
  string getOneOverString(string s);
  string getdYdXString(string y, string x);
  string getdYdPtString(string y);
  string getdYdZString(string y);
  string getVarRangeString(string var, double high);
  string getVarRangeString(double low, string var, double high);
  string getPtJetRangeString(double ptmin, double ptmax, bool addUnits);
  string getPtV0RangeString(double ptmin, double ptmax, bool addUnits);

  const string sALICE      = "ALICE";
  const string sAntikt     = "Anti-#it{k}_{T}";
  const string sCharged    = "ch";
  const string sCounts     = "Counts";
  const string sEscheme    = "#it{E}-scheme";
  const string sEta        = "#it{#eta}";
  const string sGevC       = "GeV/#it{c}";
  const string sGevCC      = "GeV/#it{c}^{2}";
  const string sJet        = "jet";
  const string sJets       = "jets";
  const string sMass       = "#it{M}";
  const string sNumber     = "#it{N}";
  const string sPtscheme   = "#it{p}_{T}-scheme";
  const string sRadius     = "#it{R}";
  const string sRatio      = "Ratio";
  const string sSigma      = "#sigma";
  const string sSqrtS      = "#sqrt{#it{s}} = 13.6 TeV";
  const string sPpData     = "pp data";
  const string sPythia     = "PYTHIA";
  const string sThisThesis = "This Thesis";

  // Particles
  const string sV0         = "V0";
  const string sK0S        = "K^{0}_{S}";
  const string sLambda     = "#Lambda";
  const string sAntiLambda = "#bar{#Lambda}";
  const string sPiPm       = "#pi^{#pm}";
  const string sPiPlus     = "#pi^{+}";
  const string sPiMinus    = "#pi^{-}";
  const string sPiZero     = "#pi^{0}";
  const string sProton     = "p";
  const string sAntiProton = "#bar{p}";

  // Strings derived from the ones above
  const string sAlicePpData = sALICE + " " + sPpData;
  const string sPythiaSim   = sPythia + " simulation pp";
  const string sThisThesisAliceData = sThisThesis + ", " + sAlicePpData;
  const string sThisThesisAliceSim  = sThisThesis + ", " + sALICE + " ";
  const string sThisThesisPythiaSim = sThisThesis + ", " + sPythiaSim;
  const string sAntiktJets = sAntikt + " " + sJets;

  const string sChJets      = sCharged + " " + sJets;
  const string sChV0Jets    = sCharged + "+" + sV0 + " " + sJets;
  const string sAliceData   = sALICE + " " + sPpData;
  const string sAlicePythia = sALICE + " " + sPythia;

  const string sNjets       = sNumber + "_{" + sJets + "}";
  const string sNevts       = sNumber + "_{evts}";
  const string sNV0         = sNumber + "_{" + sV0 + "}";
  const string sNK0S        = sNumber + "_{" + sK0S + "}";

  const string sEtaJet      = sEta + "_{" + sJet + "}";
  const string sEtaV0       = sEta + "_{" + sV0 + "}";
  const string sEtaK0S      = sEta + "_{" + sK0S + "}";

  const string sJetRadius = addSubscript(sRadius, sJet);
  const string sJetRadius04 = sJetRadius + " = 0.4";

  const string sEtaJetRange035 = "|" + sEtaJet + "| < 0.35";
  const string sEtaJetRange05  = "|" + sEtaJet + "| < 0.5";
  const string sEtaJetRange20  = "|" + sEtaJet + "| < 2.0";
  const string sEtaV0Range075  = "|" + sEtaV0 + "| < 0.75";
  const string sEtaV0Range09   = "|" + sEtaV0 + "| < 0.9";
  const string sEtaK0SRange09  = "|" + sEtaK0S + "| < 0.9";

  const string sPtJet        = getPtString(sJet);
  const string sPtV0         = getPtString(sV0);
  const string sZV0          = getZString(sV0);
  const string sPtK0S        = getPtString(sK0S);
  const string sZK0S         = getZString(sK0S);
  const string sPtLambda     = getPtString(sLambda);
  const string sZLambda      = getZString(sLambda);
  const string sPtAntiLambda = getPtString(sAntiLambda);
  const string sZAntiLambda  = getZString(sAntiLambda);

  const string sPtJetWithUnits = sPtJet + " (" + sGevC + ")";

  // Measurements
  const string sJetsPerEvent     = getOneOverString(sNevts) + " " + getdYdXString(sNjets, sPtJet);
  const string sV0PtPerEvt       = getOneOverString(sNevts) + " " + getdYdXString(sNV0, sPtV0);
  const string sK0SPtPerEvt      = getOneOverString(sNevts) + " " + getdYdXString(sNK0S, sPtK0S);
  const string sV0ZPerEvt        = getOneOverString(sNevts) + " " + getdYdXString(sNV0, sZV0);
  const string sK0SZPerEvt       = getOneOverString(sNevts) + " " + getdYdXString(sNK0S, sZK0S);

  const string sV0PtPerJet        = getOneOverString(sNjets) + " " + getdYdXString(sNV0, sPtV0);
  const string sK0SPtPerJet       = getOneOverString(sNjets) + " " + getdYdXString(sNK0S, sPtK0S);
  const string sV0ZPerJet         = getOneOverString(sNjets) + " " + getdYdXString(sNV0, sZV0);
  const string sK0SZPerJet        = getOneOverString(sNjets) + " " + getdYdXString(sNK0S, sZK0S);
  const string sLambdaZPerJet     = getOneOverString(sNjets) + " " + getdYdXString(sLambda, sZLambda);
  const string sAntiLambdaZPerJet = getOneOverString(sNjets) + " " + getdYdXString(sAntiLambda, sZAntiLambda);

  // Simulations
  const string sJetXsec = sSigma + "_{" + sJets + "}";
  const string sV0Xsec = sSigma + "_{" + sV0 + "}";
  const string sK0SXsec = sSigma + "_{" + sK0S + "}";
  const string sJetsPerXsec = getOneOverString(sSigma) + " " + getdYdXString(sJetXsec, sPtJet);
  const string sV0PtPerXsec = getOneOverString(sSigma) + " " + getdYdXString(sV0Xsec, sPtV0);
  const string sK0SPtPerXsec = getOneOverString(sSigma) + " " + getdYdXString(sK0SXsec, sPtK0S);
  const string sK0SZPerJetXsec = getOneOverString(sJetXsec) + " " + getdYdXString(sK0SXsec, sZK0S);
}

string mystrings::addSubscript(string base, string subscript) {
  return base + "_{" + subscript + "}";
}

string mystrings::addSuperscript(string base, string superscript) {
  return base + "^{" + superscript + "}";
}

// Formats the hadron name to look nice (Greek letters, sub- and superscripts)
string mystrings::formatHadronName(string hadron) {
  string had = hadron;
  if (hadron == "pi") {
    had = sPiPm;
  } else if (hadron == "piplus") {
    had = sPiPlus;
  } else if (hadron == "piminus") {
    had = sPiMinus;
  } else if (hadron == "pi0") {
    had = sPiZero;
  } else if (hadron == "p") {
    had = sProton;
  } else if (hadron == "pbar") {
    had = sAntiProton;
  } else if (hadron == "K0S") {
    had = sK0S;
  } else if (hadron == "Lambda") {
    had = sLambda;
  } else if (hadron == "AntiLambda") {
    had = sAntiLambda;
  }
  return had;
}

// Returns decay products given a hadron
string mystrings::formatHadronDaughters(string hadron) {
  string daughters = hadron;
  if ("K0S" == hadron) {
    daughters = sPiPlus + sPiMinus;
  } else if ("Lambda" == hadron) {
    daughters = sProton + sPiMinus;
  } else if ("AntiLambda" == hadron) {
    daughters = sAntiProton + sPiPlus;
  }
  return daughters;
}

string mystrings::invMassOf(string system) {
  return "#it{M}(" + system + ")";
}

string mystrings::massOf(string particle) {
  return "#it{m}(" + particle + ")";
}

string mystrings::getPtString(string subscript) {
  string p = "#it{p}";
  string T = "T";
  if (subscript.empty())
    return addSubscript(p, T);
  else
    return addSubscript(p, (T + ", " + subscript).c_str());
}

string mystrings::getZString(string subscript) {
  if (subscript.empty())
    return TString::Format("#it{z}").Data();
  else
    return TString::Format("#it{z}_{%s}", subscript.c_str()).Data();
}

string mystrings::getRatioString(string num, string den) {
  return TString::Format("#frac{%s}{%s}", num.c_str(), den.c_str()).Data();
}

string mystrings::getOneOverString(string s) {
  return (getRatioString("1", s));
}

string mystrings::getdYdXString(string y, string x) {
  return getRatioString("d" + y, "d" + x);
}

string mystrings::getdYdPtString(string y) {
  return getdYdXString(y, getPtString(y));
}

string mystrings::getdYdZString(string y) {
  return getdYdXString(y, getZString(y));
}

string mystrings::getVarRangeString(double low, string var, double high) {
  string sLow  = TString::Format("%.0f", low).Data();
  string sHigh = TString::Format("%.0f", high).Data();
  return TString::Format("%s < %s < %s", sLow.c_str(), var.c_str(), sHigh.c_str()).Data();
}

string mystrings::getVarRangeString(string var, double high) {
  string sHigh = TString::Format("%.0f", high).Data();
  return TString::Format("%s < %s", var.c_str(), sHigh.c_str()).Data();
}

string mystrings::getPtJetRangeString(double ptmin, double ptmax, bool addUnits = true) {
  string s = TString::Format("%.f < %s < %.f", ptmin, sPtJet.c_str(), ptmax).Data();
  if (addUnits)
    s += TString::Format(" %s", sGevC.c_str()).Data();

  return s;
}

string mystrings::getPtV0RangeString(double ptmin, double ptmax, bool addUnits = true) {
  string s = TString::Format("%.1f < %s < %.1f", ptmin, sPtV0.c_str(), ptmax).Data();
  if (addUnits)
    s += TString::Format(" %s", sGevC.c_str()).Data();

  return s;
}

#endif

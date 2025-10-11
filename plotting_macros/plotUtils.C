
#include "TColor.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "histUtils.C"

#ifndef PLOTUTILS_H
#define PLOTUTILS_H

namespace plotutils {

static  int      myDarkRed     = TColor::GetColor(128,0,0);
static  int      myDarkGreen   = TColor::GetColor(0,128,0);
static  int      myDarkBlue    = TColor::GetColor(0,0,128);

TH1F *DrawFrame(Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, TString xTitle = "", TString yTitle = "", bool setMargins = true) {
  if(setMargins) {
    gPad->SetLeftMargin(0.22);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    // gPad->SetRightMargin(0.15);
    // gPad->SetTopMargin(0.1);
  }

  TH1F *frame = gPad->DrawFrame(xmin,ymin,xmax,ymax);
  frame->SetXTitle(xTitle.Data());
  frame->SetYTitle(yTitle.Data());
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetXaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitleOffset(1.6);
  frame->GetXaxis()->CenterTitle(true);
  frame->GetYaxis()->CenterTitle(true);

  gPad->SetTicks(1,1);
  return frame;
}

TH1F* DrawFrame(TH1* hist, bool setMargins = true) {
  double xMinFrame = hist->GetXaxis()->GetXmin();
  double xMaxFrame = hist->GetXaxis()->GetXmax();
  double yMinFrame = histutils::getLowerBound(hist, 0) * 0.9;
  double yMaxFrame = histutils::getUpperBound(hist, 0) * 1.2;
  string xTitle = hist->GetXaxis()->GetTitle();
  string yTitle = hist->GetYaxis()->GetTitle();
  return DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle, setMargins);
}

TLegend *CreateLegend(Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, TString title = "", Double_t textSize = 0.06) {
  TLegend *leg = new TLegend(xmin,ymin,xmax,ymax,title.Data());
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(textSize);
  leg->SetTextFont(42);
  return leg;
}

TLatex* CreateLatex(Double_t x, Double_t y, TString strText = "", Double_t textSize = 0.06, Int_t color = 1, int all = 11, bool bndc = true) {
  TLatex* text = new TLatex(x, y, strText.Data());
  if(bndc) text->SetNDC();
  text->SetTextAlign(all);
  text->SetTextSize(textSize);
  text->SetTextFont(42);
  text->SetTextColor(color);
  return text;
}

void DrawLatex(Double_t x, Double_t y, TString strText = "", Double_t textSize = 0.06, Int_t color = 1, int all = 11, bool bndc = true) {
  TLatex text;
  if(bndc) text.SetNDC();
  text.SetTextAlign(all);
  text.SetTextSize(textSize);
  text.SetTextFont(42);
  text.SetTextColor(color);
  text.DrawLatex(x,y,strText.Data());
}

Int_t GetColor(Int_t i) {
  const Int_t nc = 11;
  // Int_t color[nc] = {1,kRed+1,4,kGreen+2,kAzure-2,kOrange+7,kGray+2,myDarkRed,kAzure+10,kGreen+4,kYellow+2};
  Int_t color[nc] = {1,kRed+1,4,kGreen+2,kCyan+4,kOrange+7,kGray+2,myDarkRed,kAzure+10,kGreen+4,kYellow+2};
  if(i<nc) return color[i];
  else     return i;
}

Int_t GetFillColor(Int_t i) {
  const Int_t nc = 11;
  Int_t color[nc] = {kGray,kRed-9,kAzure+1,kGreen-6,kOrange-9,kGray,kRed-6,kAzure+6,kGreen-5,kYellow-3,kGreen-9};
  if(i<nc) return color[i];
  else     return i;
}

Int_t GetMarker(Int_t i) {
  const Int_t nc = 8;
  Int_t markerStyle[nc] = {20,21,33,34,24,25,27,28};
  if(i<nc) return markerStyle[i];
  else     return 20+i;
}

// Set histogram colours and markers
void setStyle(TH1* hist, int styleNumber, double alpha = -1, int lineWidth = 3) {
  hist->SetLineWidth(lineWidth);
  hist->SetLineColor(GetColor(styleNumber));
  hist->SetMarkerStyle(GetMarker(styleNumber));
  hist->SetMarkerColor(GetColor(styleNumber));
  if (alpha > 0.) {
    hist->SetFillColorAlpha(GetColor(styleNumber), alpha);
  }
}

void setStyle(TLine* line, int styleNumber, int lineStyle = 9, int lineWidth = 3) {
  line->SetLineWidth(lineWidth);
  line->SetLineColor(GetColor(styleNumber));
  line->SetLineStyle(lineStyle);
}

void setStyle(TF1* f, int styleNumber, int lineStyle = 1, int lineWidth = 3) {
  f->SetLineWidth(lineWidth);
  f->SetLineColor(GetColor(styleNumber));
  f->SetLineStyle(lineStyle);
}

void SetPadMargins(double left, double right, double bottom, double top) {
  gPad->SetLeftMargin(left);
  gPad->SetRightMargin(right);
  gPad->SetBottomMargin(bottom);
  gPad->SetTopMargin(top);
}

struct Plotter {
  private:
    string _drawOption; // Private because it requires caution with spaces
    bool   _logPlot;
    string _outputFileName;
    double _textSize;

    TH1F*            _frame   = nullptr;
    TCanvas*         _canvas  = nullptr;
    TLegend*         _legend  = nullptr;
    vector<TH1*>     _hists   = {};
    vector<TObject*> _objects = {};

    bool isHistVectorEmpty(string origin) {
      if (_hists.empty()) {
        std::cout << "Plotter::" << origin << "(): Hist vector is empty!" << std::endl;
        return true;
      }
      return false;
    }
  public:
    Plotter(string ofn = "", bool lp = false, double ts = 0.04) : _outputFileName(ofn), _logPlot(lp), _textSize(ts) { Plotter::reset(); }

    string getDrawOption() { return _drawOption; }
    bool   getLogPlot() { return _logPlot; }
    string getOutputFileName() { return _outputFileName; }
    double getTextSize() { return _textSize; }

    void setDrawOption(string s) {
      // If the user already added a space at the start, don't add another one
      if (s.at(0) == ' ') { // Single quotes, because checking for a char
        _drawOption = s;
      } else {
        _drawOption = " " + s;
      }
    }
    void setLogPlot(bool x = true) { _logPlot = x; }
    void setOutputFileName(string s) { _outputFileName = s; }
    void setTextSize(double x) { _textSize = x; }

    // Utilities
    void addLatex(double x, double y, string s) {
      _objects.push_back(CreateLatex(x, y, s.c_str(), _textSize));
    }
    template <typename T>
    void addLegendEntry(T* h, string s, string option = "lp") {
      if (_legend)
        _legend->AddEntry(h, s.c_str(), option.c_str());
      else
        std::cout << "Plotter::addLegendEntry(): Legend does not exist!" << std::endl;
    }
    void addHistogram(TH1* h) {
      _hists.push_back(h);
    }
    void addLine(double x0, double y0, double x1, double y1, int styleNumber, int lineStyle = 9, int lineWidth = 3) {
      TLine* l = new TLine(x0, y0, x1, y1);
      setStyle(l, styleNumber, lineStyle, lineWidth);
      _objects.push_back(l);
    }
    void makeCanvas(string s = "c", double x = 800, double y = 600) {
      _canvas = new TCanvas(s.c_str(), s.c_str(), x, y);
      _canvas->SetLogy(_logPlot);
    }
    void makeFrame(double x0, double x1, double y0, double y1, string sx, string sy) {
      if (!_canvas) makeCanvas();
      _frame = DrawFrame(x0, x1, y0, y1, sx, sy);
    }
    void makeFrame(string sx, string sy) {
      // Automatically determine frame ranges from hists vector
      if (isHistVectorEmpty("makeFrame"))
        return;

      if (!_canvas) makeCanvas();
      double xMinFrame = _hists[0]->GetXaxis()->GetXmin();
      double xMaxFrame = _hists[0]->GetXaxis()->GetXmax();
      double yMinFrame = histutils::getLowerBound(_hists, 0, 0) * 0.9;
      double yMaxFrame = histutils::getUpperBound(_hists, 0, 0) * 1.2;
      makeFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, sx, sy);
    }
    void makeLegend(double x0, double x1, double y0, double y1, string s) {
      _legend = CreateLegend(x0, x1, y0, y1, s.c_str(), _textSize);
    }
    void makeRatios(TH1* baseHist) {
      TH1* baseCopy = (TH1*)baseHist->Clone("baseCopy");
      for (auto& h : _hists) h = histutils::divideWithProtection(h, baseCopy);
    }
    void makeRatios(int baseIndex = 0) {
      if (isHistVectorEmpty("makeRatios"))
        return;

      TH1* baseCopy = (TH1*)_hists[baseIndex]->Clone("baseCopy");
      makeRatios(baseCopy);
    }
    void plot() {
      if (!_canvas) makeCanvas();
      if (!_frame) makeFrame("x", "y");

      _frame->Draw();
      if (_legend) _legend->Draw("same");
      for (auto o : _objects)
        o->Draw("same");

      for (auto h : _hists)
        h->Draw(("same" + _drawOption).c_str());

      if (_outputFileName != "")
        _canvas->SaveAs(_outputFileName.c_str());
    }
    void resetCanvas()  { _canvas = nullptr; }
    void resetFrame()   { _frame = nullptr; }
    void resetLegend()  { _legend = nullptr; }
    void resetObjects() { _objects.clear(); }
    void resetHists()   { _hists.clear(); }
    void reset() {
      resetCanvas();
      resetFrame();
      resetLegend();
      resetObjects();
      resetHists();
    }
    void setHists(vector<TH1*> h) {
      _hists = h;
    }
    void setHistStyles() {
      for (unsigned int i = 0; i < _hists.size(); i++)
        setStyle(_hists[i], i);
    }
}; // struct Plotter

} // namespace plotutils

#endif // PLOTUTILS_H

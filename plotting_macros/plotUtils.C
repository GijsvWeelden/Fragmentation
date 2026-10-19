
#include "TColor.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "histUtils.C"

#ifndef PLOTUTILS_H
#define PLOTUTILS_H

namespace plotutils {

namespace defaultValues {
  const double defPadMarginLeft   = 0.22;
  const double defPadMarginRight  = 0.05;
  const double defPadMarginBottom = 0.15;
  const double defPadMarginTop    = 0.05;

  const double defTextSize         = 0.04;
  const double defLabelSize        = 0.05;
  const double defTitleSize        = 0.06;
  const double defTitleOffset      = 1.0;
  const double defTitleOffsetLarge = 1.6;

  const int defTicks = 1;
} // namespace defaultValues

static  int      myDarkRed     = TColor::GetColor(128,0,0);
static  int      myDarkGreen   = TColor::GetColor(0,128,0);
static  int      myDarkBlue    = TColor::GetColor(0,0,128);

void SetPadMargins(double left, double right, double bottom, double top) {
  gPad->SetLeftMargin(left);
  gPad->SetRightMargin(right);
  gPad->SetBottomMargin(bottom);
  gPad->SetTopMargin(top);
}

void SetFrameLabelSizes(TH1F* frame, double x, double y) {
  frame->GetXaxis()->SetLabelSize(x);
  frame->GetYaxis()->SetLabelSize(y);
}

void SetFrameTitleSizes(TH1F* frame, double x, double y) {
  frame->GetXaxis()->SetTitleSize(x);
  frame->GetYaxis()->SetTitleSize(y);
}

void SetFrameTitleOffsets(TH1F* frame, double x, double y) {
  frame->GetXaxis()->SetTitleOffset(x);
  frame->GetYaxis()->SetTitleOffset(y);
}

void SetFrameTitleCenter(TH1F* frame, bool b) {
  frame->GetXaxis()->CenterTitle(b);
  frame->GetYaxis()->CenterTitle(b);
}

TH1F *DrawFrame(Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, TString xTitle = "", TString yTitle = "", bool setMargins = true) {
  if (setMargins)
    SetPadMargins(defaultValues::defPadMarginLeft, defaultValues::defPadMarginRight, defaultValues::defPadMarginBottom, defaultValues::defPadMarginTop);
  
  TH1F *frame = gPad->DrawFrame(xmin,ymin,xmax,ymax);
  frame->SetXTitle(xTitle.Data());
  frame->SetYTitle(yTitle.Data());
  SetFrameLabelSizes(frame, defaultValues::defLabelSize, defaultValues::defLabelSize);
  SetFrameTitleSizes(frame, defaultValues::defTitleSize, defaultValues::defTitleSize);
  SetFrameTitleOffsets(frame, defaultValues::defTitleOffset, defaultValues::defTitleOffsetLarge);
  SetFrameTitleCenter(frame, true);

  gPad->SetTicks(defaultValues::defTicks, defaultValues::defTicks);
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

TLegend *CreateLegend(Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, TString title = "", Double_t textSize = defaultValues::defTextSize) {
  TLegend *leg = new TLegend(xmin,ymin,xmax,ymax,title.Data());
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(textSize);
  leg->SetTextFont(42);
  return leg;
}

TLatex* CreateLatex(Double_t x, Double_t y, TString strText = "", Double_t textSize = defaultValues::defTextSize, Int_t color = 1, int all = 11, bool bndc = true) {
  TLatex* text = new TLatex(x, y, strText.Data());
  if(bndc) text->SetNDC();
  text->SetTextAlign(all);
  text->SetTextSize(textSize);
  text->SetTextFont(42);
  text->SetTextColor(color);
  return text;
}

void DrawLatex(Double_t x, Double_t y, TString strText = "", Double_t textSize = defaultValues::defTextSize, Int_t color = 1, int all = 11, bool bndc = true) {
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
  // Int_t color[nc] = {1,kRed+1,4,kGreen+2,kCyan+4,kOrange+7,kGray+2,myDarkRed,kAzure+10,kGreen+4,kYellow+2};
  Int_t color[nc] = {1,kRed+1,4,kGreen+2,kMagenta+2,kOrange+7,kGray+2,myDarkRed,kAzure+10,kGreen+4,kYellow+2};
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

struct Plotter {
  private:
    string _drawOption = ""; // Private because it requires caution with spaces
    double _frameLabelSize    = defaultValues::defLabelSize;
    int    _frameTicks        = defaultValues::defTicks;
    bool   _frameTitleCenter  = true;
    double _frameTitleSize    = defaultValues::defTitleSize;
    double _frameTitleOffsetX = defaultValues::defTitleOffset;
    double _frameTitleOffsetY = defaultValues::defTitleOffsetLarge;
    bool   _logx = false, _logy = false, _logz = false;
    string _outputFileName = "";
    double _padMarginBottom = defaultValues::defPadMarginBottom;
    double _padMarginLeft   = defaultValues::defPadMarginLeft;
    double _padMarginRight  = defaultValues::defPadMarginRight;
    double _padMarginTop    = defaultValues::defPadMarginTop;
    double _textSize = defaultValues::defTextSize;

    TH1F*            _frame   = nullptr;
    TCanvas*         _canvas  = nullptr;
    TLegend*         _legend  = nullptr;
    vector<TH1*>     _hists   = {};
    vector<TObject*> _objects = {};

    bool isHistVectorEmpty(string origin) {
      if (_hists.empty()) {
        std::cout << "plotutils::Plotter::" << origin << "(): Hist vector is empty!" << std::endl;
        return true;
      }
      return false;
    }
  public:
    Plotter(string ofn, bool ly) : _outputFileName(ofn), _logy(ly) { Plotter::reset(); }
    Plotter(string ofn, bool lx, bool ly, bool lz) : _outputFileName(ofn), _logx(lx), _logy(ly), _logz(lz) { Plotter::reset(); }

    TCanvas* getCanvas() { return _canvas; }
    string getDrawOption() { return _drawOption; }
    TH1F* getFrame() { return _frame; }
    double getFrameLabelSize() { return _frameLabelSize; }
    int    getFrameTicks() { return _frameTicks; }
    bool   getFrameTitleCenter() { return _frameTitleCenter; }
    double getFrameTitleSize() { return _frameTitleSize; }
    double getFrameTitleOffsetX() { return _frameTitleOffsetX; }
    double getFrameTitleOffsetY() { return _frameTitleOffsetY; }
    vector<TH1*> getHists() { return _hists; }
    TLegend* getLegend() { return _legend; }
    bool getLogX() { return _logx; }
    bool getLogY() { return _logy; }
    bool getLogZ() { return _logz; }
    vector<TObject*> getObjects() { return _objects; }
    string getOutputFileName() { return _outputFileName; }
    double getPadMarginBottom() { return _padMarginBottom; }
    double getPadMarginLeft() { return _padMarginLeft; }
    double getPadMarginRight() { return _padMarginRight; }
    double getPadMarginTop() { return _padMarginTop; }
    double getTextSize() { return _textSize; }

    void setDrawOption(string s) {
      // If the user already added a space at the start, don't add another one
      if (s.at(0) == ' ') { // Single quotes, because checking for a char
        _drawOption = s;
      } else {
        _drawOption = " " + s;
      }
    }
    void setFrameLabelSize(double x) { _frameLabelSize = x; }
    void setFrameTitleCenter(bool b) { _frameTitleCenter = b; }
    void setFrameTicks(int x) { _frameTicks = x; }
    void setFrameTitleSize(double x) { _frameTitleSize = x; }
    void setFrameTitleOffsetX(double x) { _frameTitleOffsetX = x; }
    void setFrameTitleOffsetY(double y) { _frameTitleOffsetY = y; }
    void setLogX(bool b = true) { _logx = b; }
    void setLogY(bool b = true) { _logy = b; }
    void setLogZ(bool b = true) { _logz = b; }
    void setOutputFileName(string s) { _outputFileName = s; }
    void setPadMargins(double left, double right, double bottom, double top) {
      _padMarginLeft = left;
      _padMarginRight = right;
      _padMarginBottom = bottom;
      _padMarginTop = top;
    }
    void setTextSize(double x) { _textSize = x; }
    void setZAxisRange(double min, double max) {
      if (isHistVectorEmpty("setZAxisRange"))
        return;

      for (auto h : _hists)
        histutils::setHistMinMax(h, min, max);
    }    

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
    void addHists(vector<TH1*> v) {
      for (auto h : v) {
        _hists.push_back(h);
      }
    }
    void addLine(double x0, double y0, double x1, double y1, int styleNumber, int lineStyle = 9, int lineWidth = 3) {
      TLine* l = new TLine(x0, y0, x1, y1);
      setStyle(l, styleNumber, lineStyle, lineWidth);
      _objects.push_back(l);
    }
    void addObject(TObject* obj) {
      _objects.push_back(obj);
    }
    void makeCanvas(string s = "c", double x = 800, double y = 600) {
      _canvas = new TCanvas(s.c_str(), s.c_str(), x, y);
      _canvas->SetLogx(_logx);
      _canvas->SetLogy(_logy);
      _canvas->SetLogz(_logz);
    }
    void makeFrame(double x0, double x1, double y0, double y1, string sx, string sy) {
      if (!_canvas) makeCanvas();
      _frame = DrawFrame(x0, x1, y0, y1, sx, sy);
      
      SetFrameLabelSizes(_frame, _frameLabelSize, _frameLabelSize);
      SetFrameTitleSizes(_frame, _frameTitleSize, _frameTitleSize);
      SetFrameTitleOffsets(_frame, _frameTitleOffsetX, _frameTitleOffsetY);
      SetFrameTitleCenter(_frame, _frameTitleCenter);
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
      if (!_canvas)
        makeCanvas();
      if (!_frame)
        makeFrame("x", "y");

      SetPadMargins(_padMarginLeft, _padMarginRight, _padMarginBottom, _padMarginTop);
      gPad->SetTicks(_frameTicks, _frameTicks);

      _frame->Draw();
      if (_legend)
        _legend->Draw("same");
      else
        cout << "Warning: plotutils::Plotter::plot() no legend\n";

      for (auto o : _objects)
        o->Draw("same");

      for (auto h : _hists)
        h->Draw(("same" + _drawOption).c_str());

      if (_outputFileName != "")
        _canvas->SaveAs(_outputFileName.c_str());
      else
        cout << "Warning: plotutils::Plotter::plot() No outputfilename selected, plot not saved\n";
    }
    void resetCanvas()  { _canvas = nullptr; }
    void resetFrame()   { _frame = nullptr; }
    void resetLegend()  { _legend = nullptr; }
    void resetObjects() { _objects.clear(); }
    void resetHists()   { _hists.clear(); }
    void reset() {
      gStyle->SetNdivisions(505, "xy");
      resetCanvas();
      resetFrame();
      resetLegend();
      resetObjects();
      resetHists();
    }
    void selfNormaliseHists() {
      if (isHistVectorEmpty("selfNormaliseHists"))
        return;

      for (auto& h : _hists) {
        double integral = h->Integral("width");
        if (std::isnan(integral))
          cout << "Integral is nan! Skipping hist " << h->GetName() << endl;
        else if (std::abs(integral) < 1e-25)
          cout << "Integral is < 1e-25! Skipping hist " << h->GetName() << endl;
        else
          h->Scale(1. / integral);
      }
    }
    void setHists(vector<TH1*> h) {
      _hists = h;
    }
    void setHistStyles() {
      for (unsigned int i = 0; i < _hists.size(); i++)
        setStyle(_hists[i], i);
    }
}; // struct Plotter

void ColorTester() {
  double x0 = 0., x1 = 1., y0 = 0., y1 = 1.;
  TCanvas* c = new TCanvas("c", "c", 800, 600);

  TLine* l1 = new TLine(x0, y0, x1, y1);
  TLine* l2 = new TLine(x0, 0.2*y1, x1, 0.8*y1);
  TLine* l3 = new TLine(x0, 0.4*y1, x1, 0.6*y1);
  TLine* l4 = new TLine(x0, 0.6*y1, x1, 0.4*y1);
  TLine* l5 = new TLine(x0, 0.8*y1, x1, 0.2*y1);
  TLine* l6 = new TLine(x0, y1, x1, y0);

  l1->SetLineColor(kBlack);
  l2->SetLineColor(kRed+1);
  l3->SetLineColor(4);
  l4->SetLineColor(kGreen+2);
  l5->SetLineColor(kMagenta);
  l6->SetLineColor(kMagenta+2);

  l1->SetLineWidth(3);
  l2->SetLineWidth(3);
  l3->SetLineWidth(3);
  l4->SetLineWidth(3);
  l5->SetLineWidth(3);
  l6->SetLineWidth(3);

  TH1* h = new TH1F("h", "", 1, x0, x1);
  h->SetStats(0);
  vector<TLine*> lines = {l1, l2, l3, l4, l5, l6};

  h->Draw();
  for (auto& line : lines)
    line->Draw("same");

  gPad->SetLeftMargin(0);
  gPad->SetRightMargin(0);
  gPad->SetBottomMargin(0);
  gPad->SetTopMargin(0);

  c->SaveAs("colortest.pdf");
}

} // namespace plotutils

#endif // PLOTUTILS_H

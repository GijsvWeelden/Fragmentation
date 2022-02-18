/*
 *   Some utilities for plotting legends in ROOT
 *   By Marco van Leeuwen
 *
 */


Int_t __legend_font = 62;
Float_t __legend_size = 0.05;

void set_legend_font(Int_t font, Float_t size) {
  __legend_font = font;
  if (font%10 == 3) {
    if (size < 5) {
      cout << "WARNING: font pixel size too small, using default=12 instead" << endl;
      size = 16;
    }
  }
  else {
    if (size > 0.5) {
      cout << "WARNING: font scalable too large, using default=0.05 instead" << endl;
      size = 0.05;
    }
  }
  __legend_size = size;
}

Float_t ndc_x(Float_t x) {
  return (1-TVirtualPad::Pad()->GetLeftMargin()-TVirtualPad::Pad()->GetRightMargin())*x+TVirtualPad::Pad()->GetLeftMargin();
}

Float_t ndc_y(Float_t y) {
  return (1-TVirtualPad::Pad()->GetTopMargin()-TVirtualPad::Pad()->GetBottomMargin())*y+TVirtualPad::Pad()->GetBottomMargin();
}

void draw_fill_graph(TGraph *gr, const Char_t *opt="p") {
  TGraph *gr_fill=(TGraph*) gr->Clone();
  gr_fill->SetMarkerStyle(gr->GetMarkerStyle()-4);
  gr_fill->SetMarkerColor(10);
  Char_t optx[25];
  sprintf(optx,"%sx",opt);
  gr_fill->Draw(optx);
  gr->Draw(opt);
}

void draw_legend_b(Float_t x, Float_t y, TAttLine *att_l, TAttFill *att_f, const Char_t *label, Int_t i_ent, Int_t colored=1, Size_t txt_size=-1) {
  //y-=0.07*i_ent;
  if (txt_size < 0) 
    txt_size = __legend_size; // Get default (global var)

  Float_t pix_y_siz=TVirtualPad::Pad()->AbsPixeltoY(1)-TVirtualPad::Pad()->AbsPixeltoY(2);
  Float_t y_min_real=TVirtualPad::Pad()->GetUymin()+y*(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin())-0.4*pix_y_siz*txt_size;
  Float_t y_max_real=TVirtualPad::Pad()->GetUymin()+y*(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin())+0.4*pix_y_siz*txt_size;

  if (__legend_font%10 == 3) {
    y-=pix_y_siz*txt_size*1.1*i_ent/(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin());
  }
  else {
    y-=txt_size*1.2*i_ent;
  }
  Float_t y_real=TVirtualPad::Pad()->GetUymin()+y*(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin());

  Float_t x_min_real=TVirtualPad::Pad()->GetUxmin()+(x+0.01)*(TVirtualPad::Pad()->GetUxmax()-TVirtualPad::Pad()->GetUxmin());
  Float_t x_max_real=TVirtualPad::Pad()->GetUxmin()+(x+0.08)*(TVirtualPad::Pad()->GetUxmax()-TVirtualPad::Pad()->GetUxmin());
  if (TVirtualPad::Pad()->GetLogx()) {
    x_min_real=exp(log(10)*x_min_real);
    x_max_real=exp(log(10)*x_max_real);
  }
  if (TVirtualPad::Pad()->GetLogy()) {
    y_real=exp(log(10)*y_real);
  }
  if (__legend_font%10 != 3)
    /*  {
    y_min_real -= pix_y_siz*txt_size*1.1*i_ent/(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin());
    y_max_real += pix_y_siz*txt_size*1.1*i_ent/(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin());
  }
  else*/
    {
      y_min_real = y_real-0.5*txt_size*(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin());
      y_max_real = y_real+0.5*txt_size*(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin());
  }
  
  if (att_f) {
    TBox *box_f=new TBox(x_min_real,y_min_real,x_max_real,y_max_real);
    box_f->SetFillStyle(att_f->GetFillStyle());
    if (colored!=0)
      box_f->SetFillColor(att_f->GetFillColor());
    box_f->Draw();
  }
  if (att_l) {
    TBox *box_l=new TBox(x_min_real,y_min_real,x_max_real,y_max_real);
    box_l->SetLineStyle(att_l->GetLineStyle());
    box_l->SetLineWidth(att_l->GetLineWidth());
    if (colored!=0)
      box_l->SetLineColor(att_l->GetLineColor());
    box_l->SetFillStyle(0);
    box_l->Draw();
  }
  TLatex *ltx=new TLatex();
  ltx->SetTextFont(__legend_font);
  ltx->SetTextSize(txt_size);
  ltx->SetTextAlign(12);
  if (colored==2)
    ltx->SetTextColor(att_l->GetLineColor());
  ltx->SetNDC();
  ltx->DrawLatex(ndc_x(x+0.1),ndc_y(y),label);
}

void draw_legend_b2(Float_t x, Float_t y, TAttLine *att_l, TAttFill *att_f, const Char_t *label, Int_t i_ent, Int_t colored=1, Size_t txt_size=16) {
  //y-=0.07*i_ent;
  Float_t pix_y_siz=TVirtualPad::Pad()->AbsPixeltoY(1)-TVirtualPad::Pad()->AbsPixeltoY(2);
  y-=pix_y_siz*txt_size*1.1*i_ent/(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin());

  Float_t x_min_real=TVirtualPad::Pad()->GetUxmin()+(x+0.01)*(TVirtualPad::Pad()->GetUxmax()-TVirtualPad::Pad()->GetUxmin());
  Float_t x_max_real=TVirtualPad::Pad()->GetUxmin()+(x+0.08)*(TVirtualPad::Pad()->GetUxmax()-TVirtualPad::Pad()->GetUxmin());
  Float_t y_min_real=TVirtualPad::Pad()->GetUymin()+y*(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin())-0.4*pix_y_siz*txt_size;
  Float_t y_max_real=TVirtualPad::Pad()->GetUymin()+y*(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin())+0.4*pix_y_siz*txt_size;
  if (att_f) {
    TBox *box_f=new TBox(x_min_real,y_min_real,x_max_real,y_max_real);
    box_f->SetFillStyle(att_f->GetFillStyle());
    if (colored!=0)
      box_f->SetFillColor(att_f->GetFillColor());
    box_f->Draw();
  }
  if (att_l) {
    TLine *box_l=new TLine(x_min_real,y_min_real,x_max_real,y_min_real);
    box_l->SetLineStyle(att_l->GetLineStyle());
    box_l->SetLineWidth(att_l->GetLineWidth());
    if (colored!=0)
      box_l->SetLineColor(att_l->GetLineColor());
    box_l->Draw();
    box_l->DrawLine(x_min_real,y_max_real,x_max_real,y_max_real);
  }
  TLatex *ltx=new TLatex();
  ltx->SetTextFont(133);
  ltx->SetTextSize(txt_size);
  ltx->SetTextAlign(12);
  if (colored==2)
    ltx->SetTextColor(att_l->GetLineColor());
  ltx->SetNDC();
  ltx->DrawLatex(ndc_x(x+0.1),ndc_y(y),label);
}

void draw_legend_l(Float_t x, Float_t y, TAttLine *att_l, const Char_t *label, Int_t i_ent, Int_t colored = 1, Size_t txt_size=-1, Int_t txt_font=-1) {
  //y-=0.07*i_ent;
  Float_t pix_y_siz=TVirtualPad::Pad()->AbsPixeltoY(1)-TVirtualPad::Pad()->AbsPixeltoY(2);
  if (txt_font < 0) 
    txt_font = __legend_font; // Get default (global var)
  if (txt_size < 0) 
    txt_size = __legend_size; // Get default (global var)

  if (txt_font%10 == 3) {
    y-=pix_y_siz*txt_size*1.1*i_ent/(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin());
  }
  else {
    y-=txt_size*1.2*i_ent;
  }
  Float_t x_min_real=TVirtualPad::Pad()->GetUxmin()+(x+0.01)*(TVirtualPad::Pad()->GetUxmax()-TVirtualPad::Pad()->GetUxmin());
  Float_t x_max_real=TVirtualPad::Pad()->GetUxmin()+(x+0.08)*(TVirtualPad::Pad()->GetUxmax()-TVirtualPad::Pad()->GetUxmin());
  Float_t y_real=TVirtualPad::Pad()->GetUymin()+y*(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin());
  Int_t l_style=att_l->GetLineStyle();
  if (TVirtualPad::Pad()->GetLogx()) {
    x_min_real=exp(log(10)*x_min_real);
    x_max_real=exp(log(10)*x_max_real);
  }
  if (TVirtualPad::Pad()->GetLogy())
    y_real=exp(log(10)*y_real);
  TLine *line=new TLine(x_min_real,y_real,x_max_real,y_real);
  Float_t l_siz=att_l->GetLineWidth();
  line->SetLineWidth(l_siz);
  line->SetLineStyle(l_style);
  line->SetLineColor(att_l->GetLineColor());
  line->Draw();
  TLatex *ltx=new TLatex();
  ltx->SetTextFont(txt_font);
  ltx->SetTextSize(txt_size);
  ltx->SetTextAlign(12);
  if (colored==2)
    ltx->SetTextColor(att_l->GetLineColor());
  ltx->SetNDC();
  ltx->DrawLatex(ndc_x(x+0.1),ndc_y(y),label);
}

void draw_legend_m(Float_t x, Float_t y, TAttMarker *att_m, const Char_t *label, Float_t i_ent, Int_t colored=1, Size_t txt_size=-1, Int_t txt_font=-1) {
  if (txt_font < 0) 
    txt_font = __legend_font; // Get default (global var)
  if (txt_size < 0) 
    txt_size = __legend_size; // Get default (global var)

  Float_t pix_y_siz=TVirtualPad::Pad()->AbsPixeltoY(1)-TVirtualPad::Pad()->AbsPixeltoY(2);
  if (txt_font%10 == 3) {
    y-=pix_y_siz*txt_size*1.1*i_ent/(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin());
  }
  else {
    y-=txt_size*1.2*i_ent;
  }

  Float_t x_real=TVirtualPad::Pad()->GetUxmin()+(x+0.035)*(TVirtualPad::Pad()->GetUxmax()-TVirtualPad::Pad()->GetUxmin());
  Float_t y_real=TVirtualPad::Pad()->GetUymin()+y*(TVirtualPad::Pad()->GetUymax()-TVirtualPad::Pad()->GetUymin());
  Int_t m_style=att_m->GetMarkerStyle();
  if (TVirtualPad::Pad()->GetLogy())
    y_real=exp(log(10)*y_real);
  if (TVirtualPad::Pad()->GetLogx())
    x_real=exp(log(10)*x_real);
  // cout << "uymin " << TVirtualPad::Pad()->GetUymin() << " max " 
  //     << TVirtualPad::Pad()->GetUymax() << " y " << y << " y_real " << y_real << endl; 
  TMarker *mrk=new TMarker(x_real,y_real,m_style);
  Float_t m_siz=att_m->GetMarkerSize();
  mrk->SetMarkerSize(m_siz);
  if (colored!=0)
    mrk->SetMarkerColor(att_m->GetMarkerColor());
  mrk->Draw();
  TLatex *ltx=new TLatex();
  ltx->SetTextFont(txt_font);
  ltx->SetTextSize(txt_size);
  ltx->SetTextAlign(12);
  if (colored==2)
    ltx->SetTextColor(att_m->GetMarkerColor());
  ltx->SetNDC();
  ltx->DrawLatex(ndc_x(x+0.08),ndc_y(y),label);
}

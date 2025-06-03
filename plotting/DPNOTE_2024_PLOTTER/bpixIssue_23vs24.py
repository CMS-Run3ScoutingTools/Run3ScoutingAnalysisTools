#!/usr/bin/env python3

# ***************************************
# usage: 
#    python3 EXO-23-007_scoutVSoffline.py
# ***************************************

import ROOT, array, random, copy
from ROOT import TCanvas, TFile, TH1, TH1F, TF1, gSystem
from ROOT import *
import ROOT, array, CMSGraphics, CMS_lumi, random, copy
from ROOT import RooCmdArg, RooArgSet, kFALSE, RooLinkedList, kBlue, kRed, kBlack, kOpenStar, kWhite, kGray
from ROOT import gStyle, TStyle, TGraph, TGraphErrors, TMath, TMultiGraph, TLine, gPad, TGaxis, TLegend, TText, TLatex, TColor, TPaveText
from ROOT import TAttFill, TLegend, TRatioPlot, TPad, THStack, TFileCollection
from ROOT import kBlue, kRed, kBlack, kWhite, kAzure, kOrange, kPink, kGreen, kYellow, kCyan
from array import array
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import argparse
import sys

argparser = argparse.ArgumentParser(description='Parser used for non default arguments', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=True)
argparser.add_argument('--outdir', dest='outdir', default='/eos/user/e/elfontan/www/CMS_SCOUTING/2024/DIMUON/DP_NOTE/', help='Output directory')
args = argparser.parse_args()
outputdir = args.outdir

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


######################################
# List of files and output directory #
######################################
def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]
#files_2024 = list_full_paths("/eos/user/e/elfontan/2024_SCOUTING/ScoutingPFRun3/")
#files_2024 = list_full_paths("/eos/user/e/elfontan/2024_SCOUTING/ScoutingPFMonitor/")
#files_2023C = list_full_paths("/eos/user/e/elfontan/2023_SCOUTING/ScoutingPFMonitor/")
#files_2023D = list_full_paths("/eos/user/e/elfontan/2023_SCOUTING/ScoutingPFMonitor/")
#files_2024 = files_2024[:-2] 
files_2024 = ROOT.TFile.Open("/eos/user/e/elfontan/2024_SCOUTING/noVtxMu_ScoutingPFMonitor/scoutMonitor_2024E_noVtxMu.root")
#files_2024 = ROOT.TFile.Open("/eos/user/e/elfontan/2024_SCOUTING/ScoutingPFMonitor/scoutMonitor_2024Ev1.root")
files_2023C = ROOT.TFile.Open("/eos/user/e/elfontan/2023_SCOUTING/ScoutingPFMonitor/scoutMonitor_2023Cv1.root")
files_2023D = ROOT.TFile.Open("/eos/user/e/elfontan/2023_SCOUTING/ScoutingPFMonitor/scoutMonitor_2023Dv1.root")


########################
# Variables and histos #
########################
h_2024_pt    = TH1F("h_2024_pt", "h_2024_pt", 100, 0, 50)
h_2024_eta    = TH1F("h_2024_eta", "h_2024_eta", 80, -4, 4)
h_2024_phi    = TH1F("h_2024_phi", "h_2024_phi", 80, -4, 4)
h_2024_etaphi = TH2F("h_2024_etaphi", "h_2024_etaphi", 80, -4, 4, 80, -4, 4)

h_2023C_pt    = TH1F("h_2023C_pt", "h_2023C_pt", 100, 0, 50)
h_2023C_eta    = TH1F("h_2023C_eta", "h_2023C_eta", 80, -4, 4)
h_2023C_phi    = TH1F("h_2023C_phi", "h_2023C_phi", 80, -4, 4)
h_2023C_etaphi = TH2F("h_2023C_etaphi", "h_2023C_etaphi", 80, -4, 4, 80, -4, 4)

h_2023D_pt    = TH1F("h_2023D_pt", "h_2023D_pt", 100, 0, 50)
h_2023D_eta    = TH1F("h_2023D_eta", "h_2023D_eta", 80, -4, 4)
h_2023D_phi    = TH1F("h_2023D_phi", "h_2023D_phi", 80, -4, 4)
h_2023D_etaphi = TH2F("h_2023D_etaphi", "h_2023D_etaphi", 80, -4, 4, 80, -4, 4)

nev=0
t_scoutMuon_24 = files_2024.Get('scoutingTree/tree')
for ev in t_scoutMuon_24:
    nev+=1
    #if (nev==100): break
    h_2024_pt.Fill(ev.pt1_scout)
    h_2024_pt.Fill(ev.pt2_scout)
    h_2024_eta.Fill(ev.eta1_scout)
    h_2024_eta.Fill(ev.eta2_scout)
    h_2024_phi.Fill(ev.phi1_scout)
    h_2024_phi.Fill(ev.phi2_scout)
    h_2024_etaphi.Fill(ev.eta1_scout,ev.phi1_scout)
    h_2024_etaphi.Fill(ev.eta2_scout,ev.phi2_scout)

nev=0
t_scoutMuon_23C = files_2023C.Get('scoutingTree/tree')
for ev in t_scoutMuon_23C:
    nev+=1
    #if (nev==100): break
    h_2023C_pt.Fill(ev.pt1_scout)
    h_2023C_pt.Fill(ev.pt2_scout)
    h_2023C_eta.Fill(ev.eta1_scout)
    h_2023C_eta.Fill(ev.eta2_scout)
    h_2023C_phi.Fill(ev.phi1_scout)
    h_2023C_phi.Fill(ev.phi2_scout)
    h_2023C_etaphi.Fill(ev.eta1_scout,ev.phi1_scout)
    h_2023C_etaphi.Fill(ev.eta2_scout,ev.phi2_scout)

nev=0
t_scoutMuon_23D = files_2023D.Get('scoutingTree/tree')
for ev in t_scoutMuon_23D:
    nev+=1
    #if (nev==100): break
    h_2023D_pt.Fill(ev.pt1_scout)
    h_2023D_pt.Fill(ev.pt2_scout)
    h_2023D_eta.Fill(ev.eta1_scout)
    h_2023D_eta.Fill(ev.eta2_scout)
    h_2023D_phi.Fill(ev.phi1_scout)
    h_2023D_phi.Fill(ev.phi2_scout)
    h_2023D_etaphi.Fill(ev.eta1_scout,ev.phi1_scout)
    h_2023D_etaphi.Fill(ev.eta2_scout,ev.phi2_scout)


# Scaling
# -------
print("h_2023C_eta.Integral() = ", h_2023C_eta.Integral())
print("h_2023D_eta.Integral() = ", h_2023D_eta.Integral())
print("h_2024_eta.Integral() = ", h_2024_eta.Integral())
h_2023C_eta.Scale(h_2024_eta.Integral()/h_2023C_eta.Integral())
h_2023C_phi.Scale(h_2024_phi.Integral()/h_2023C_phi.Integral())
h_2023C_pt.Scale(h_2024_pt.Integral()/h_2023C_pt.Integral())
h_2023D_eta.Scale(h_2024_eta.Integral()/h_2023D_eta.Integral())
h_2023D_phi.Scale(h_2024_phi.Integral()/h_2023D_phi.Integral())
h_2023D_pt.Scale(h_2024_pt.Integral()/h_2023D_pt.Integral())

# Legend
# -----------------
legend = ROOT.TLegend (0.65, 0.65, 0.88, 0.86)
legend.SetTextSize (0.045)
legend.AddEntry (h_2024_eta, "2024E", "L")
legend.AddEntry (h_2023C_eta, "2023C", "L")
legend.AddEntry (h_2023D_eta, "2023D", "L")
#legend.AddEntry (h_2024_eta, "2024B", "NDC")
#legend.AddEntry (h_2023_eta, "2023D", "NDC")
legend.SetLineWidth (0)

# CMS lumi info 
# -----------------
CMS_lumi.writeExtraText = True
CMS_lumi.extraText      = "Preliminary"
CMS_lumi.lumi_sqrtS      = "2023C-D && 2024E (13.6 TeV)"                                                                                   
CMS_lumi.cmsTextSize    = 0.6
CMS_lumi.lumiTextSize   = 0.46
CMS_lumi.extraOverCmsTextSize = 0.75
CMS_lumi.relPosX = 0.06
#CMS_lumi.relPosX = 0.12


# --------------------------------------------------------------------
#gr_text1 = ROOT.TPaveText(0.12, 0.83, 0.7, 0.87, "NDC")
#gr_text1.AddText("Events with two matched muons, #DeltaR_{#mu#mu} > 0.2 and p_{T}^{#mu} > 3 GeV")
#gr_text1.SetTextSize(0.032)
#gr_text1.SetFillColor(0)

c_eta = ROOT.TCanvas("c_eta", "c_eta", 1000, 1000)
c_eta.cd()    
#c_eta.SetLeftMargin(0.11)
c_eta.SetBottomMargin(0.17)

# Create main pad
pad_eta_main = TPad("pad_eta_main", "pad_eta_main", 0.0, 0.3, 1.0, 1.0)
pad_eta_main.SetBottomMargin(0.02)
pad_eta_main.SetLeftMargin(0.15)
#pad_eta_main.SetLogx()    
#pad_eta_main.SetLogy()    
pad_eta_main.Draw()
pad_eta_main.cd()

h_2024_eta.GetXaxis().SetLabelOffset(0.4)
h_2024_eta.GetYaxis().SetLabelSize(0.045)
#h_2024_eta.GetYaxis().SetLabelOffset(0.01)
h_2024_eta.GetYaxis().SetTitleOffset(1.12)
h_2024_eta.GetYaxis().SetTitleSize(0.06)
h_2024_eta.SetLineWidth(3)
#h_2024_eta.SetLineStyle(8)
h_2024_eta.SetLineColor(kMagenta-7)
h_2024_eta.SetFillColorAlpha(kMagenta-9,0.65)
h_2024_eta.SetFillStyle(3015)
h_2023C_eta.SetLineWidth(3)
h_2023C_eta.SetLineColor(kAzure-4-4)
h_2023D_eta.SetLineWidth(3)
h_2023D_eta.SetLineColor(kOrange-3)

h_2024_eta.SetMaximum(1.4*h_2023C_eta.GetMaximum())
h_2023D_eta.SetMaximum(1.4*h_2023C_eta.GetMaximum())

h_2024_eta.Draw("hist")
h_2023C_eta.Draw("same hist")
h_2023D_eta.Draw("same hist")
h_2024_eta.GetYaxis().SetTitle("A.U.")
#h_2024_eta.GetYaxis().SetTitle("Events / 0.1")
#h_2023C_eta.GetXaxis().SetTitle("Muon #eta")

legend.Draw()

# CMS Info                                                            
# -----------------                                                                                                              
#CMS_lumi.CMS_lumi(pad_eta_main, 0, 0)
latex = TLatex();                                                                                                                            
latex.SetTextSize(0.06);                                                                                                                    
latex.SetTextAlign(13);                                                                                                                         
latex.SetTextFont(62)                                                                                                                       
#latex.DrawLatexNDC(.12,.95,"CMS"); #OutOfFrame                                                                                             
latex.DrawLatexNDC(.18,.86,"CMS");                                                                                                 
latex.SetTextFont(52)                                                                                                                    
#latex.DrawLatexNDC(.20,.95, " Preliminary"); #OutOfFrame                                                                                      
latex.DrawLatexNDC(.27,.86, " Preliminary");                                                                                        
latex.SetTextFont(42)                                                                                                                   
latex.SetTextSize(0.05);                                                                                             
latex.DrawLatexNDC(.59,.95,"2023-2024 (13.6 TeV)");
#latex.DrawLatexNDC(.45,.95,"2023C-D && 2024E (13.6 TeV)");


# Create ratio pad
c_eta.cd()  # Go back to the main canvas
pad_eta_ratio = TPad("pad_eta_ratio", "pad_eta_ratio", 0.0, 0.0, 1.0, 0.3)
pad_eta_ratio.SetTopMargin(0.05)
pad_eta_ratio.SetLeftMargin(0.15)
pad_eta_ratio.SetBottomMargin(0.4)
#pad_eta_ratio.SetLogx()    
pad_eta_ratio.Draw()
pad_eta_ratio.cd()

# Compute and draw ratio histogram
h_ratio = h_2024_eta.Clone()
h_ratio.Divide(h_2023C_eta)
h_ratio.SetLineColor(kBlack)
h_ratio.SetLineWidth(1)
h_ratio.SetMarkerStyle(20)
h_ratio.SetMarkerSize(0.9)

h_ratio.GetXaxis().SetLabelOffset(0.02)
h_ratio.GetXaxis().SetLabelSize(0.09)
h_ratio.GetYaxis().SetLabelSize(0.09)
h_ratio.GetXaxis().SetTitleSize(0.12)
h_ratio.GetYaxis().SetTitleSize(0.11)
h_ratio.GetYaxis().SetRangeUser(0.7,1.3)
h_ratio.GetYaxis().SetNdivisions(303) 
h_ratio.GetXaxis().SetTitle("Scouting muon #eta")
h_ratio.GetXaxis().SetTitleOffset(1.1)
h_ratio.GetYaxis().SetTitle("2024 / 2023C")
h_ratio.GetYaxis().SetTitleOffset(0.5)
h_ratio.Draw("ep")

# Draw a horizontal line at y=1 for reference
line_at_one = TLine(h_ratio.GetXaxis().GetXmin(), 1, h_ratio.GetXaxis().GetXmax(), 1)
line_at_one.SetLineStyle(2)
line_at_one.Draw("same")

# Update and save canvas
c_eta.Update()
#c_eta.SaveAs(outputdir + "/etaCompScoutMu_23vs24.png")
#c_eta.SaveAs(outputdir + "/etaCompScoutMu_23vs24.pdf")
#c_eta.SaveAs(outputdir + "/etaCompScoutMu_23vs24.root")


c_phi = ROOT.TCanvas("c_phi", "c_phi", 1000, 1000)
c_phi.cd()    
#c_phi.SetLeftMargin(0.11)
c_phi.SetBottomMargin(0.4)

# Create main pad
pad_phi_main = TPad("pad_phi_main", "pad_phi_main", 0.0, 0.3, 1.0, 1.0)
pad_phi_main.SetBottomMargin(0.02)
pad_phi_main.SetLeftMargin(0.15)
pad_phi_main.Draw()
pad_phi_main.cd()

h_2024_phi.GetXaxis().SetLabelOffset(0.4)
h_2024_phi.GetYaxis().SetLabelSize(0.045)
h_2024_phi.GetYaxis().SetTitleOffset(1.12)
h_2024_phi.GetYaxis().SetTitleSize(0.06)
h_2024_phi.SetLineWidth(3)
h_2024_phi.SetLineColor(kMagenta-7)
h_2024_phi.SetFillColorAlpha(kMagenta-9,0.65)
h_2024_phi.SetFillStyle(3015)
h_2023C_phi.SetLineWidth(3)
h_2023C_phi.SetLineColor(kAzure-4)
h_2023D_phi.SetLineWidth(3)
h_2023D_phi.SetLineColor(kOrange-3)

print("h_2023C_phi.Integral() = ", h_2023C_phi.Integral())
print("h_2024_phi.Integral() = ", h_2024_phi.Integral())

h_2024_phi.SetMaximum(1.5*h_2023C_phi.GetMaximum())
h_2023D_phi.SetMaximum(1.5*h_2023C_phi.GetMaximum())

h_2024_phi.Draw("hist")
h_2023C_phi.Draw("same hist")
h_2023D_phi.Draw("same hist")
h_2024_phi.GetYaxis().SetTitle("A.U.")
#h_2024_phi.GetYaxis().SetTitle("Events / 0.1")
#h_2023C_phi.GetXaxis().SetTitle("Muon #eta")

legend.Draw()

# CMS Info                                                            
# -----------------                                                                                                              
#CMS_lumi.CMS_lumi(pad_phi_main, 0, 0)
latex = TLatex();                                                                                                                            
latex.SetTextSize(0.06);                                                                                                                    
latex.SetTextAlign(13);                                                                                                                         
latex.SetTextFont(62)                                                                                                                       
#latex.DrawLatexNDC(.12,.95,"CMS"); #OutOfFrame                                                                                             
latex.DrawLatexNDC(.18,.86,"CMS");                                                                                                 
latex.SetTextFont(52)                                                                                                                    
#latex.DrawLatexNDC(.20,.95, " Preliminary"); #OutOfFrame                                                                                      
latex.DrawLatexNDC(.27,.86, " Preliminary");                                                                                        
latex.SetTextFont(42)                                                                                                                   
latex.SetTextSize(0.05);                                                                                             
latex.DrawLatexNDC(.59,.95,"2023-2024 (13.6 TeV)");
#latex.DrawLatexNDC(.45,.95,"2023C-D && 2024E (13.6 TeV)");


# Create ratio pad
c_phi.cd()  # Go back to the main canvas
pad_phi_ratio = TPad("pad_phi_ratio", "pad_phi_ratio", 0.0, 0.0, 1.0, 0.3)
pad_phi_ratio.SetTopMargin(0.05)
pad_phi_ratio.SetLeftMargin(0.15)
pad_phi_ratio.SetBottomMargin(0.4)
pad_phi_ratio.Draw()
pad_phi_ratio.cd()

# Compute and draw ratio histogram
h_ratio = h_2024_phi.Clone()
h_ratio.Divide(h_2023C_phi)
h_ratio.SetLineColor(kBlack)
h_ratio.SetLineWidth(1)
h_ratio.SetMarkerStyle(20)
h_ratio.SetMarkerSize(0.9)

h_ratio.GetXaxis().SetLabelOffset(0.02)
h_ratio.GetXaxis().SetLabelSize(0.09)
h_ratio.GetYaxis().SetTitleSize(0.12)
h_ratio.GetYaxis().SetLabelSize(0.09)
h_ratio.GetXaxis().SetTitleSize(0.12)
h_ratio.GetYaxis().SetTitleSize(0.11)
h_ratio.GetYaxis().SetRangeUser(0.7,1.3)
h_ratio.GetYaxis().SetNdivisions(303) 
h_ratio.GetXaxis().SetTitle("Scouting muon #phi")
h_ratio.GetXaxis().SetTitleOffset(1.1)
h_ratio.GetYaxis().SetTitle("2024 / 2023C")
h_ratio.GetYaxis().SetTitleOffset(0.5)
h_ratio.Draw("ep")

# Draw a horizontal line at y=1 for reference
line_at_one = TLine(h_ratio.GetXaxis().GetXmin(), 1, h_ratio.GetXaxis().GetXmax(), 1)
line_at_one.SetLineStyle(2)
line_at_one.Draw("same")

# Update and save canvas
c_phi.Update()
c_phi.SaveAs(outputdir + "/phiCompScoutMu_23vs24.png")
c_phi.SaveAs(outputdir + "/phiCompScoutMu_23vs24.pdf")
c_phi.SaveAs(outputdir + "/phiCompScoutMu_23vs24.root")


c_pt = ROOT.TCanvas("c_pt", "c_pt", 1000, 1000)
c_pt.cd()    
#c_pt.SetLeftMargin(0.11)
c_pt.SetBottomMargin(0.17)

# Create main pad
pad_pt_main = TPad("pad_pt_main", "pad_pt_main", 0.0, 0.3, 1.0, 1.0)
pad_pt_main.SetLeftMargin(0.15)
pad_pt_main.SetBottomMargin(0.02)
#pad_pt_main.SetLogx()    
#pad_pt_main.SetLogy()    
pad_pt_main.Draw()
pad_pt_main.cd()

h_2024_pt.GetXaxis().SetLabelOffset(0.4)
h_2024_pt.GetYaxis().SetLabelSize(0.045)
#h_2024_pt.GetYaxis().SetLabelOffset(0.01)
h_2024_pt.GetYaxis().SetTitleOffset(1.12)
h_2024_pt.GetYaxis().SetTitleSize(0.06)
h_2024_pt.SetLineWidth(3)
#h_2024_pt.SetLineStyle(8)
h_2024_pt.SetLineColor(kMagenta-7)
h_2024_pt.SetFillColorAlpha(kMagenta-9,0.65)
h_2024_pt.SetFillStyle(3015)
h_2023C_pt.SetLineWidth(3)
h_2023C_pt.SetLineColor(kAzure-4)
h_2023D_pt.SetLineWidth(3)
h_2023D_pt.SetLineColor(kOrange-3)

# Scaling
print("h_2023C_pt.Integral() = ", h_2023C_pt.Integral())
print("h_2024_pt.Integral() = ", h_2024_pt.Integral())

h_2024_pt.SetMaximum(1.5*h_2023C_pt.GetMaximum())
h_2023D_pt.SetMaximum(1.5*h_2023C_pt.GetMaximum())

h_2024_pt.Draw("hist")
h_2023C_pt.Draw("same hist")
h_2023D_pt.Draw("same hist")
h_2024_pt.GetYaxis().SetTitle("A.U.")
#h_2024_pt.GetYaxis().SetTitle("Events / 0.1")
#h_2023C_pt.GetXaxis().SetTitle("Muon #pt")

legend.Draw()

# CMS Info                                                            
# -----------------                                                                                                              
#CMS_lumi.CMS_lumi(pad_pt_main, 0, 0)
latex = TLatex();                                                                                                                            
latex.SetTextSize(0.06);                                                                                                                    
latex.SetTextAlign(13);                                                                                                                         
latex.SetTextFont(62)                                                                                                                       
#latex.DrawLatexNDC(.12,.95,"CMS"); #OutOfFrame                                                                                             
latex.DrawLatexNDC(.18,.86,"CMS");                                                                                                 
latex.SetTextFont(52)                                                                                                                    
#latex.DrawLatexNDC(.20,.95, " Preliminary"); #OutOfFrame                                                                                      
latex.DrawLatexNDC(.27,.86, " Preliminary");                                                                                        
latex.SetTextFont(42)                                                                                                                   
latex.SetTextSize(0.05);                                                                                             
latex.DrawLatexNDC(.595,.95,"2023-2024 (13.6 TeV)");
#latex.DrawLatexNDC(.45,.95,"2023C-D && 2024E (13.6 TeV)");


# Create ratio pad
c_pt.cd()  # Go back to the main canvas
pad_pt_ratio = TPad("pad_pt_ratio", "pad_pt_ratio", 0.0, 0.0, 1.0, 0.3)
pad_pt_ratio.SetTopMargin(0.05)
pad_pt_ratio.SetLeftMargin(0.15)
pad_pt_ratio.SetBottomMargin(0.3)
#pad_pt_ratio.SetLogx()    
pad_pt_ratio.Draw()
pad_pt_ratio.cd()

# Compute and draw ratio histogram
h_ratio = h_2024_pt.Clone()
h_ratio.Divide(h_2023C_pt)
h_ratio.SetLineColor(kBlack)
h_ratio.SetLineWidth(1)
h_ratio.SetMarkerStyle(20)
h_ratio.SetMarkerSize(0.9)

h_ratio.GetXaxis().SetLabelOffset(0.02)
h_ratio.GetXaxis().SetLabelSize(0.08)
h_ratio.GetYaxis().SetTitleSize(0.12)
h_ratio.GetYaxis().SetLabelSize(0.07)
h_ratio.GetXaxis().SetTitleSize(0.12)
h_ratio.GetYaxis().SetTitleSize(0.11)
h_ratio.GetYaxis().SetRangeUser(0.7,1.3)
h_ratio.GetYaxis().SetNdivisions(505) 
h_ratio.GetXaxis().SetTitle("Scouting muon p_{T}")
h_ratio.GetXaxis().SetTitleOffset(1.1)
h_ratio.GetYaxis().SetTitle("2024 / 2023C")
h_ratio.GetYaxis().SetTitleOffset(0.5)
h_ratio.Draw("ep")

# Draw a horizontal line at y=1 for reference
line_at_one = TLine(h_ratio.GetXaxis().GetXmin(), 1, h_ratio.GetXaxis().GetXmax(), 1)
line_at_one.SetLineStyle(2)
line_at_one.Draw("same")

# Update and save canvas
c_pt.Update()
#c_pt.SaveAs(outputdir + "/ptCompScoutMu_23vs24.png")
#c_pt.SaveAs(outputdir + "/ptCompScoutMu_23vs24.pdf")
#c_pt.SaveAs(outputdir + "/ptCompScoutMu_23vs24.root")

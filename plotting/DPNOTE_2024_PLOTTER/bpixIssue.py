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
argparser.add_argument('--log', dest='log', default='False', help='Output directory')
args = argparser.parse_args()
outputdir = args.outdir
log = args.log

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

######################################
# List of files and output directory #
######################################
def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]
f_23C = ROOT.TFile.Open("/eos/user/e/elfontan/2023_SCOUTING/ScoutingPFMonitor/scoutMonitor_2023Cv1.root")
f_23D = ROOT.TFile.Open("/eos/user/e/elfontan/2023_SCOUTING/ScoutingPFMonitor/scoutMonitor_2023Dv1.root")
#files = files[:-5] 


########################
# Variables and histos #
########################
h_23C_etaphi = TH2F("h_23C_etaphi", "h_23C_etaphi", 60, -3, 3, 80, -4, 4)
h_23D_etaphi = TH2F("h_23D_etaphi", "h_23D_etaphi", 60, -3, 3, 80, -4, 4)

t_scoutMuon_23C = f_23C.Get('scoutingTree/tree')
t_scoutMuon_23D = f_23D.Get('scoutingTree/tree')

# Check if the TTree was extracted successfully
if not t_scoutMuon_23C:
    print('Error: failed to extract TTree from {f_23C}')
    f_23C.Close()

print("Running on 23C tree...")
nEv=0
for ev in t_scoutMuon_23C:
    nEv+=1
    #if (nEv==10): break
    h_23C_etaphi.Fill(ev.eta1,ev.phi1_scout)
    h_23C_etaphi.Fill(ev.eta2,ev.phi2_scout)

# Check if the TTree was extracted successfully
if not t_scoutMuon_23D:
    print('Error: failed to extract TTree from {f_23D}')
    f_23D.Close()

print("Running on 23D tree...")
nEv=0
for ev in t_scoutMuon_23D:
    nEv+=1
    #if (nEv==10): break
    h_23D_etaphi.Fill(ev.eta1,ev.phi1_scout)
    h_23D_etaphi.Fill(ev.eta2,ev.phi2_scout)

    
ROOT.gStyle.SetPalette(112)

################
### 2D 2023C ###
################
c_23C = ROOT.TCanvas("c_23C", "c_23C", 1200, 900)
c_23C.cd()    
c_23C.SetRightMargin(0.19)
c_23C.SetBottomMargin(0.15)

h_23C_etaphi.GetXaxis().SetLabelSize(0.05)
h_23C_etaphi.GetYaxis().SetLabelSize(0.05)
h_23C_etaphi.GetZaxis().SetLabelSize(0.05)
h_23C_etaphi.GetXaxis().SetTitleOffset(1.1)
h_23C_etaphi.GetYaxis().SetTitleOffset(0.9)
h_23C_etaphi.GetZaxis().SetTitleOffset(1.1)
h_23C_etaphi.GetXaxis().SetTitleSize(0.055)
h_23C_etaphi.GetYaxis().SetTitleSize(0.055)
h_23C_etaphi.GetZaxis().SetTitleSize(0.055)
h_23C_etaphi.GetXaxis().SetTitle("Scouting muon #eta")
h_23C_etaphi.GetYaxis().SetTitle("Scouting muon #phi")
h_23C_etaphi.GetZaxis().SetTitle("Number of events")
h_23C_etaphi.Draw("colz")

latex = TLatex();                                                                                                                  
latex.SetTextSize(0.055);                                                                                                         
latex.SetTextAlign(13);                                                                                                                
latex.SetTextFont(62)                                                                                                                
latex.DrawLatexNDC(.105,.95,"CMS"); #OutOfFrame                                                                                     
#latex.DrawLatexNDC(.13,.86,"CMS");                                                                                                 
latex.SetTextFont(52)                                                                                                             
latex.DrawLatexNDC(.19,.95, " Preliminary"); #OutOfFrame                                                                             
#latex.DrawLatexNDC(.20,.86, " Preliminary");                                                                                    
latex.SetTextFont(42)                                        
latex.SetTextSize(0.05);                                                                                
latex.DrawLatexNDC(.6,.95,"2023C (13.6 TeV)");


c_23C.SaveAs(outputdir + "/etaphi_2023C.png")
c_23C.SaveAs(outputdir + "/etaphi_2023C.pdf")

################
### 2D 2023D ###
################
c_23D = ROOT.TCanvas("c_23D", "c_23D", 1200, 900)
c_23D.cd()    
c_23D.SetRightMargin(0.19)
c_23D.SetBottomMargin(0.15)

h_23D_etaphi.GetXaxis().SetLabelSize(0.05)
h_23D_etaphi.GetYaxis().SetLabelSize(0.05)
h_23D_etaphi.GetZaxis().SetLabelSize(0.05)
h_23D_etaphi.GetXaxis().SetTitleOffset(1.1)
h_23D_etaphi.GetYaxis().SetTitleOffset(0.9)
h_23D_etaphi.GetZaxis().SetTitleOffset(1.1)
h_23D_etaphi.GetXaxis().SetTitleSize(0.055)
h_23D_etaphi.GetYaxis().SetTitleSize(0.055)
h_23D_etaphi.GetZaxis().SetTitleSize(0.055)
h_23D_etaphi.GetXaxis().SetTitle("Scouting muon #eta")
h_23D_etaphi.GetYaxis().SetTitle("Scouting muon #phi")
h_23D_etaphi.GetZaxis().SetTitle("Number of events")
h_23D_etaphi.Draw("colz")

latex = TLatex();                                                                                                                  
latex.SetTextSize(0.055);                                                                                                         
latex.SetTextAlign(13);                                                                                                                
latex.SetTextFont(62)                                                                                                                
latex.DrawLatexNDC(.105,.95,"CMS"); #OutOfFrame                                                                                     
latex.SetTextFont(52)                                                                                                             
latex.DrawLatexNDC(.19,.95, " Preliminary"); #OutOfFrame                                                                             
#latex.DrawLatexNDC(.20,.86, " Preliminary");                                                                                    
latex.SetTextFont(42)                                        
latex.SetTextSize(0.05);                                                                                
latex.DrawLatexNDC(.6,.95,"2023D (13.6 TeV)");


c_23D.SaveAs(outputdir + "/etaphi_2023D.png")
c_23D.SaveAs(outputdir + "/etaphi_2023D.pdf")

#!/usr/bin/env python3
# ***************************************
# usage: 
#    python3 scoutVSoffline.py
# ***************************************

import ROOT, array, random, copy
from ROOT import TCanvas, TFile, TH1, TH1F, TF1, gSystem
#from ROOT import *
import ROOT, array, CMSGraphics, CMS_lumi, random, copy
from ROOT import RooCmdArg, RooArgSet, kFALSE, RooLinkedList, kBlue, kRed, kBlack, kOpenStar, kWhite, kGray
from ROOT import gStyle, TStyle, TGraph, TGraphErrors, TMath, TMultiGraph, TLine, gPad, TGaxis, TLegend, TText, TLatex, TColor, TPaveText
from ROOT import TAttFill, TLegend, TRatioPlot, TPad, THStack, TFileCollection
from ROOT import kBlue, kRed, kBlack, kWhite, kAzure, kOrange, kPink, kGreen, kYellow, kCyan, kMagenta
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

#lowmass = False
#fullmass = True
lowmass = True
fullmass = False

vtx = False
novtx = True

######################################
# List of files and output directory #
######################################
def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]
#files = list_full_paths("/eos/user/e/elfontan/2024_SCOUTING/vtxMu_ScoutingPFMonitor/")
#files = files[-1:] 

if (vtx):
    files = [
        "/eos/user/e/elfontan/2024_SCOUTING/vtxMu_ScoutingPFMonitor/scoutMonitor_2024D_vtxMu.root",
        "/eos/user/e/elfontan/2024_SCOUTING/vtxMu_ScoutingPFMonitor/scoutMonitor_2024E_vtxMu.root",
        "/eos/user/e/elfontan/2024_SCOUTING/vtxMu_ScoutingPFMonitor/scoutMonitor_2024F_vtxMu.root",
        "/eos/user/e/elfontan/2024_SCOUTING/vtxMu_ScoutingPFMonitor/scoutMonitor_2024G_vtxMu.root",
        "/eos/user/e/elfontan/2024_SCOUTING/vtxMu_ScoutingPFMonitor/scoutMonitor_2024H_vtxMu.root",
        "/eos/user/e/elfontan/2024_SCOUTING/vtxMu_ScoutingPFMonitor/scoutMonitor_2024I_vtxMu.root",
    ]
elif (novtx):
    files = [
        "/eos/user/e/elfontan/2024_SCOUTING/noVtxMu_ScoutingPFMonitor/scoutMonitor_2024D_noVtxMu.root",
        "/eos/user/e/elfontan/2024_SCOUTING/noVtxMu_ScoutingPFMonitor/scoutMonitor_2024E_noVtxMu.root",
        "/eos/user/e/elfontan/2024_SCOUTING/noVtxMu_ScoutingPFMonitor/scoutMonitor_2024F_noVtxMu.root",
        "/eos/user/e/elfontan/2024_SCOUTING/noVtxMu_ScoutingPFMonitor/scoutMonitor_2024G_noVtxMu.root",
        "/eos/user/e/elfontan/2024_SCOUTING/noVtxMu_ScoutingPFMonitor/scoutMonitor_2024H_noVtxMu.root",
        "/eos/user/e/elfontan/2024_SCOUTING/noVtxMu_ScoutingPFMonitor/scoutMonitor_2024I_noVtxMu.root",
    ]


########################
# Variables and histos #
########################
h_pt_res_zoom   = TH1F("h_pt_res_zoom", "h_pt_res_zoom", 200, -0.1, 0.1)
h_mass_res_zoom = TH1F("h_mass_res_zoom", "h_mass_res_zoom", 200, -0.1, 0.1)
h_pt_res        = TH1F("h_pt_res", "h_pt_res", 400, -0.4, 0.4) 
h_mass_res      = TH1F("h_mass_res", "h_mass_res", 400, -0.4, 0.4) 

if (fullmass):
    xbins = [0.215]
    while (xbins[-1]<250):
        xbins.append(1.01*xbins[-1])
    #print("xbins", xbins)
    xbins_rebin = [0.215]
    while (xbins_rebin[-1]<250):
        xbins_rebin.append(1.05*xbins_rebin[-1])
if (lowmass):
    xbins = [0.215]
    while (xbins[-1]<20):
        xbins.append(1.01*xbins[-1])
    print("xbins", xbins)
    #print("len(xbins)", len(xbins))
    xbins_rebin = [0.215]
    while (xbins_rebin[-1]<20):
        xbins_rebin.append(1.05*xbins_rebin[-1])

h_mass_offline = TH1F("h_mass_offline", "h_mass_offline", len(xbins)-1,array('f',xbins)) 
h_mass_scout = TH1F("h_mass_scout", "h_mass_scout", len(xbins)-1,array('f',xbins)) 
h_mass_offline_reb = TH1F("h_mass_offline_reb", "h_mass_offline_reb", len(xbins_rebin)-1,array('f',xbins_rebin))
h_mass_scout_reb = TH1F("h_mass_scout_reb", "h_mass_scout_reb", len(xbins_rebin)-1,array('f',xbins_rebin))

if (fullmass):
    frame = TH1F("frame","",1000,0.18,250.2)
    f_ratio = TH1F("f_ratio","",1000,0.18,250.2)
elif (lowmass):
    frame = TH1F("frame","",1000,0.18,20.2)
    f_ratio = TH1F("f_ratio","",1000,0.18,20.2)

# Loop over the files and fill the histo
# --------------------------------------
print(">>>>>> READING...")
print(">>>>>> List of files:")
for filename in files:
    root_file = ROOT.TFile.Open(filename)
    print(filename)

    t_scoutMuon = root_file.Get('scoutingTree/tree')

    if not t_scoutMuon:
        print('Error: failed to extract TTree from {file_name}')
        root_file.Close()
        continue

    for ev in t_scoutMuon:
      #print("nScoutingMuons = ", ev.nScoutingMuons )
      if (not(ev.nScoutingMuons == 2)): continue
      #if (ev.pt1_scout == ev.pt2_scout): continue
      #if (ev.drmm < 0.2): continue 
      #if (ev.drmm < 0.2 or ev.drmm_scout < 0.2): continue 
      if (ev.dr_matching_1 > 0.2 or ev.dr_matching_2 > 0.2): continue

      if (ev.ndvtx > 0 and ev.pt1_scout > 3 and ev.pt2_scout > 3):
          h_pt_res_zoom.Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
          h_pt_res_zoom.Fill((ev.pt2_scout - ev.pt2)/ev.pt2)
          h_pt_res.Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
          h_pt_res.Fill((ev.pt2_scout - ev.pt2)/ev.pt2)
          h_mass_res_zoom.Fill((ev.mass_scout - ev.mass)/ev.mass)
          h_mass_res.Fill((ev.mass_scout - ev.mass)/ev.mass)
          h_mass_offline.Fill(ev.mass)
          h_mass_scout.Fill(ev.mass_scout)
          h_mass_offline_reb.Fill(ev.mass)
          h_mass_scout_reb.Fill(ev.mass_scout)
      else:
          continue

    root_file.Close()


print("h_mass_offline.Integral() = ", h_mass_offline.Integral())
print("h_mass_scout.Integral() = ", h_mass_scout.Integral())
legend = ROOT.TLegend (0.6, 0.6, 0.86, 0.86)
legend.SetTextSize (0.03)
legend.AddEntry (h_pt_res, "Uncorrected muons", "NDC")
legend.SetLineWidth (0)

gr_text1 = ROOT.TPaveText(0.11, 0.73, 0.7, 0.76, "NDC")
gr_text1.AddText("Events with two matched muons, #DeltaR_{#mu#mu} > 0.2 and p_{T}^{#mu} > 3 GeV")
gr_text1.SetTextSize(0.032)
gr_text1.SetFillColor(0)

leg_mass = ROOT.TLegend (0.7, 0.62, 0.85, 0.79)
leg_mass.SetTextSize(0.037)
leg_mass.AddEntry (h_mass_offline, "Offline", "F")
leg_mass.AddEntry (h_mass_scout, "Scouting", "F")
leg_mass.SetLineWidth (0)

labels = TLatex()
masses = {
    "#bf{#eta}": 0.546862,
    "#bf{#rho,#omega}": 0.780,
    "#bf{#phi}": 1.019,
    "#bf{J/#Psi}": 3.096,
    "#bf{#Psi'}": 3.686,
    "#bf{#Upsilon(nS)}": 9.460,
    "#bf{Z}": 91.1876,
}
labels.SetTextSize(0.04)
labels.SetTextAlign(21)

if (fullmass):
    c_mass = ROOT.TCanvas("c_mass", "c_mass", 1200, 1000)
    c_mass.cd()    
    c_mass.SetLeftMargin(0.13)
    c_mass.SetBottomMargin(0.17)
    
    pad_main = TPad("pad_main", "pad_main", 0.0, 0.3, 1.0, 1.0)
    pad_main.SetBottomMargin(0.02)
    pad_main.SetLogx()    
    pad_main.SetLogy()    
    pad_main.Draw()
    pad_main.cd()
    
    frame.SetMinimum(10)
    frame.SetMaximum(1000000000)
    frame.GetXaxis().SetLabelOffset(0.2)
    frame.GetXaxis().SetTitleOffset(2.1)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    frame.GetYaxis().SetTitle("Events / MeV")

    frame.Draw()

    h_mass_offline.SetLineWidth(2)
    h_mass_offline.SetLineColor(kBlue-3)
    h_mass_scout.SetLineWidth(3)
    h_mass_scout.SetFillColor(kMagenta-9)
    h_mass_scout.SetLineColor(kMagenta-9)
    h_mass_scout.SetFillStyle(3015)    
    h_mass_scout.Scale(1., "width")
    h_mass_offline.Scale(1., "width")
    h_mass_scout.Draw("same hist")
    h_mass_offline.Draw("same hist")
    frame.Draw("same axis")
    
    [ labels.DrawLatex(masses[m], 1.5*h_mass_offline.GetBinContent(h_mass_scout.FindBin(masses[m])), m) for m in masses ]
    labels.Draw("same")
    gr_text1.Draw("same")
    leg_mass.Draw ("same")

    latex = TLatex();                                                                                                                         
    latex.SetTextSize(0.055);                                                                                                                       
    latex.SetTextAlign(13);                                                                                                                         
    latex.SetTextFont(62)                                                                                                                
    latex.DrawLatexNDC(.13,.86,"CMS");                                                                                               
    latex.SetTextFont(52)                                                                                                            
    latex.DrawLatexNDC(.19,.86, " Preliminary");                                                                             
    latex.SetTextFont(42)                                                                                                   
    latex.SetTextSize(0.055);                                                                                         
    latex.DrawLatexNDC(.7,.96,"2024 (13.6 TeV)");

    
    # Create ratio pad
    # ----------------
    c_mass.cd()  
    pad_ratio = TPad("pad_ratio", "pad_ratio", 0.0, 0.0, 1.0, 0.3)
    pad_ratio.SetTopMargin(0.05)
    pad_ratio.SetBottomMargin(0.35)
    pad_ratio.SetTicks()    
    pad_ratio.SetLogx()    
    pad_ratio.Draw()
    pad_ratio.cd()
        
    # Compute and draw ratio histogram                     
    # --------------------------------
    h_mass_offline_reb.Scale(1., "width")
    h_mass_scout_reb.Scale(1., "width")
    h_ratio = h_mass_offline_reb.Clone()
    h_ratio.Divide(h_mass_scout_reb)                                                                                    
    
    h_ratio.SetFillColor(kGray)
    f_ratio.GetYaxis().SetRangeUser(0.8, 1.2)
    f_ratio.GetXaxis().SetLabelSize(0.12)
    f_ratio.GetYaxis().SetLabelSize(0.09)
    f_ratio.GetXaxis().SetTitleSize(0.12)
    f_ratio.GetYaxis().SetTitleSize(0.1)
    f_ratio.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    f_ratio.GetXaxis().SetTitleOffset(1.3)
    f_ratio.GetYaxis().SetTitle("Offline / Scouting")
    f_ratio.GetYaxis().SetTitleOffset(0.45)
    h_ratio.SetBinContent(0, 0)

    f_ratio.GetYaxis().SetNdivisions(505)
    f_ratio.GetXaxis().SetTickSize(0.06)

    f_ratio.Draw("")
    h_ratio.Draw("E4 same")

    line_at_one = TLine(0.215, 1, 250.1, 1)
    line_at_one.SetLineStyle(2)
    line_at_one.Draw("same")

    textRatio = TLatex();                                                                                                                         
    textRatio.SetTextSize(0.08);                                                                                                                       
    textRatio.SetTextAlign(13);                                                                                                                         
    textRatio.SetTextFont(62)                                                                                                                
    if (vtx):
        textRatio.DrawLatexNDC(.125,.88,"Scouting Vtx muon reconstruction");                                                                                               
    elif (novtx):
        textRatio.DrawLatexNDC(.125,.88,"Scouting NoVtx muon reconstruction");                                                                                               
    
    # Update and save canvas
    # ----------------------
    c_mass.Update()
    if (vtx):
        c_mass.SaveAs(outputdir + "/dimuScoutingVtx_fullmass_2024.png")
        c_mass.SaveAs(outputdir + "/dimuScoutingVtx_fullmass_2024.C")
        c_mass.SaveAs(outputdir + "/dimuScoutingVtx_fullmass_2024.pdf")
    elif (novtx):
        c_mass.SaveAs(outputdir + "/dimuScoutingNoVtx_fullmass_2024.png")
        c_mass.SaveAs(outputdir + "/dimuScoutingNoVtx_fullmass_2024.C")
        c_mass.SaveAs(outputdir + "/dimuScoutingNoVtx_fullmass_2024.pdf")
    
if (lowmass):
    c_lowmass = ROOT.TCanvas("c_lowmass", "c_lowmass", 1200, 1000)
    c_lowmass.cd()    
    c_lowmass.SetLeftMargin(0.13)
    c_lowmass.SetBottomMargin(0.17)

    # Create main pad
    # ---------------
    pad_main = TPad("pad_main", "pad_main", 0.0, 0.3, 1.0, 1.0)
    pad_main.SetBottomMargin(0.02)
    pad_main.SetLogx()    
    pad_main.SetLogy()    
    pad_main.SetTicks()    
    pad_main.Draw()
    pad_main.cd()

    frame.SetMinimum(500)
    frame.SetMaximum(180000000)
    frame.GetXaxis().SetLabelOffset(0.2)
    frame.GetXaxis().SetTitleOffset(2.1)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    frame.GetYaxis().SetTitle("Events / MeV")
    frame.Draw()

    h_mass_offline.SetLineWidth(2)
    h_mass_offline.SetLineColor(kBlue-3)
    h_mass_scout.SetLineWidth(3)
    h_mass_scout.SetFillColor(kMagenta-9)
    h_mass_scout.SetLineColor(kMagenta-9)
    h_mass_scout.SetFillStyle(3015)    
    h_mass_scout.Scale(1., "width")
    h_mass_offline.Scale(1., "width")
    h_mass_scout.Draw("same hist")
    h_mass_offline.Draw("same hist")
    frame.Draw("same axis")
    
    [ labels.DrawLatex(masses[m], 1.5*h_mass_offline.GetBinContent(h_mass_scout.FindBin(masses[m])), m) for m in masses ]
    labels.Draw("same")
    gr_text1.Draw("same")
    leg_mass.Draw ("same")
    latex = TLatex();                                                                                                                         
    latex.SetTextSize(0.055);                                                                                                                       
    latex.SetTextAlign(13);                                                                                                                         
    latex.SetTextFont(62)                                                                                                                
    latex.DrawLatexNDC(.13,.86,"CMS");                                                                                               
    latex.SetTextFont(52)                                                                                                            
    latex.DrawLatexNDC(.19,.86, " Preliminary");                                                                             
    latex.SetTextFont(42)                                                                                                   
    latex.SetTextSize(0.055);                                                                                         
    latex.DrawLatexNDC(.7,.96,"2024 (13.6 TeV)");

    # Create ratio pad
    # ----------------
    c_lowmass.cd()  
    pad_ratio = TPad("pad_ratio", "pad_ratio", 0.0, 0.0, 1.0, 0.3)
    pad_ratio.SetTopMargin(0.05)
    pad_ratio.SetBottomMargin(0.35)
    pad_ratio.SetTicks()    
    pad_ratio.SetLogx()    
    pad_ratio.Draw()
    pad_ratio.cd()
        
    # Compute and draw ratio histogram                     
    # --------------------------------
    h_mass_offline_reb.Scale(1., "width")
    h_mass_scout_reb.Scale(1., "width")
    h_ratio = h_mass_offline_reb.Clone()
    h_ratio.Divide(h_mass_scout_reb)                                                                                    
    
    h_ratio.SetFillColor(kGray)
    f_ratio.GetYaxis().SetRangeUser(0.8, 1.2)
    f_ratio.GetXaxis().SetLabelSize(0.12)
    f_ratio.GetYaxis().SetLabelSize(0.09)
    f_ratio.GetXaxis().SetTitleSize(0.12)
    f_ratio.GetYaxis().SetTitleSize(0.1)
    f_ratio.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    f_ratio.GetXaxis().SetTitleOffset(1.3)
    f_ratio.GetYaxis().SetTitle("Offline / Scouting")
    f_ratio.GetYaxis().SetTitleOffset(0.45)
    h_ratio.SetBinContent(0, 0)

    f_ratio.GetYaxis().SetNdivisions(505)
    f_ratio.GetXaxis().SetTickSize(0.06)

    f_ratio.Draw("")
    h_ratio.Draw("E4 same")

    line_at_one = TLine(0.215, 1, 20.1, 1)
    line_at_one.SetLineStyle(2)
    line_at_one.Draw("same")

    textRatio = TLatex();                                                                                                                         
    textRatio.SetTextSize(0.08);                                                                                                                       
    textRatio.SetTextAlign(13);                                                                                                                         
    textRatio.SetTextFont(62)                                                                                                                
    if (vtx):
        textRatio.DrawLatexNDC(.125,.88,"Scouting Vtx muon reconstruction");                                                                                               
    elif (novtx):
        textRatio.DrawLatexNDC(.125,.88,"Scouting NoVtx muon reconstruction");                                                                                               
    
    c_lowmass.Update()
    if (vtx):
        c_lowmass.SaveAs(outputdir + "/dimuScoutingVtx_lowmass_2024.png")
        c_lowmass.SaveAs(outputdir + "/dimuScoutingVtx_lowmass_2024.pdf")
        c_lowmass.SaveAs(outputdir + "/dimuScoutingVtx_lowmass_2024.C")
    elif (novtx):
        c_lowmass.SaveAs(outputdir + "/dimuScoutingNoVtx_lowmass_2024.png")
        c_lowmass.SaveAs(outputdir + "/dimuScoutingNoVtx_lowmass_2024.pdf")
        c_lowmass.SaveAs(outputdir + "/dimuScoutingNoVtx_lowmass_2024.C")
    

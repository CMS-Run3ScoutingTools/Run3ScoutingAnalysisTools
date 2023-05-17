#!/usr/bin/env python3

# *****************************
# usage: 
#    python3 scoutVSoffline.py
# *****************************

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

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

outputdir = "/eos/user/e/elfontan/www/ScoutingParkingPaper/2022_scoutMon_dimuonTriggerSel"

########################
# Variables and histos #
########################
h_pt_res = TH1F("h_pt_res", "h_pt_res", 200, -0.1, 0.1)
h_mass_res = TH1F("h_mass_res", "h_mass_res", 200, -0.1, 0.1)
h_mass_offline = TH1F("h_mass_offline", "h_mass_offline", 280, 0., 140)
h_mass_scout = TH1F("h_mass_scout", "h_mass_scout", 280, 0., 140)
#h_pt_res = TH1F("h_pt_res", "h_pt_res", 400, -0.4, 0.4) #LOG
#h_mass_res = TH1F("h_mass_res", "h_mass_res", 400, -0.4, 0.4) #LOG
#h_mass_offline = TH1F("h_mass_offline", "h_mass_offline", 500, 1., 251) #LOG
#h_mass_scout = TH1F("h_mass_scout", "h_mass_scout", 500, 1., 251) #LOG

def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]
files = list_full_paths("/eos/user/e/elfontan/ScoutingParkingPaper/2022_scoutMon_dimuonTriggerSel")
files = files[:-1] 

# Loop over the files and fill the histo
print(">>>>>> READING...")
print(">>>>>> List of files:")
for filename in files:
    root_file = ROOT.TFile.Open(filename)
    print(filename)

    # Extract the TTree from the ROOT file
    t_scoutMuon = root_file.Get('scoutingTree/tree')

    # Check if the TTree was extracted successfully
    if not t_scoutMuon:
        print('Error: failed to extract TTree from {file_name}')
        root_file.Close()
        continue

    #t_scoutMuon.Draw("pt1>>h_pt_res")
    for ev in t_scoutMuon:
        h_pt_res.Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
        h_mass_res.Fill((ev.mass_scout - ev.mass)/ev.mass)
        h_mass_offline.Fill(ev.mass)
        h_mass_scout.Fill(ev.mass_scout)

    # Close the ROOT file
    root_file.Close()


legend = ROOT.TLegend (0.6, 0.6, 0.86, 0.86)
legend.SetTextSize (0.03)
legend.AddEntry (h_pt_res, "Uncorrected muons", "L")
legend.SetLineWidth (0)

# Create a title for the legend
#title = ROOT.TPaveText(0.65, 0.74, 0.85, 0.92, "NDC")
#title.AddText("Scouting VS Offline")
#title.SetFillColor(0)
#title.SetTextAlign(12)
#title.SetTextSize(0.04)

#draw CMS and lumi text                                                                                                                                           
CMS_lumi.writeExtraText = True
CMS_lumi.lumi_sqrtS      = "2022 (13.6 TeV)"                                                                                                                 
#CMS_lumi.lumi_sqrtS     = lumiText + " (13.6 TeV)"                                                                                                                 
CMS_lumi.cmsTextSize    = 0.6
CMS_lumi.lumiTextSize   = 0.46
CMS_lumi.extraOverCmsTextSize = 0.75
CMS_lumi.relPosX = 0.12

c_pt_res = ROOT.TCanvas("c_pt_res", "c_pt_res", 1000, 800)
c_pt_res.cd()    
#c_pt_res.SetLogy()    
c_pt_res.SetBottomMargin(0.17)
h_pt_res.GetXaxis().SetLabelOffset(0.02)
h_pt_res.GetXaxis().SetTitleOffset(1.9)
#h_pt_res.GetXaxis().SetTitle("(p_{T}^{scout}-p_{T}^{off})/p_{T}^{off} [GeV]")
h_pt_res.GetXaxis().SetTitle("#frac{p_{T}^{scout}-p_{T}^{off}}{p_{T}^{off}} [GeV]")
h_pt_res.GetYaxis().SetTitleOffset(1.5)
h_pt_res.GetYaxis().SetTitle("A.U.")
#h_pt_res.GetYaxis().SetTitle("0.1 GeV/Event")
h_pt_res.SetLineWidth(2)
h_pt_res.SetLineColor(kAzure-4)
h_pt_res.SetFillColorAlpha(kAzure-9,0.35)
h_pt_res.SetFillStyle(3011)
h_pt_res.DrawNormalized()
legend.Draw ("same")
CMS_lumi.CMS_lumi(c_pt_res, 0, 0)
c_pt_res.Update()
#c_pt_res.SaveAs(outputdir + "/pt_res_log.png")
#c_pt_res.SaveAs(outputdir + "/pt_res_log.pdf")
c_pt_res.SaveAs(outputdir + "/pt_res_zoom.png")
c_pt_res.SaveAs(outputdir + "/pt_res_zoom.pdf")

c_mass_res = ROOT.TCanvas("c_mass_res", "c_mass_res", 1000, 800)
c_mass_res.cd()    
#c_mass_res.SetLogy()    
c_mass_res.SetBottomMargin(0.17)
h_mass_res.GetXaxis().SetLabelOffset(0.02)
h_mass_res.GetXaxis().SetTitleOffset(1.9)
#h_mass_res.GetXaxis().SetTitle("(mass^{scout}-mass^{off})/mass^{off} [GeV]")
h_mass_res.GetXaxis().SetTitle("#frac{m_{#mu#mu}^{scout}-m_{#mu#mu}^{off}}{m_{#mu#mu}^{off}} [GeV]")
h_mass_res.GetYaxis().SetTitleOffset(1.5)
h_mass_res.GetYaxis().SetTitle("A.U.")
#h_mass_res.GetYaxis().SetTitle("0.1 GeV/Event")
h_mass_res.SetLineWidth(2)
h_mass_res.SetLineColor(kAzure-4)
h_mass_res.SetFillColorAlpha(kAzure-9,0.35)
h_mass_res.SetFillStyle(3011)
h_mass_res.DrawNormalized()
legend.Draw ("same")
CMS_lumi.CMS_lumi(c_mass_res, 0, 0)
c_mass_res.Update()
#c_mass_res.SaveAs(outputdir + "/mass_res_log.png")
#c_mass_res.SaveAs(outputdir + "/mass_res_log.pdf")
c_mass_res.SaveAs(outputdir + "/mass_res_zoom.png")
c_mass_res.SaveAs(outputdir + "/mass_res_zoom.pdf")

leg_mass = ROOT.TLegend (0.65, 0.7, 0.85, 0.88)
leg_mass.AddEntry (h_mass_offline, "Offline muons", "L")
leg_mass.AddEntry (h_mass_scout, "Scouting muons", "L")
leg_mass.SetLineWidth (0)

c_mass = ROOT.TCanvas("c_mass", "c_mass", 1000, 800)
c_mass.cd()    
#c_mass.SetLogx()    
c_mass.SetLogy()    
c_mass.SetBottomMargin(0.17)
h_mass_offline.GetXaxis().SetLabelOffset(0.02)
h_mass_offline.GetXaxis().SetTitleOffset(1.9)
h_mass_offline.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
h_mass_offline.GetYaxis().SetTitleOffset(1.5)
h_mass_offline.GetYaxis().SetTitle("A.U.")
#h_mass_offline.GetYaxis().SetTitle("0.5 GeV/Event")
h_mass_offline.SetLineWidth(2)
h_mass_offline.SetLineColor(kOrange-3)
#h_mass_offline.SetFillColorAlpha(kAzure-9,0.35)
h_mass_scout.SetLineWidth(2)
h_mass_scout.SetLineColor(kMagenta-7)
h_mass_scout.SetFillColorAlpha(kMagenta-9,0.35)
h_mass_scout.SetFillStyle(3011)
#h_mass_scout.Scale(h_mass_offline.Integral()/h_mass_scout.Integral())
h_mass_offline.DrawNormalized("hist")
h_mass_scout.DrawNormalized("same hist")
leg_mass.Draw ("same")
CMS_lumi.CMS_lumi(c_mass, 0, 0)
c_mass.Update()
c_mass.SaveAs(outputdir + "/mass_logY.png")
c_mass.SaveAs(outputdir + "/mass_logY.pdf")
#c_mass.SaveAs(outputdir + "/mass_logXY.png")
#c_mass.SaveAs(outputdir + "/mass_logXY.pdf")


###################
# Files and trees #
###################
#f_scoutMuon   = TFile.Open(
#"root://cmsxrootd.fnal.gov//store/user/elfontan/ScoutingPFMonitor/monitorSkim_13Feb2023_2022/230505_183721/0000/scoutMonitor_998.root"
#)

#with open('test_scoutMuon.txt', 'r') as f:
#    # Read the contents of the text file
#    f_scoutMuon_contents = f.read()

# Split the file contents into a list of file names
#f_scoutMuon_names = f_scoutMuon_contents.split(',')

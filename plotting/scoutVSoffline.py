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

lowmass = True
fullmass = False

######################################
# List of files and output directory #
######################################
def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]
files = list_full_paths("/eos/user/e/elfontan/ScoutingParkingPaper/lxy_vtxInfo_Jun26_2022scoutMon")
#files = list_full_paths("/eos/user/e/elfontan/ScoutingParkingPaper/lxy_May25_2022scoutMon") #DEFAULT
#files = files[:-1] 

outputdir = "/eos/user/e/elfontan/www/ScoutingParkingPaper/lxy_May25_2022scoutMon"

########################
# Variables and histos #
########################
h_pt_res_zoom   = TH1F("h_pt_res_zoom", "h_pt_res_zoom", 200, -0.1, 0.1)
h_mass_res_zoom = TH1F("h_mass_res_zoom", "h_mass_res_zoom", 200, -0.1, 0.1)
h_pt_res        = TH1F("h_pt_res", "h_pt_res", 400, -0.4, 0.4) 
h_mass_res      = TH1F("h_mass_res", "h_mass_res", 400, -0.4, 0.4) 

if (fullmass):
    #xbins = [0.215]
    xbins = [0.3]
    while (xbins[-1]<250):
        xbins.append(1.01*xbins[-1])
    #print("xbins", xbins)
if (lowmass):
    xbins = [0.3]
    while (xbins[-1]<20):
        xbins.append(1.01*xbins[-1])
    #print("xbins", xbins)

h_mass_offline = TH1F("h_mass_offline", "h_mass_offline", len(xbins)-1,array('f',xbins)) #LOG
h_mass_scout = TH1F("h_mass_scout", "h_mass_scout", len(xbins)-1,array('f',xbins)) #LOG

if (fullmass):
    frame = TH1F("frame","",1000,0.3,250.3)
elif (lowmass):
    frame = TH1F("frame","",1000,0.3,20.3)

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
      #print("nScoutingMuons = ", ev.nScoutingMuons )
      if (not(ev.nScoutingMuons == 2)): continue
      #if (ev.pt1_scout == ev.pt2_scout): continue
      #if (ev.drmm < 0.2): continue 
      #if (ev.drmm < 0.2 or ev.drmm_scout < 0.2): continue 
      if (ev.dr_matching_1 > 0.1 or ev.dr_matching_2 > 0.1): continue

      if ((ev.l1Result[0]==1 or ev.l1Result[1]==1 or ev.l1Result[2]==1 or ev.l1Result[3]==1 or ev.l1Result[4]==1 or ev.l1Result[5]==1) and ev.ndvtx > 0 ):
      #if ((ev.l1Result[0]==1 or ev.l1Result[1]==1 or ev.l1Result[2]==1 or ev.l1Result[3]==1 or ev.l1Result[4]==1 or ev.l1Result[5]==1) and ev.lxy > 0.0):
          #if (ev.mu1_ID[0] and ev.mu2_ID[0] and ev.pfIso1 < 0.25 and ev.pfIso2 < 0.25):
          #print("pt1 = ", ev.pt1, " and pt1_scout = ", ev.pt1_scout)
          h_pt_res_zoom.Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
          h_pt_res_zoom.Fill((ev.pt2_scout - ev.pt2)/ev.pt2)
          h_pt_res.Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
          h_pt_res.Fill((ev.pt2_scout - ev.pt2)/ev.pt2)
          h_mass_res_zoom.Fill((ev.mass_scout - ev.mass)/ev.mass)
          h_mass_res.Fill((ev.mass_scout - ev.mass)/ev.mass)
          h_mass_offline.Fill(ev.mass)
          h_mass_scout.Fill(ev.mass_scout)
      else:
          continue
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


# --------------------------------------------------------------------
#leg_mass = ROOT.TLegend (0.15, 0.7, 0.45, 0.88)
leg_mass = ROOT.TLegend (0.58, 0.7, 0.88, 0.88)
leg_mass.SetTextSize(0.05)
leg_mass.AddEntry (h_mass_offline, "Offline muons", "L")
leg_mass.AddEntry (h_mass_scout, "Scouting muons", "L")
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
    c_mass = ROOT.TCanvas("c_mass", "c_mass", 1200, 800)
    c_mass.cd()    
    c_mass.SetLogx()    
    c_mass.SetLogy()    
    c_mass.SetLeftMargin(0.11)
    c_mass.SetBottomMargin(0.17)
    frame.SetMinimum(10)
    frame.SetMaximum(10000000)
    frame.GetXaxis().SetRangeUser(0.31,250.)
    #h_mass_offline.GetXaxis().SetRangeUser(0.215, 250)
    #h_mass_scout.GetXaxis().SetRangeUser(0.215,250)
    #h_mass_offline.GetYaxis().SetRangeUser(10, 10000000)
    h_mass_offline.GetXaxis().SetLabelOffset(0.02)
    h_mass_offline.GetXaxis().SetTitleOffset(1.9)
    h_mass_offline.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    h_mass_offline.GetYaxis().SetTitleOffset(1.5)
    h_mass_offline.GetYaxis().SetTitle("A.U.")

    frame.GetXaxis().SetLabelOffset(0.02)
    frame.GetXaxis().SetTitleOffset(1.9)
    frame.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    frame.GetYaxis().SetTitleOffset(1.5)
    frame.GetYaxis().SetTitle("A.U.")
    frame.Draw()

    #h_mass_offline.GetYaxis().SetTitle("0.5 GeV/Event")
    h_mass_offline.SetLineWidth(4)
    h_mass_offline.SetLineStyle(8)
    h_mass_offline.SetLineColor(kOrange-3)
    #h_mass_offline.SetFillColorAlpha(kAzure-9,0.35)
    h_mass_scout.SetLineWidth(2)
    h_mass_scout.SetLineColor(kMagenta-7)
    h_mass_scout.SetFillColorAlpha(kMagenta-9,0.35)
    h_mass_scout.SetFillStyle(3011)
    h_mass_scout.Scale(1, "width")
    h_mass_offline.Scale(1, "width")
    h_mass_scout.Draw("same hist")
    h_mass_offline.Draw("same hist")

    [ labels.DrawLatex(masses[m], 1.5*h_mass_offline.GetBinContent(h_mass_offline.FindBin(masses[m])), m) for m in masses ]
    labels.Draw("same")
    
    leg_mass.Draw ("same")
    CMS_lumi.CMS_lumi(c_mass, 0, 0)
    c_mass.Update()
    c_mass.SaveAs(outputdir + "/mass_logXYwidth_TEST.png")
    c_mass.SaveAs(outputdir + "/mass_logXYwidth_TEST.pdf")
    
if (lowmass):
    c_lowmass = ROOT.TCanvas("c_lowmass", "c_lowmass", 1200, 800)
    c_lowmass.cd()    
    c_lowmass.SetLogx()    
    c_lowmass.SetLogy()    
    c_lowmass.SetLeftMargin(0.11)
    c_lowmass.SetBottomMargin(0.17)
    frame.SetMinimum(100)
    frame.SetMaximum(10000000)
    frame.GetXaxis().SetRangeUser(0.35,20.)
    frame.GetXaxis().SetLabelOffset(0.02)
    frame.GetXaxis().SetTitleOffset(1.9)
    frame.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    frame.GetYaxis().SetTitleOffset(1.5)
    frame.GetYaxis().SetTitle("A.U.")
    frame.Draw()

    #h_mass_scout.GetYaxis().SetRangeUser(10, 1000000000)
    h_mass_offline.SetLineWidth(4)
    h_mass_offline.SetLineStyle(8)
    h_mass_offline.SetLineColor(kOrange-3)
    #h_mass_offline.SetFillColorAlpha(kAzure-9,0.35)
    h_mass_scout.SetLineWidth(2)
    h_mass_scout.SetLineColor(kMagenta-7)
    h_mass_scout.SetFillColorAlpha(kMagenta-9,0.35)
    h_mass_scout.SetFillStyle(3011)
    h_mass_scout.Scale(1., "width")
    h_mass_offline.Scale(1., "width")
    h_mass_scout.Draw("same hist")
    h_mass_offline.Draw("same hist")
    
    [ labels.DrawLatex(masses[m], 1.5*h_mass_offline.GetBinContent(h_mass_scout.FindBin(masses[m])), m) for m in masses ]
    labels.Draw("same")
    
    leg_mass.Draw ("same")
    CMS_lumi.CMS_lumi(c_lowmass, 0, 0)
    c_lowmass.Update()
    c_lowmass.SaveAs(outputdir + "/lowmass_logXYwidth_TEST.png")
    c_lowmass.SaveAs(outputdir + "/lowmass_logXYwidth_TEST.pdf")


# --------------------------------------------------------------------
c_pt_res = ROOT.TCanvas("c_pt_res", "c_pt_res", 1000, 800)
c_pt_res.cd()    
c_pt_res.SetLogy()    
c_pt_res.SetBottomMargin(0.17)
h_pt_res.GetXaxis().SetLabelOffset(0.02)
h_pt_res.GetXaxis().SetTitleOffset(1.9)
h_pt_res.GetXaxis().SetTitle("#frac{p_{T}^{scout}-p_{T}^{off}}{p_{T}^{off}} [GeV]")
h_pt_res.GetYaxis().SetTitleOffset(1.5)
#h_pt_res.GetYaxis().SetTitle("A.U.")
h_pt_res.GetYaxis().SetTitle("0.1 GeV/Event")
h_pt_res.SetLineWidth(2)
h_pt_res.SetLineColor(kAzure-4)
h_pt_res.SetFillColorAlpha(kAzure-9,0.35)
h_pt_res.SetFillStyle(3011)
h_pt_res.Draw()
legend.Draw ("same")
CMS_lumi.CMS_lumi(c_pt_res, 0, 0)
c_pt_res.Update()
#c_pt_res.SaveAs(outputdir + "/pt_res_Matching0p1_minLxy_log.png")
#c_pt_res.SaveAs(outputdir + "/pt_res_Matching0p1_minLxy_log.pdf")

c_mass_res = ROOT.TCanvas("c_mass_res", "c_mass_res", 1000, 800)
c_mass_res.cd()    
c_mass_res.SetLogy()    
c_mass_res.SetBottomMargin(0.17)
h_mass_res.GetXaxis().SetLabelOffset(0.02)
h_mass_res.GetXaxis().SetTitleOffset(1.9)
h_mass_res.GetXaxis().SetTitle("#frac{m_{#mu#mu}^{scout}-m_{#mu#mu}^{off}}{m_{#mu#mu}^{off}} [GeV]")
#h_mass_res.GetYaxis().SetTitleOffset(1.5)
h_mass_res.GetYaxis().SetTitle("0.1 GeV/Event")
h_mass_res.SetLineWidth(2)
h_mass_res.SetLineColor(kAzure-4)
h_mass_res.SetFillColorAlpha(kAzure-9,0.35)
h_mass_res.SetFillStyle(3011)
h_mass_res.Draw()
legend.Draw ("same")
CMS_lumi.CMS_lumi(c_mass_res, 0, 0)
c_mass_res.Update()
#c_mass_res.SaveAs(outputdir + "/mass_res_Matching0p1_minLxy_log.png")
#c_mass_res.SaveAs(outputdir + "/mass_res_Matching0p1_minLxy_log.pdf")

# --------------------------------------------------------------------
c_pt_res_zoom = ROOT.TCanvas("c_pt_res_zoom", "c_pt_res_zoom", 1000, 800)
c_pt_res_zoom.cd()    
c_pt_res_zoom.SetLogy()    
c_pt_res_zoom.SetBottomMargin(0.17)
h_pt_res_zoom.GetXaxis().SetLabelOffset(0.02)
h_pt_res_zoom.GetXaxis().SetTitleOffset(1.9)
h_pt_res_zoom.GetXaxis().SetTitle("#frac{p_{T}^{scout}-p_{T}^{off}}{p_{T}^{off}} [GeV]")
h_pt_res_zoom.GetYaxis().SetTitleOffset(1.5)
#h_pt_res_zoom.GetYaxis().SetTitle("A.U.")
h_pt_res_zoom.GetYaxis().SetTitle("0.1 GeV/Event")
h_pt_res_zoom.SetLineWidth(2)
h_pt_res_zoom.SetLineColor(kAzure-4)
h_pt_res_zoom.SetFillColorAlpha(kAzure-9,0.35)
h_pt_res_zoom.SetFillStyle(3011)
h_pt_res_zoom.Draw()
legend.Draw ("same")
CMS_lumi.CMS_lumi(c_pt_res_zoom, 0, 0)
c_pt_res_zoom.Update()
#c_pt_res_zoom.SaveAs(outputdir + "/pt_res_zoom_Matching0p1_minLxy_log.png")
#c_pt_res_zoom.SaveAs(outputdir + "/pt_res_zoom_Matching0p1_minLxy_log.pdf")

c_mass_res_zoom = ROOT.TCanvas("c_mass_res_zoom", "c_mass_res_zoom", 1000, 800)
c_mass_res_zoom.cd()    
c_mass_res_zoom.SetLogy()    
c_mass_res_zoom.SetBottomMargin(0.17)
h_mass_res_zoom.GetXaxis().SetLabelOffset(0.02)
h_mass_res_zoom.GetXaxis().SetTitleOffset(1.9)
h_mass_res_zoom.GetXaxis().SetTitle("#frac{m_{#mu#mu}^{scout}-m_{#mu#mu}^{off}}{m_{#mu#mu}^{off}} [GeV]")
h_mass_res_zoom.GetYaxis().SetTitleOffset(1.5)
#h_mass_res_zoom.GetYaxis().SetTitle("A.U.")
h_mass_res_zoom.GetYaxis().SetTitle("0.1 GeV/Event")
h_mass_res_zoom.SetLineWidth(2)
h_mass_res_zoom.SetLineColor(kAzure-4)
h_mass_res_zoom.SetFillColorAlpha(kAzure-9,0.35)
h_mass_res_zoom.SetFillStyle(3011)
h_mass_res_zoom.Draw()
legend.Draw ("same")
CMS_lumi.CMS_lumi(c_mass_res_zoom, 0, 0)
c_mass_res_zoom.Update()
#c_mass_res_zoom.SaveAs(outputdir + "/mass_res_all_Matching0p1_minLxy_log.png")
#c_mass_res_zoom.SaveAs(outputdir + "/mass_res_all_Matching0p1_minLxy_log.pdf")


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

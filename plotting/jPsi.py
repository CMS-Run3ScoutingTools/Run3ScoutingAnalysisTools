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

######################################
# List of files and output directory #
######################################
def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]
files = list_full_paths("/eos/user/e/elfontan/ScoutingParkingPaper/bkp_lxy_May25_2022scoutMon")
#files = list_full_paths("/eos/user/e/elfontan/ScoutingParkingPaper/2022_scoutMon_dRmatching")
#files = list_full_paths("/eos/user/e/elfontan/ScoutingParkingPaper/2022_scoutMon_FIX_noID_IsoSaved_noDimuonTriggerSel")
#files = list_full_paths("/eos/user/e/elfontan/ScoutingParkingPaper/2022_scoutMon_noID_IsoSaved_dimuonTriggerSel")
files = files[:-1] 

outputdir = "/eos/user/e/elfontan/www/ScoutingParkingPaper/lxy_May25_2022scoutMon"
#outputdir = "/eos/user/e/elfontan/www/ScoutingParkingPaper/2022_scoutMon_dRmatching"
#outputdir = "/eos/user/e/elfontan/www/ScoutingParkingPaper/2022_scoutMon_FIX_noID_IsoSaved_noDimuonTriggerSel"

########################
# Variables and histos #
########################
#xbins = [0.215]
#while (xbins[-1]<250):
#  xbins.append(1.01*xbins[-1])
#h_mass_offline = TH1F("h_mass_offline", "h_mass_offline", len(xbins)-1,array('f',xbins)) #LOG
#h_mass_scout = TH1F("h_mass_scout", "h_mass_scout", len(xbins)-1,array('f',xbins)) #LOG
h_mass_offline = TH1F("h_mass_offline", "h_mass_offline", 200, 2., 4.) 
h_mass_scout = TH1F("h_mass_scout", "h_mass_scout", 200, 2., 4.)

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
      if (ev.pt1_scout == ev.pt2_scout): continue
      if ((ev.drmm < 0.2 or ev.drmm_scout < 0.2)): continue
      if (ev.dr_matching_1 > 0.1 or ev.dr_matching_2 > 0.1): continue
      if ((ev.l1Result[0]==1 or ev.l1Result[1]==1 or ev.l1Result[2]==1 or ev.l1Result[3]==1 or ev.l1Result[4]==1 or ev.l1Result[5]==1)):
          if (ev.lxy < 0.00001): continue
          if (ev.mass > 2.0 and ev.mass < 4.0):
            #if (ev.mu1_ID[1] and ev.mu2_ID[1] and ev.pfIso1 < 0.25 and ev.pfIso2 < 0.25):
            h_mass_offline.Fill(ev.mass)
            h_mass_scout.Fill(ev.mass_scout)        
    # Close the ROOT file
    root_file.Close()

#draw CMS and lumi text                                                                                                                                           
CMS_lumi.writeExtraText = True
CMS_lumi.lumi_sqrtS      = "2022 (13.6 TeV)"                                                                                                                 
#CMS_lumi.lumi_sqrtS     = lumiText + " (13.6 TeV)"                                                                                                                 
CMS_lumi.cmsTextSize    = 0.6
CMS_lumi.lumiTextSize   = 0.46
CMS_lumi.extraOverCmsTextSize = 0.75
CMS_lumi.relPosX = 0.12

#leg_mass = ROOT.TLegend (0.15, 0.7, 0.45, 0.88)
leg_mass = ROOT.TLegend (0.65, 0.7, 0.85, 0.88)
leg_mass.AddEntry (h_mass_offline, "Offline muons", "L")
leg_mass.AddEntry (h_mass_scout, "Scouting muons", "L")
leg_mass.SetLineWidth (0)

c_jpsi = ROOT.TCanvas("c_jpsi", "c_jpsi", 1000, 800)
c_jpsi.cd()    
#c_jpsi.SetLogx()    
c_jpsi.SetLogy()    
c_jpsi.SetBottomMargin(0.17)
h_mass_offline.GetXaxis().SetRangeUser(2.,4.)
h_mass_offline.GetXaxis().SetLabelOffset(0.02)
h_mass_offline.GetXaxis().SetTitleOffset(1.9)
h_mass_offline.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
h_mass_offline.GetYaxis().SetTitleOffset(1.5)
#h_mass_offline.GetYaxis().SetTitle("A.U.")
h_mass_offline.GetYaxis().SetTitle("0.05 GeV/Event")
h_mass_offline.SetLineWidth(2)
h_mass_offline.SetLineColor(kOrange-3)
#h_mass_offline.SetFillColorAlpha(kAzure-9,0.35)
h_mass_scout.SetLineWidth(2)
h_mass_scout.SetLineColor(kMagenta-7)
h_mass_scout.SetFillColorAlpha(kMagenta-9,0.35)
h_mass_scout.SetFillStyle(3011)
#h_mass_scout.Scale(1., "width")
#h_mass_offline.Scale(1., "width")
h_mass_offline.Draw("EP")
h_mass_scout.Draw("same EP")
leg_mass.Draw ("same")
CMS_lumi.CMS_lumi(c_jpsi, 0, 0)
c_jpsi.Update()
c_jpsi.SaveAs(outputdir + "/jpsi_mass_dR0p1_drmm0p2_lxy_less0p00001.png")
c_jpsi.SaveAs(outputdir + "/jpsi_mass_dR0p1_drmm0p2_lxy_less0p00001.pdf")


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

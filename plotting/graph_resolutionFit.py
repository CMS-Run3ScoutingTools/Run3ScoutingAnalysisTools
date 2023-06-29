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
import numpy as np
import os
import mplhep

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

######################################
# List of files and output directory #
######################################
def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]
files = list_full_paths("/eos/user/e/elfontan/ScoutingParkingPaper/lxy_vtxInfo_Jun26_2022scoutMon")
files = files[:-1] 

outputdir = "/eos/user/e/elfontan/www/ScoutingParkingPaper/resPlots_2022scoutMon/profile"

########################
# Variables and histos #
########################
h_list_B = []
h_name_list_B = []
f_list_B = []
f_name_list_B = []
h_list_E = []
h_name_list_E = []
f_list_E = []
f_name_list_E = []
min_list = []
max_list = []

h_higherPt_list_B = []
h_higherPt_list_E = []
f_higherPt_list_B = []
f_higherPt_list_E = []

for v in range(0,60):
    minval = str(v)
    maxval = str(v+1)
    min_list.append(v)
    max_list.append(v+1)
    h_name_B = "h_res_B_pt_"+minval+"_"+maxval 
    h_name_list_B.append(h_name_B)
    h_res_pt_B =  TH1F(h_name_B, h_name_B, 400, -0.05, 0.05)
    f_name_B = "f_gaus_B_pt_"+minval+"_"+maxval
    f_res_pt_B = ROOT.TF1(f_name_B, "gaus")
    h_list_B.append(h_res_pt_B)
    f_list_B.append(f_res_pt_B)
    f_name_list_B.append(f_name_B)
    h_name_E = "h_res_E_pt_"+minval+"_"+maxval 
    h_name_list_E.append(h_name_E)
    h_res_pt_E =  TH1F(h_name_E, h_name_E, 400, -0.05, 0.05)
    f_name_E = "f_gaus_E_pt_"+minval+"_"+maxval
    f_res_pt_E = ROOT.TF1(f_name_E, "gaus")
    h_list_E.append(h_res_pt_E)
    f_list_E.append(f_res_pt_E)
    f_name_list_E.append(f_name_E)

h_res_pt_60_70_B =  TH1F("h_res_pt_60_70_B", "h_res_pt_60_70_B", 400, -0.05, 0.05)
h_res_pt_70_80_B =  TH1F("h_res_pt_70_80_B", "h_res_pt_70_80_B", 400, -0.05, 0.05)
h_res_pt_80_90_B =  TH1F("h_res_pt_80_90_B", "h_res_pt_80_90_B", 400, -0.05, 0.05)
h_res_pt_90_100_B =  TH1F("h_res_pt_90_100_B", "h_res_pt_90_100_B", 400, -0.05, 0.05)
f_res_pt_60_70_B = ROOT.TF1("f_res_pt_60_70_B", "gaus")
f_res_pt_70_80_B = ROOT.TF1("f_res_pt_70_80_B", "gaus")
f_res_pt_80_90_B = ROOT.TF1("f_res_pt_80_90_B", "gaus")
f_res_pt_90_100_B = ROOT.TF1("f_res_pt_90_100_B", "gaus")

h_res_pt_60_70_E =  TH1F("h_res_pt_60_70_E", "h_res_pt_60_70_E", 400, -0.05, 0.05)
h_res_pt_70_80_E =  TH1F("h_res_pt_70_80_E", "h_res_pt_70_80_E", 400, -0.05, 0.05)
h_res_pt_80_90_E =  TH1F("h_res_pt_80_90_E", "h_res_pt_80_90_E", 400, -0.05, 0.05)
h_res_pt_90_100_E =  TH1F("h_res_pt_90_100_E", "h_res_pt_90_100_E", 400, -0.05, 0.05)
f_res_pt_60_70_E = ROOT.TF1("f_res_pt_60_70_E", "gaus")
f_res_pt_70_80_E = ROOT.TF1("f_res_pt_70_80_E", "gaus")
f_res_pt_80_90_E = ROOT.TF1("f_res_pt_80_90_E", "gaus")
f_res_pt_90_100_E = ROOT.TF1("f_res_pt_90_100_E", "gaus")

h_mass_res_zoom_B = TH1F("h_mass_res_zoom_B", "h_mass_res_zoom_B", 400, -0.05, 0.05)    
h_mass_res_zoom_E = TH1F("h_mass_res_zoom_E", "h_mass_res_zoom_E", 400, -0.05, 0.05)    
f_res_mass_B = ROOT.TF1("f_res_mass_B", "gaus")
f_res_mass_E = ROOT.TF1("f_res_mass_E", "gaus")

x = ROOT.RooRealVar("x", "x", -0.1, 0.1) 
mean = ROOT.RooRealVar("mean", "mean", 0, -0.05, 0.05)
sigma1 = ROOT.RooRealVar("sigma1", "sigma1", 0.01, -0.05, 0.05)
alpha1 = ROOT.RooRealVar("alpha1", "alpha1", 1, 0, 3)
n1 = ROOT.RooRealVar("n1", "n1", 1, 0, 3)
sigma2 = ROOT.RooRealVar("sigma2", "sigma2", 0.01, -0.05, 0.05)
alpha2 = ROOT.RooRealVar("alpha2", "alpha2", -1, -3, 0)
n2 = ROOT.RooRealVar("n2", "n2", 1, 0, 3)
cb1 = ROOT.RooCBShape("cb1", "cb1", x, mean, sigma1, alpha1, n1)
extend_cb1 = ROOT.RooExtendPdf("extend_cb1", "extend_cb1", cb1, n1)
cb2 = ROOT.RooCBShape("cb2", "cb2", x, mean, sigma2, alpha2, n2)
extend_cb2 = ROOT.RooExtendPdf("extend_cb2", "extend_cb2", cb2, n2)
dcbs = ROOT.RooAddPdf("dcbs", "dcbs", ROOT.RooArgList(extend_cb1, extend_cb2))

#print(h_name_list_B)
#print(h_name_list_E)

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
      if (ev.drmm < 0.2 or ev.drmm_scout < 0.2): continue                                                                                                            
      if (ev.dr_matching_1 > 0.1 or ev.dr_matching_2 > 0.1): continue

      if ((ev.l1Result[0]==1 or ev.l1Result[1]==1 or ev.l1Result[2]==1 or ev.l1Result[3]==1 or ev.l1Result[4]==1 or ev.l1Result[5]==1) and ndvtx > 0):
        if (abs(ev.eta1) > 0.9 and abs(ev.eta1) < 1.9 and abs(ev.eta2) > 0.9 and abs(ev.eta2) < 1.9):
            h_mass_res_zoom_E.Fill((ev.mass_scout - ev.mass)/ev.mass)
        elif(abs(ev.eta1) < 0.9 and abs(ev.eta2) < 0.9):
            h_mass_res_zoom_B.Fill((ev.mass_scout - ev.mass)/ev.mass)

        for h in range(len(h_list_B)):
            #print(h)
            if (ev.pt1_scout > min_list[h] and ev.pt1_scout <= max_list[h] and abs(ev.eta1) < 0.9):
                #print("min = ", min_list[h])
                #print("max = ", max_list[h])
                h_list_B[h].Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
            if (ev.pt2_scout > min_list[h] and ev.pt2_scout <= max_list[h] and abs(ev.eta2) < 0.9):
                h_list_B[h].Fill((ev.pt2_scout - ev.pt2)/ev.pt2)
                    
            if (ev.pt1_scout > min_list[h] and ev.pt1_scout <= max_list[h] and abs(ev.eta1) > 0.9 and abs(ev.eta1) < 1.9):
                h_list_E[h].Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
            if (ev.pt2_scout > min_list[h] and ev.pt2_scout <= max_list[h] and abs(ev.eta2) > 0.9 and abs(ev.eta2) < 1.9):
                h_list_E[h].Fill((ev.pt2_scout - ev.pt2)/ev.pt2)

        if (ev.pt1_scout > 60 and ev.pt1_scout <= 70 and abs(ev.eta1) < 0.9):
            h_res_pt_60_70_B.Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
        if (ev.pt2_scout > 60 and ev.pt2_scout <= 70 and abs(ev.eta2) < 0.9):
            h_res_pt_60_70_B.Fill((ev.pt2_scout - ev.pt2)/ev.pt2)
        if (ev.pt1_scout > 60 and ev.pt1_scout <= 70 and abs(ev.eta1) > 0.9 and abs(ev.eta1) < 1.9):
            h_res_pt_60_70_E.Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
        if (ev.pt2_scout > 60 and ev.pt2_scout <= 70 and abs(ev.eta2) > 0.9 and abs(ev.eta2) < 1.9):
            h_res_pt_60_70_E.Fill((ev.pt2_scout - ev.pt2)/ev.pt2)

        if (ev.pt1_scout > 70 and ev.pt1_scout <= 80 and abs(ev.eta1) < 0.9):
            h_res_pt_70_80_B.Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
        if (ev.pt2_scout > 70 and ev.pt2_scout <= 80 and abs(ev.eta2) < 0.9):
            h_res_pt_70_80_B.Fill((ev.pt2_scout - ev.pt2)/ev.pt2)
        if (ev.pt1_scout > 70 and ev.pt1_scout <= 80 and abs(ev.eta1) > 0.9 and abs(ev.eta1) < 1.9):
            h_res_pt_70_80_E.Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
        if (ev.pt2_scout > 70 and ev.pt2_scout <= 80 and abs(ev.eta2) > 0.9 and abs(ev.eta2) < 1.9):
            h_res_pt_70_80_E.Fill((ev.pt2_scout - ev.pt2)/ev.pt2)

        if (ev.pt1_scout > 80 and ev.pt1_scout <= 90 and abs(ev.eta1) < 0.9):
            h_res_pt_80_90_B.Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
        if (ev.pt2_scout > 80 and ev.pt2_scout <= 90 and abs(ev.eta2) < 0.9):
            h_res_pt_80_90_B.Fill((ev.pt2_scout - ev.pt2)/ev.pt2)
        if (ev.pt1_scout > 80 and ev.pt1_scout <= 90 and abs(ev.eta1) > 0.9 and abs(ev.eta1) < 1.9):
            h_res_pt_80_90_E.Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
        if (ev.pt2_scout > 80 and ev.pt2_scout <= 90 and abs(ev.eta2) > 0.9 and abs(ev.eta2) < 1.9):
            h_res_pt_80_90_E.Fill((ev.pt2_scout - ev.pt2)/ev.pt2)

        if (ev.pt1_scout > 90 and ev.pt1_scout <= 100 and abs(ev.eta1) < 0.9):
            h_res_pt_90_100_B.Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
        if (ev.pt2_scout > 90 and ev.pt2_scout <= 100 and abs(ev.eta2) < 0.9):
            h_res_pt_90_100_B.Fill((ev.pt2_scout - ev.pt2)/ev.pt2)
        if (ev.pt1_scout > 90 and ev.pt1_scout <= 100 and abs(ev.eta1) > 0.9 and abs(ev.eta1) < 1.9):
            h_res_pt_90_100_E.Fill((ev.pt1_scout - ev.pt1)/ev.pt1)
        if (ev.pt2_scout > 90 and ev.pt2_scout <= 100 and abs(ev.eta2) > 0.9 and abs(ev.eta2) < 1.9):
            h_res_pt_90_100_E.Fill((ev.pt2_scout - ev.pt2)/ev.pt2)

    h_higherPt_list_B.append(h_res_pt_60_70_B)
    h_higherPt_list_B.append(h_res_pt_70_80_B)
    h_higherPt_list_B.append(h_res_pt_80_90_B)
    h_higherPt_list_B.append(h_res_pt_90_100_B)
    h_higherPt_list_E.append(h_res_pt_60_70_E)
    h_higherPt_list_E.append(h_res_pt_70_80_E)
    h_higherPt_list_E.append(h_res_pt_80_90_E)
    h_higherPt_list_E.append(h_res_pt_90_100_E)

    # Close the ROOT file
    root_file.Close()


# ----------------------------------------------------------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------------------------------------------------------

# Legend and title
# ----------------
legend_B = ROOT.TLegend (0.56, 0.6, 0.86, 0.86)
legend_B.SetTextSize (0.032)
legend_B.AddEntry (h_mass_res_zoom_B, "Uncorrected muons (B)", "L")
legend_B.SetLineWidth (0)
legend_E = ROOT.TLegend (0.56, 0.6, 0.86, 0.86)
legend_E.SetTextSize (0.032)
legend_E.AddEntry (h_mass_res_zoom_E, "Uncorrected muons (E)", "L")
legend_E.SetLineWidth (0)

CMS_lumi.writeExtraText = True
CMS_lumi.lumi_sqrtS      = "2022 (13.6 TeV)"                                                                                                                 
CMS_lumi.cmsTextSize    = 0.6
CMS_lumi.lumiTextSize   = 0.46
CMS_lumi.extraOverCmsTextSize = 0.75
CMS_lumi.relPosX = 0.12

rms_list_B = []
rms_list_E = []
rmserr_list_B = []
rmserr_list_E = []

# PT resolution for 1 GeV bins in the barrel region
# -------------------------------------------------
for h in range(len(h_list_B)):
    #print(h_list_B[h])
    #print(h_list_E[h])
    c_name_B = "c_res_pt_B_"+str(min_list[h])+"_"+str(max_list[h])
    c_res_pt_B = ROOT.TCanvas(c_name_B, c_name_B, 1000, 800)
    c_res_pt_B.cd()
    #c_res_pt_B.SetLogy()

    c_res_pt_B.SetBottomMargin(0.07)
    h_list_B[h].GetXaxis().SetLabelOffset(0.02)
    h_list_B[h].GetXaxis().SetTitleOffset(1.9)
    h_list_B[h].GetXaxis().SetTitle("#frac{p_{T}^{scout}-p_{T}^{off}}{p_{T}^{off}}")
    h_list_B[h].GetYaxis().SetTitleOffset(1.5)
    if (min_list[h] <= 10):
        h_list_B[h].GetXaxis().SetRangeUser(-0.05, 0.05)
        h_list_B[h].GetYaxis().SetTitle("0.25 MeV/Event")
    else:
        h_list_B[h].GetYaxis().SetTitle("0.5 MeV/Event")
    #h_list_B[h].GetYaxis().SetTitle("A.U.")
    h_list_B[h].SetLineWidth(2)
    h_list_B[h].SetLineColor(kAzure-4)
    h_list_B[h].SetFillColorAlpha(kAzure-9,0.35)
    h_list_B[h].SetFillStyle(3011)
    h_list_B[h].Draw()
    
    f_list_B[h].SetLineWidth(2)
    f_list_B[h].SetLineColor(kBlue)
    if (min_list[h] <= 10):
        fit_result_B = h_list_B[h].Fit(f_list_B[h], "", "", -0.05, 0.05)
    else: 
        fit_result_B = h_list_B[h].Fit(f_list_B[h], "", "", -0.1, 0.1)

    print("########################")
    fit_pt_results_B = f_list_B[h].GetParameters()
    fit_pt_errors_B = f_list_B[h].GetParErrors()
    print("Fit results Pt - B:")
    print("Mean:", fit_pt_results_B[1], "+/-", fit_pt_errors_B[1])
    print("Standard Deviation:", fit_pt_results_B[2], "+/-", fit_pt_errors_B[2])
    # Calculate the RMS
    rms_value_B = fit_pt_results_B[2] * ROOT.TMath.Sqrt(2)
    rms_list_B.append(rms_value_B)    
    rmserr_list_B.append(fit_pt_errors_B[2]*ROOT.TMath.Sqrt(2))
    print("-----------------")    

    text_box_B = ROOT.TPaveText(0.68, 0.6, 0.83, 0.7, "NDC")
    text_box_B.AddText("RMS = %.3f" % rms_value_B)
    text_box_B.SetTextFont(42)
    text_box_B.SetTextSize(0.03)
    text_box_B.SetFillColor(0)

    #fit_pt0_10_results = f_list[h].GetParameters()
    #fit_pt0_10_errors = f_list[h].GetParErrors()
    f_list_B[h].Draw("same")
    
    legend_B.Draw ("same")
    text_box_B.Draw("same")
    CMS_lumi.CMS_lumi(c_res_pt_B, 0, 0)
    c_res_pt_B.Update()
    c_res_pt_B.SaveAs(outputdir + "/pt"+str(min_list[h])+"_"+str(max_list[h])+"_res_zoom_fit_B.png")
    c_res_pt_B.SaveAs(outputdir + "/pt"+str(min_list[h])+"_"+str(max_list[h])+"_res_zoom_fit_B.pdf")
    
# PT resolution for 1 GeV bins in the endcap region
# -------------------------------------------------
for h in range(len(h_list_E)):
    c_name_E = "c_res_pt_E_"+str(min_list[h])+"_"+str(max_list[h])
    c_res_pt_E = ROOT.TCanvas(c_name_E, c_name_E, 1000, 800)
    c_res_pt_E.cd()
    #c_res_pt_E.SetLogy()

    c_res_pt_E.SetBottomMargin(0.17)

    h_list_E[h].GetXaxis().SetLabelOffset(0.02)
    h_list_E[h].GetXaxis().SetTitleOffset(1.9)
    h_list_E[h].GetXaxis().SetTitle("#frac{p_{T}^{scout}-p_{T}^{off}}{p_{T}^{off}}")
    h_list_E[h].GetYaxis().SetTitleOffset(1.5)
    if (min_list[h] <= 10):
        h_list_E[h].GetXaxis().SetRangeUser(-0.05, 0.05)
        h_list_E[h].GetYaxis().SetTitle("0.25 MeV/Event")
    else:
        h_list_E[h].GetYaxis().SetTitle("0.5 MeV/Event")
    #h_list_E[h].GetYaxis().SetTitle("A.U.")
    h_list_E[h].SetLineWidth(2)
    h_list_E[h].SetLineColor(kAzure-4)
    h_list_E[h].SetFillColorAlpha(kAzure-9,0.35)
    h_list_E[h].SetFillStyle(3011)
    h_list_E[h].Draw()
    
    f_list_E[h].SetLineWidth(2)
    f_list_E[h].SetLineColor(kOrange)
    if (min_list[h] <= 10):
        fit_result_E = h_list_E[h].Fit(f_list_E[h], "", "", -0.05, 0.05)
    else: 
        fit_result_E = h_list_E[h].Fit(f_list_E[h], "", "", -0.1, 0.1)

    print("########################")
    fit_pt_results_E = f_list_E[h].GetParameters()
    fit_pt_errors_E = f_list_E[h].GetParErrors()
    print("Fit results Pt - B:")
    print("Mean:", fit_pt_results_E[1], "+/-", fit_pt_errors_E[1])
    print("Standard Deviation:", fit_pt_results_E[2], "+/-", fit_pt_errors_E[2])
    # Calculate the RMS
    rms_value_E = fit_pt_results_E[2] * ROOT.TMath.Sqrt(2)
    rms_list_E.append(rms_value_E)
    rmserr_list_E.append(fit_pt_errors_E[2]*ROOT.TMath.Sqrt(2))
    print("-----------------")    

    text_box_E = ROOT.TPaveText(0.68, 0.6, 0.83, 0.7, "NDC")
    text_box_E.AddText("RMS = %.3f" % rms_value_E)
    text_box_E.SetTextFont(42)
    text_box_E.SetTextSize(0.03)
    text_box_E.SetFillColor(0)

    f_list_E[h].Draw("same")
    
    legend_E.Draw ("same")
    text_box_E.Draw("same")
    CMS_lumi.CMS_lumi(c_res_pt_E, 0, 0)
    c_res_pt_E.Update()
    c_res_pt_E.SaveAs(outputdir + "/pt"+str(min_list[h])+"_"+str(max_list[h])+"_res_zoom_fit_E.png")
    c_res_pt_E.SaveAs(outputdir + "/pt"+str(min_list[h])+"_"+str(max_list[h])+"_res_zoom_fit_E.pdf")


# Higher Pt: Barrel and endcap histos
# -----------------------------------
higherPt_min_list = [60, 70, 80, 90]
higherPt_max_list = [70, 80, 90, 100]
rms_higherPt_list_B = []
rms_higherPt_list_E = []
rmserr_higherPt_list_B = []
rmserr_higherPt_list_E = []

for h in range(len(higherPt_min_list)):
    f_higherPt_name_B = "f_gaus_B_pt_"+str(higherPt_min_list[h])+"_"+str(higherPt_max_list[h])
    f_res_pt_B = ROOT.TF1(f_higherPt_name_B, "gaus")
    f_higherPt_list_B.append(f_res_pt_B)
    f_higherPt_name_E = "f_gaus_E_pt_"+str(higherPt_min_list[h])+"_"+str(higherPt_max_list[h])
    f_res_pt_E = ROOT.TF1(f_higherPt_name_E, "gaus")
    f_higherPt_list_E.append(f_res_pt_E)

    c_name_B = "c_res_pt_B_"+str(higherPt_min_list[h])+"_"+str(higherPt_max_list[h])
    c_res_pt_B = ROOT.TCanvas(c_name_B, c_name_B, 1000, 800)
    c_res_pt_B.cd()
    c_res_pt_B.SetBottomMargin(0.17)
    h_higherPt_list_B[h].GetXaxis().SetLabelOffset(0.02)
    h_higherPt_list_B[h].GetXaxis().SetTitleOffset(1.9)
    h_higherPt_list_B[h].GetXaxis().SetTitle("#frac{p_{T}^{scout}-p_{T}^{off}}{p_{T}^{off}}")
    h_higherPt_list_B[h].GetYaxis().SetTitleOffset(1.5)
    h_higherPt_list_B[h].GetYaxis().SetTitle("0.5 MeV/Event")
    h_higherPt_list_B[h].SetLineWidth(2)
    h_higherPt_list_B[h].SetLineColor(kAzure-4)
    h_higherPt_list_B[h].SetFillColorAlpha(kAzure-9,0.35)
    h_higherPt_list_B[h].SetFillStyle(3011)
    h_higherPt_list_B[h].Draw()
    
    f_higherPt_list_B[h].SetLineWidth(2)
    f_higherPt_list_B[h].SetLineColor(kOrange)
    fit_result_B = h_higherPt_list_B[h].Fit(f_higherPt_list_B[h], "", "", -0.1, 0.1)
    fit_pt_results_B = f_higherPt_list_B[h].GetParameters()
    fit_pt_errors_B = f_higherPt_list_B[h].GetParErrors()
    rms_value_B = fit_pt_results_B[2] * ROOT.TMath.Sqrt(2)
    rms_higherPt_list_B.append(rms_value_B)
    rmserr_higherPt_list_B.append(fit_pt_errors_B[2]*ROOT.TMath.Sqrt(2))

    text_box_B = ROOT.TPaveText(0.68, 0.6, 0.83, 0.7, "NDC")
    text_box_B.AddText("RMS = %.3f" % rms_value_B)
    text_box_B.SetTextFont(42)
    text_box_B.SetTextSize(0.03)
    text_box_B.SetFillColor(0)
    
    f_higherPt_list_B[h].Draw("same")
    
    legend_B.Draw ("same")
    text_box_B.Draw("same")
    CMS_lumi.CMS_lumi(c_res_pt_B, 0, 0)
    c_res_pt_B.Update()
    c_res_pt_B.SaveAs(outputdir + "/pt"+str(higherPt_min_list[h])+"_"+str(higherPt_max_list[h])+"_res_zoom_fit_B.png")
    c_res_pt_B.SaveAs(outputdir + "/pt"+str(higherPt_min_list[h])+"_"+str(higherPt_max_list[h])+"_res_zoom_fit_B.pdf")


    c_name_E = "c_res_pt_E_"+str(higherPt_min_list[h])+"_"+str(higherPt_max_list[h])
    c_res_pt_E = ROOT.TCanvas(c_name_E, c_name_E, 1000, 800)
    c_res_pt_E.cd()
    c_res_pt_E.SetBottomMargin(0.17)
    h_higherPt_list_E[h].GetXaxis().SetLabelOffset(0.02)
    h_higherPt_list_E[h].GetXaxis().SetTitleOffset(1.9)
    h_higherPt_list_E[h].GetXaxis().SetTitle("#frac{p_{T}^{scout}-p_{T}^{off}}{p_{T}^{off}}")
    h_higherPt_list_E[h].GetYaxis().SetTitleOffset(1.5)
    h_higherPt_list_E[h].GetYaxis().SetTitle("0.5 MeV/Event")
    h_higherPt_list_E[h].SetLineWidth(2)
    h_higherPt_list_E[h].SetLineColor(kAzure-4)
    h_higherPt_list_E[h].SetFillColorAlpha(kAzure-9,0.35)
    h_higherPt_list_E[h].SetFillStyle(3011)
    h_higherPt_list_E[h].Draw()
    
    f_higherPt_list_E[h].SetLineWidth(2)
    f_higherPt_list_E[h].SetLineColor(kOrange)
    fit_result_E = h_higherPt_list_E[h].Fit(f_higherPt_list_E[h], "", "", -0.1, 0.1)
    fit_pt_results_E = f_higherPt_list_E[h].GetParameters()
    fit_pt_errors_E = f_higherPt_list_E[h].GetParErrors()
    rms_value_E = fit_pt_results_E[2] * ROOT.TMath.Sqrt(2)
    rms_higherPt_list_E.append(rms_value_E)
    rmserr_higherPt_list_E.append(fit_pt_errors_E[2]*ROOT.TMath.Sqrt(2))

    text_box_E = ROOT.TPaveText(0.68, 0.6, 0.83, 0.7, "NDC")
    text_box_E.AddText("RMS = %.3f" % rms_value_E)
    text_box_E.SetTextFont(42)
    text_box_E.SetTextSize(0.03)
    text_box_E.SetFillColor(0)
    
    f_higherPt_list_E[h].Draw("same")
    
    legend_E.Draw ("same")
    text_box_E.Draw("same")
    CMS_lumi.CMS_lumi(c_res_pt_E, 0, 0)
    c_res_pt_E.Update()
    c_res_pt_E.SaveAs(outputdir + "/pt"+str(higherPt_min_list[h])+"_"+str(higherPt_max_list[h])+"_res_zoom_fit_E.png")
    c_res_pt_E.SaveAs(outputdir + "/pt"+str(higherPt_min_list[h])+"_"+str(higherPt_max_list[h])+"_res_zoom_fit_E.pdf")
    

# Create a TGraphError for resolution
# ----------------------------------------------
x_values = []
for i in range(0,100):
    x_values.append(i)

h_ptres_profile_B = ROOT.TGraphErrors(len(x_values))
h_ptres_profile_E = ROOT.TGraphErrors(len(x_values))
for bin in range(3,60):
    h_ptres_profile_B.SetPoint(bin, x_values[bin], rms_list_B[bin])
    h_ptres_profile_B.SetPointError(bin, 1, rmserr_list_B[bin])
    h_ptres_profile_E.SetPoint(bin, x_values[bin], rms_list_E[bin])
    h_ptres_profile_E.SetPointError(bin, 1, rmserr_list_E[bin])
h_ptres_profile_B.SetPoint(65, 65, rms_higherPt_list_B[0])
h_ptres_profile_B.SetPointError(65, 5, rmserr_higherPt_list_B[0])
h_ptres_profile_B.SetPoint(75, 75, rms_higherPt_list_B[1])
h_ptres_profile_B.SetPointError(75, 5, rmserr_higherPt_list_B[1])
h_ptres_profile_B.SetPoint(85, 85, rms_higherPt_list_B[2])
h_ptres_profile_B.SetPointError(85, 5, rmserr_higherPt_list_B[2])
h_ptres_profile_B.SetPoint(95, 95, rms_higherPt_list_B[3])
h_ptres_profile_B.SetPointError(95, 5, rmserr_higherPt_list_B[3])

h_ptres_profile_E.SetPoint(65, 65, rms_higherPt_list_E[0])
h_ptres_profile_E.SetPointError(65, 5, rmserr_higherPt_list_E[0])
h_ptres_profile_E.SetPoint(75, 75, rms_higherPt_list_E[1])
h_ptres_profile_E.SetPointError(75, 5, rmserr_higherPt_list_E[1])
h_ptres_profile_E.SetPoint(85, 85, rms_higherPt_list_E[2])
h_ptres_profile_E.SetPointError(85, 5, rmserr_higherPt_list_E[2])
h_ptres_profile_E.SetPoint(95, 95, rms_higherPt_list_E[3])
h_ptres_profile_E.SetPointError(95, 5, rmserr_higherPt_list_E[3])

# Create a canvas and draw the histogram
c_graph = ROOT.TCanvas("c_graph", "c_graph", 1000, 800)
c_graph.SetLeftMargin(0.17)
c_graph.SetBottomMargin(0.17)
h_ptres_profile_B.GetXaxis().SetRangeUser(2,100)
h_ptres_profile_B.GetYaxis().SetRangeUser(0,0.026)
h_ptres_profile_E.GetYaxis().SetRangeUser(0,0.026)
h_ptres_profile_E.GetXaxis().SetRangeUser(2,100)
h_ptres_profile_B.GetXaxis().SetTitle("p_{T}^{#mu} [GeV]")
h_ptres_profile_B.GetYaxis().SetTitle("RMS ( #frac{p_{T}^{scout}-p_{T}^{off}}{p_{T}^{off}} )")
h_ptres_profile_B.GetYaxis().SetTitleOffset(1.9)
h_ptres_profile_B.GetXaxis().SetLabelOffset(0.02)
h_ptres_profile_B.GetXaxis().SetTitleOffset(1.9)
    
# Set marker style and color for value1
h_ptres_profile_B.SetMarkerStyle(20)
h_ptres_profile_B.SetLineColor(kAzure-3)
h_ptres_profile_B.SetMarkerColor(kAzure-3)

# Set marker style and color for value2
h_ptres_profile_E.SetMarkerStyle(23)
h_ptres_profile_E.SetLineColor(kBlue-9)
h_ptres_profile_E.SetMarkerColor(kBlue-9)

h_ptres_profile_B.Draw("AP")  # "P" option for marker plotting
h_ptres_profile_E.Draw("P same")  # "P" option for marker plotting

# Customize the legend
gr_text = ROOT.TPaveText(0.18, 0.78, 0.65, 0.83, "NDC")
gr_text.AddText("Dimuon events, dR_{#mu#mu} > 0.2, p_{T}^{#mu} > 3 GeV")
gr_text.SetTextFont(62)
gr_text.SetTextSize(0.032)
gr_text.SetFillColor(0)


leg_profile = ROOT.TLegend(0.2, 0.68, 0.5, 0.78)
leg_profile.SetTextSize (0.032)
leg_profile.SetLineWidth(0)
leg_profile.AddEntry(h_ptres_profile_B, "Barrel", "P")
leg_profile.AddEntry(h_ptres_profile_E, "Endcap", "P")
leg_profile.Draw()
gr_text.Draw("same")

c_graph.Draw()
CMS_lumi.CMS_lumi(c_graph, 0, 0)
c_graph.Update()
    
c_graph.SaveAs(outputdir + "/ptres_graph_BE.png")
c_graph.SaveAs(outputdir + "/ptres_graph_BE.pdf")

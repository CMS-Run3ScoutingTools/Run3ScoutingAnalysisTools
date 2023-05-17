import os, sys, imp
import ROOT
from array import array
import CMS_lumi

def makeCMSCanvas(name="canvas",title="canvas",width=800,height=600):
    canvas = ROOT.TCanvas(name,name,50,50,width,height)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    T = 0.08*600
    B = 0.12*600
    L = 0.12*800
    R = 0.04*800
    canvas.SetLeftMargin( L/width )
    canvas.SetRightMargin( R/width )
    canvas.SetTopMargin( T/height )
    canvas.SetBottomMargin( B/height )
    canvas.SetTickx(0)
    canvas.SetTicky(0)
    if width == 800: ROOT.gStyle.SetTitleYOffset(1)
    return canvas

def makeLegend(nentries=1):
    height = nentries*0.06
    leg = ROOT.TLegend(0.56,(0.88 - height),0.92,0.88)
    leg.SetMargin(0.)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)
    leg.SetTextFont(20)
    return leg

def printLumiPrelLeft(canvas, lumitext="13 TeV"):
    #change the CMS_lumi variables (see CMS_lumi.py)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_sqrtS = lumitext # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    
    iPos = 11 # inside frame left, default
    
    if ( iPos==0 ): CMS_lumi.relPosX = 0.12
    iPeriod = 4
    
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)


def printLumiPrelOut(canvas, lumitext="13 TeV"):
    #change the CMS_lumi variables (see CMS_lumi.py)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_sqrtS = lumitext # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    
    iPos = 0 # outside frame left
    
    if ( iPos==0 ): CMS_lumi.relPosX = 0.12
    iPeriod = 4
    
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)

def printLumiLeft(canvas, lumitext="13 TeV"):
    #change the CMS_lumi variables (see CMS_lumi.py)
    CMS_lumi.writeExtraText = 0
    CMS_lumi.extraText = ""
    CMS_lumi.lumi_sqrtS = lumitext # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    
    iPos = 11 # inside frame left, default
    
    if ( iPos==0 ): CMS_lumi.relPosX = 0.12
    iPeriod = 4
    
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)


def printLumiOut(canvas, lumitext="13 TeV"):
    #change the CMS_lumi variables (see CMS_lumi.py)
    CMS_lumi.writeExtraText = 0
    CMS_lumi.extraText = ""
    CMS_lumi.lumi_sqrtS = lumitext # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    
    iPos = 0 # outside frame left
    
    if ( iPos==0 ): CMS_lumi.relPosX = 0.12
    iPeriod = 4
    
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)

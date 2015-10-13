import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ()

options.register ('dataset',
				  'none',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "Dataset to process")
options.parseArguments()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )


print 'dataset=',options.dataset

process.load("Analysis.Presel.Presel_cfi")

process.presel.hltPaths = cms.vstring ( "HLT_PFHT900_v.*")

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    #fileNames = cms.untracked.vstring(

     #   'file:/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/eschmitz/b2g/Spring15/TToTH/TprimeTToTH_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/TprimeTToTH_M-1200_LH/150620_211013/0000/B2GEDMNtuple_1.root',
 #   )
    fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/group/phys_b2g/B2GAnaFW/TprimeBToTH_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_v74x_v6p1_25ns/150930_172819/0000/B2GEDMNtuple_1.root')


)

process.TFileService = cms.Service("TFileService", fileName = cms.string(options.dataset+".root") )


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("PreselEvents_"+options.dataset+".root"),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('p')
      ),
    outputCommands = cms.untracked.vstring(
      "drop *",
      "keep *_eventInfo_*_*",
      "keep *_met_*_*",
      "keep *_jetsAK4_*_*",
      "keep *_jetsAK8_*_*",
      "keep *_AK8Jets_*_*",
      "keep *_subjetsAK8_*_*",
      "keep *_subjetsCmsTopTag_*_*",
      "keep *_hbb_*_*",
      "keep *_presel_*_*",
      "keep *_genPart_*_*",
      "keep *_TriggerUserData_*_*"
      )
    )

#process.p = cms.Path(process.hbb*process.presel)
process.p = cms.Path(process.presel)

process.outpath = cms.EndPath(process.out)

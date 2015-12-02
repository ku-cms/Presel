from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'TPrime1200_LH'
config.section_('JobType')
config.JobType.psetName = 'ConfFileNew_cfg.py'
config.JobType.pyCfgParams = ['dataset=TPrime1200_LH']
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['PreselEvents_TPrime1200_LH.root','TPrime1200_LH.root']
config.JobType.disableAutomaticOutputCollection = True
config.section_('Data')
config.Data.inputDataset = '/TprimeBToTH_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/phys_b2g-B2GAnaFW_v74x_v6p1_25ns-d5e3976e154b1b4043d031aaa9c2809b/USER'
config.Data.publication = False
config.Data.unitsPerJob = 20 
#config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.splitting = 'FileBased'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/'
config.Data.outLFNDirBase = '/store/user/stringer/TPrime'
#config.Data.publishDataName = 'my_publication_name'
config.section_('Site')
#config.Site.blacklist = ['T2_IT_Legnaro']
#config.Site.whitelist = ['T2_IT_Bari']
config.Site.storageSite = 'T2_US_Nebraska'

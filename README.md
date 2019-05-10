# ProtoTyping UPS product
Collection of LarSoft analysers to test out things. This is a stand-alone ups-like product, that is tested and working in Uboonecode v07_11_00 and higher.

## NuCC
Charged current neutrino selection using advanced slice ID neutrino ID and pandora consolidated reconstruction.

## CosmicStudies
This analyser stores PMT and TPC information and is used to run on CORSICA and BNBext files.
It includes reco-truth matching for both slimmed and non-slimmed files.
Momentum of tracks is stored using MCS. 
Flashes are stored for all combnations of op/simple reconstruction and cosmic/beam discriminators.
Truth level information is stored of all MC particles above 100MeV, crossing or not crossing the TPC.

### CosmicStudies: CRTstudies
Module to store the CRThit info in data and MC

## CosmicStudies:job
Contains filelists of different samples and fcl file to run over them. The fcl file combine the two analysers (CosmicStudies and CRTstudies) in one project. A set of xml files to run over MC and BNB extrenal are provided.



### Contact
For further questions:
woutervdp@g.harvard.edu

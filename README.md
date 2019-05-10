# ProtoTyping UPS product
Collection of LarSoft analysers to test out things. This is a stand-alone ups-like product, that is tested and working in Uboonecode v08_00_00_12 and higher.

## NuCC
Charged current neutrino selection using advanced slice neutrino ID and pandora consolidated reconstruction.

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


### Setup instructions

```
setup uboonecode v08_00_00_12 -q e17:prof
mrb newDev
mrbsetenv
```
From now on, always do:
```
source /uboone/app/users/xxxx/Binaries/Uboonecode/ubc_v08_00_00_12/localProducts_larsoft_v08_05_00_04a_e17_prof/setup
setup ninja v1_8_2
```
When initiating a new interactive session.

Go to `srcs` directory:
```
mrb g larpandora
cd larpandora
git checkout origin/wvdp_larpandora_overlay_v12

cd ..
mrb g ubcrt
cd ubcrt
git checkout origin/mstancar_crtreco_tag_v08_00_00_12

cd..
mrb g ubreco
cd ubreco
git checkout origin/wvdp_improved_cuts_and_chi

cd..
git clone https://github.com/Wouter-VDP/ProtoTyping.git
mrb uc

```
Try to compile for the first time, go to the build directory and do:
```
mrb uc
mrb i --generator=ninja
mrbslp
mrbslp
```
For compilation after this first time, do inside the build directory:
```
ninja install
```

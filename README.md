Program for Monte Carlo simulation of polymer nanofiber crystallization.
In the foder code-run, 
[1] Create INPUT.DAT file to set up the simulation system 
(Example: bulk system for 100 chains of C60H122; 30 monomer beads): 
--------------------------------------------------------------------------
Input file for MC simulation of polymer on the 2nnd lattice
POLYETHYLENE       ; polymer              (A18)
JACS 1966,88,631   ; literature           (A47)
      30           ; chainlength          (I8)
      24           ; latticesize_X        (I8)
      24           ; latticesize_Y        (I8)
      24           ; latticesize_Z        (I8)
     100           ; chainnumber          (I8)
 1000000           ; chainmoves, MCS      (I8)
   10000           ; storeinterval        (I8)
       0           ; isequence (0=random) (I8)
   90321           ; initseed             (I8)
     473           ; temperature (K)      (I8)
       1           ; ifexcludedvolume     (I8)
       0           ; ifring               (I8)
       1           ; nos                  (I8)
       0           ; iformfilm            (I8)
  15.048           ; energyinterfirstshell(F8.3)
   0.620           ; energyintersecondshel(F8.3)
  -0.625           ; energyinterthirdshell(F8.3)
  --------------------------------------------------------------------------
[2] Run simulation 1 million MCS to get CONFORMATIONS file
[3] Then cut the final snapshot to create CONFORMATION.START file to generate nanofilm (set iformfilm = 1), run another 1 million MCS
     1           ; iformfilm        (I8)
[4] Then cut the final snapshot to create CONFORMATION.START file to generate nanofiber (set iformfilm = 2), run another 1 million MCS
     2           ; iformfilm        (I8)
[5] After than the production for 100 million MCS at each tempertures n the range of 300-370K
Change the last three lines for LJ Intermolecular interaction 
energyinterfirstshell(F8.3)  energyintersecondshel(F8.3) energyinterthirdshell(F8.3)
370 11.907 0.288 -0.633
360 11.594 0.254 -0.634
350 11.281 0.220 -0.635
340 10.966 0.186 -0.636
330 10.649 0.151 -0.637
320 10.331 0.116 -0.639
310 10.011 0.081 -0.640
300  9.689 0.045 -0.641
=============================================================================
[6] After 100 million MCS run, the data is in CONFORMATIONS and then do data analysis
For example, convector-ok.f used for calculation of Ltrans (Fig.3)
timesg-ok.f and timeorder-ok.fis used for SG and Sb (Fig.4), rgxyz-ok.f is used for RE^2 and Rg^2 (Fig.6)
ENERGY_INTRA_INTER generated during MC run is used for Fig.2)
=============================================================================
Example data and analysis for B30 nanofiber with 100 million MCS at 300 K are given.

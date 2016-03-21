(Anthony) calculates the fake rates using
fakeratecalc.C

and

fakeratio.C

Other files needed:

runonefake.C
hadd_fake.sh
interface/fakeratecalc.h
interface/controlpannel_fakerate.h

These .h files don’t do much other than separate fakeratecalc.C from the main calculator.


fakeratecalc.C makes counts, pt, and eta distributions of electrons and muons in a QCD control region labeled CR1.
fakeratecalc.C is run using runonefake.C.
After it finishes running, run hadd_fake.sh, which just hadd’s the output root files into one file.
Then run fakeratio.C then takes counts/histograms of loose and tight and divides to get the fake rate.
The result is a file called fakeevent_result.root

Prompt rates:
Prompt rates are calculated using an event topology within my normal analysis code, tpdlmaincalc.C, that looks at the Z (with a b-jet veto,) and makes distributions of dilepton mass for OSee and OSmm when the leptons are tight-tight, tight-loose, and loose-loose.
The results are here:
eventloop_Data2lepAll.root
Where the relevant plots are those following the naming convention: h_Mll[SS,OS]DL0b[mm, ee][ ,1T, 0T]
For example:
        h_MllOSDL0bmm1T is a plot of Mll for OSDL muons where 1 is tight and the other is loose-but-not-tight.
        h_MllSSDL0bee   is a plot of Mll for SSDL electrons where both electrons are tight.


All these are processed by LJMet so if the problem is LJMet constricting the loose region, it won’t be apparent in the above code. I suspect this is the case because such a problem would explain both the disagreement in fake rate and why I see lower statistics than T 5/3.
For reference: Here are the raw counts of tight and loose-but-not-tight dileptons without any selection or filtration (except for the duplicate event filter an the CSC bad-event filter):

Counts of ssee: LL, TL, TT:     83      496     794
Counts of ssmm: LL, TL, TT:     160     457     4521
Counts of ssem: LL:     32      Loose ele+tight µ:      53      tight ele+loose µ:      484     TT:     1107

Counts of osee: LL, TL, TT:     183     3025    22710
Counts of osmm: LL, TL, TT:     2082    15657   429439
Counts of osem: LL:     40      Loose ele+tight µ:      233     tight ele+loose µ:      1094    TT:     9864# LepFakeRate

# 2HDMC-NLO
NLO improved 2HDMC. The original 2HDMC code can be found at https://2hdmc.hepforge.org/ . This project is based on 2HDMC 1.7.0.

## NLO Calculation
Currently, we only improve the on-shell decay calculation with NLO corrections. The NLO corrections are calculated using FeynArts/FormCalc using the 2HDM NLO model files in https://github.com/ycwu1030/THDMNLO_FA.

## NLO vs. LO
In the modification, we keep the original LO 2HDMC working. If you wish just using LO results in this code, please comment the line `CFLAGS+=-DNLOCOUPLINGS` in Makefile. Then `make distclean; make` to rebuild the program.

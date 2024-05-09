1. Download the asptranslate package from https://github.com/tomijanhunen.
2. Copy acyc2mip-FVS.c and acyc2mip-VE.c to asptranslate/translators
3. Modify asptranslate/translators/Makefile to include compiling acyc2mip-FVS.c and acyc2mip-VE.c
4. Compile according to the instruction given for the asptranslate package
5. Install gringo
6. USAGE:     gringo INPUT_FILE1 INPUT_FILE2 --output=smodels | ./lp2acyc -w | ./acyc2mip-VE --mip
         or   gringo INPUT_FILE1 INPUT_FILE2 --output=smodels | ./lp2acyc -w | ./acyc2mip-FSV --mip


#!/bin/sh
inf=1
while [ "$inf" = 1 ]
do
ans=$(zenity --title " 3D_DFT AUTOMATED PROGRAM MENU " --list --text "\n \n 
=======================================================\n
                WELCOME TO THE 3D-DFT AUTOMATED PROGRAM MENU \n
           HERE YOU CAN CHOOSE DIFFERENT PARAMETERS TO BE CHANGED  \n
        AUTOMATICALLY.THE OUTPUTS WILL BE ARCHIVED IN SEPERATE NAMES\n 
=======================================================\n\n\n" --radiolist  --column "choose" --column "                         TASKS" TRUE "AUTOMATIC MESH GENERATION UTILITY "  FALSE "AUTOMATIC ELECTRON NUMBER GENERATER"  --height=360 --width=400 )

    		
case "$ans"
in
	"AUTOMATIC MESH GENERATION UTILITY ")./ggen_meshinput.sh;;
	"AUTOMATIC ELECTRON NUMBER GENERATER")./ggen_electron.sh;;
	*)inf=2;;
esac
done

#!/bin/sh
inf=1
while [ "$inf" = 1 ]
do
dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title " 3D_DFT AUTOMATED PROGRAM MENU " --colors --infobox "\n \n \n\Zb\Z1
===============================================================================\n
                WELCOME TO THE 3D-DFT AUTOMATED PROGRAM MENU \n
           HERE YOU CAN CHOOSE DIFFERENT PARAMETERS TO BE CHANGED  \n
        AUTOMATICALLY.THE OUTPUTS WILL BE ARCHIVED IN SEPERATE NAMES\n 
===============================================================================\n\n\n
\Z0                      WHAT DO YOU WANT TO USE ???? \n
\n
    		=======================================================\n
    		|	(1)---AUTOMATIC MESH GENERATION UTILITY              |\n
    		|	(2)---AUTOMATIC ELECTRON NUMBER GENERATER            |\n 
    		 =======================================================\n \n
                        ENTER THE NO ::\Zb\Z1\Zu(will not be displayed)\ZU \n\n\n
                              \Zu TO EXIT TYPE  "q" \n \n \n \n \n " 27 88
read number
case "$number"
in
	1)./gen_meshinput.sh;;
	2)./gen_electron.sh;;
	q)inf=2;dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --colors --timeout 1 --msgbox " \n \Zb\Z1\Zu GOING BACK TO MAIN MENU " 10 40;;
esac
done

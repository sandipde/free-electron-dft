#!/bin/sh
inf=1
while [ "$inf" = 1 ]
do
dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "INPUT EDITOR FOR 3D-DFT PROGRAM" --colors --infobox "\n\n\n
\Zb\Z1 ==============================================================================\n
       THIS PROGRAM WILL LET YOU EDIT THE INPUTS OF THE 3DDFT PROGRAM \n
 ===============================================================================\n\n\n
                 \Z0   WHAT DO YOU WANT TO EDIT ???? \n
\n
    		=======================================================\n
    		|	(1)---MAIN PARAMETERS                                |\n
    		|	(2)---DAVIDSON'S PARAMETERS                          |\n
    		|	(3)---MESH PARAMETERS                                |\n 
    		|	(4)---EXTERNAL POTENTIAL                             |\n
    		 =======================================================\n \n
                        ENTER THE NO ::\Zb\Z1\Zu(will not be displayed)\ZU \n\n\n
                              \Zu TO EXIT TYPE  "q" \n \n \n \n \n 
         " 27 88
read number
case "$number"
in
	1) filename=main.in ;;
	2) filename=dvdson.in ;;
	3) filename=mesh.in;;
	4) filename=ext_pot.in ;;
	q) inf=2;;
esac
#printf "$filename \n"
case "$inf"
in
  1) #dialog --title "YOU ARE EDITING $filename FILE" --editbox $filename 27 88   ;;
	vi $filename;;
  2) dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --colors --timeout 1 --msgbox " \n \Zb\Z1\Zu GOING BACK TO MAIN MENU " 10 40
esac
done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# script by sandip de 19.02.09	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

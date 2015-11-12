#!/bin/sh
inf=1
while [ "$inf" = 1 ]
do
ans=0
ans=$(zenity  --list --title "INPUT EDITOR FOR 3D-DFT PROGRAM" --text "\n\n========================================================\n
       THIS PROGRAM WILL LET YOU EDIT THE INPUTS OF THE 3DDFT PROGRAM \n
 =======================================================\n\n" --radiolist --column "choose" --column "                         TASKS" TRUE "MAIN PARAMETERS " FALSE "DAVIDSON'S PARAMETERS"  FALSE "MESH PARAMETERS" FALSE "EXTERNAL POTENTIAL" --height=300 --width=400 )
#echo $ans    		 
case "$ans" 
in
	"MAIN PARAMETERS ") filename=main.in ;;
	"DAVIDSON'S PARAMETERS") filename=dvdson.in ;;
	"MESH PARAMETERS") filename=mesh.in;;
	"EXTERNAL POTENTIAL") filename=ext_pot.in ;;
	*) inf=2;;
esac
#printf "$filename \n"
case "$inf"
in
  1) #dialog --title "YOU ARE EDITING $filename FILE" --editbox $filename 27 88   ;;
	 zenity --text-info --title "YOU ARE EDITING $ans | ANY CHANGE YOU MADE WILL BE SAVED " --editable --filename=$filename --height=280 --width=640 >temp.in;mv temp.in $filename ;;
#  2) dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --colors --timeout 1 --msgbox " \n \Zb\Z1\Zu GOING BACK TO MAIN MENU " 10 40
esac
done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# script by sandip de 19.02.09	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

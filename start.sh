#!/bin/sh
printf "\n\n ================================================================\nCHECKING FOR ZENITY...."
(for i in `seq 1 40 100 `
do
echo $i
echo " # CHECKING FOR ZENITY.... $i % done"
sleep 0.5 ; done )|zenity --progress --text "CHECKING FOR ZENITY..... " --auto-close
case $?
in
0) printf ".....ok \n "
   printf "STARTING GRAPHICAL MENU >>>>>\n"
   ./gmenu.sh;;
*) printf ".....failed \n"
   printf "CHECKING FOR DIALOG ................" 
   dialog --infobox "DIALOG FOUND !!!" 5 30 
	case $?
	in
	0) dialog --colors --msgbox  "\Z0\Zb YOU APPEARED NOT TO HAVE \Z1zenity\Z0 \n    INSTALLED IN YOUR COMPUTER \n YOU WILL NOT BE ABLE TO USE THE\n         \Z1GRAPHICAL MENU " 10 40 
	   ./menu.sh ;;
	*) printf ".........failed \n" 
	   printf "COMPATIBILITY NOT SATISFIED >>>> EXITING !!!\n\n=================================================================\n";;
	esac;;
esac
printf "==============================EXITING !!! =================================\n"
printf "========================PROGRAM BY SANDIP DE===============================\n\n"



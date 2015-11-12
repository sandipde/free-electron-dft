#!/bin/sh
inf=1
while [ "$inf" = 1 ]
do
ans=$(zenity  --list --title " 3D_DFT VISUALIZER MENU " --text  "\n 
   ======================================================\n
                          WELCOME TO THE 3D-DFT PROGRAM \n
                   THIS IS THE VISUALIZER MENU OF THE PROGRAM \n
                       YOU CAN SEE ALL THE OUTPUTS HERE \n
   ======================================================\n\n\n" --radiolist --column "choose" --column "                         RESULTS" TRUE "SCREEN OUTPUTS" FALSE "RUN INFO" FALSE "ENERGY_ITERATION LOG " FALSE "EIGEN ENERGIES " FALSE "POISSON EQUATION SOLVER LOG" FALSE "DAVIDSON'S LOG" FALSE "FERMI-OCCUPANCY LOG " FALSE "RHO_CHECKING" FALSE "DISASTER MANAGEMENT LOG"  --height=500 --width=200 )
case "$ans"
in
	"SCREEN OUTPUTS") filename=screen.out;;
	"RUN INFO") filename=run.info;;
	"ENERGY_ITERATION LOG ") filename=energy.info;;
	"EIGEN ENERGIES ") filename=eigen.info;;
	"POISSON EQUATION SOLVER LOG") filename=poisson.info;;
	"DAVIDSON'S LOG") filename=davidson.info;;
	"FERMI-OCCUPANCY LOG ") filename=fermi_occupancy.info;;
	"RHO_CHECKING") filename=rho_check.info;;
	"DISASTER MANAGEMENT LOG") filename=disaster_management.info;;
	*) inf=2;;
esac
case "$inf" 
in
	1) zenity --text-info --title " YOU ARE WATCHING $ans" --filename=$filename --height=480 --width=640  ;;
#	2) dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --colors --timeout 1 --msgbox " \n \Zb\Z1\Zu GOING BACK TO MAIN MENU " 10 40;;
esac
done

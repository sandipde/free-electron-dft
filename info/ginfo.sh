#!/bin/sh
zenity --question --text "THIS PROGRAM CAN BE SEEN ONLY IN TERMINAL \n IF YOU ARE NOT RUNNING THE PROGRAM FROM TERMINAL \n press CANCEL press OK TO PROCEED" 
case $? 
in 
0)inf=1
while [ "$inf" = 1 ]
do
dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title " 3D_DFT VISUALIZER MENU " --colors --infobox "\n \n \n\Zb\Z1
===============================================================================\n
                WELCOME TO THE 3D-DFT PROGRAM \n
             THIS IS THE VISUALIZER MENU OF THE PROGRAM \n
             YOU CAN SEE ALL THE OUTPUTS HERE ....LIVE !!!\n
===============================================================================\n\n\n
\Z0            WHAT DO YOU WANT TO SEE ???? \n
\n
    		=======================================================\n
    		|	(1)---SCREEN OUTPUTS                                 |\n
    		|	(2)---RUN INFO                                       |\n
    		|	(3)---ENERGY_ITERATION LOG                           |\n 
    		|	(4)---STARTING EIGEN ENERGIES                        |\n
    		|	(5)---EIGEN ENERGIES                                 |\n
    		|	(6)---POISSON EQUATION SOLVER LOG                    |\n
    		|	(7)---DAVIDSON'S LOG                                 |\n 
    		|	(8)---FERMI-OCCUPANCY LOG                            |\n
    		|	(9)---RHO-CHECKING                                   |\n
    		|	(10)--DISASTER MANAGEMENT LOG                        |\n 
    		 =======================================================\n \n
                        ENTER THE NO ::\Zb\Z1\Zu(will not be displayed)\ZU \n\n\n
                              \Zu TO EXIT TYPE  "q" \n \n \n \n \n " 37 88
read number
case "$number"
in
	1) filename=screen.out;;
	2) filename=run.info;;
	3) filename=energy.info;;
	4) filename=1steigen.info;;
	5) filename=eigen.info;;
	6) filename=poisson.info;;
	7) filename=davidson.info;;
	8) filename=fermi_occupancy.info;;
	9)filename=rho_check.info;;
	10) filename=disaster_management.info;;
	q) inf=2;;
esac
case "$inf" 
in
	1) dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title " YOU ARE WATCHING $filename FILE " --tailbox $filename 45 140;;
	2) dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --colors --timeout 1 --msgbox " \n \Zb\Z1\Zu GOING BACK TO MAIN MENU " 10 40;;
esac
done;;
esac
clear

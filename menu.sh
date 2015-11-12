#!/bin/sh
inf=1
while [ "$inf" = 1 ]
do
dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title " 3D_DFT PROGRAM MAIN MENU " --colors --infobox "\n \n \n\Zb\Z1
===============================================================================\n
                         WELCOME TO THE 3D-DFT PROGRAM \n
                     THIS IS THE MAIN MENU OF THE PROGRAM \n
 ===============================================================================\n\n\n
\Z0            WHAT DO YOU WANT TO DO ???? \n
\n
    		=======================================================\n
    		|	(1)---\Z4EDIT INPUTS\Z0                                    |\n
    		|	(2)---\Z4COMPILE & RUN\Z0                                  |\n
    		|	(3)---\Z4VIEW THE OUTPUTS\Z0                               |\n 
    		|	(4)---\Z0ONLY COMPILE\Z0                                   |\n
    		|	(5)---\Z4SEPERATE WAVE FUNCTIONS FROM COMBINED RESULT\Z0   |\n
    		|	(6)---\Z2AUTOMATED RUN !!!\Z0                              |\n 
    		|	(7)---\Z1CLEAN UP\Z0(ALL RESULTS WILL BE DELEATED)         |\n
    		|	(8)---\Z4ARCHIVE THE RESULTS\Z0                            |\n 
    		 =======================================================\n \n
                        ENTER THE NO ::\Zb\Z1\Zu(will not be displayed)\ZU \n\n\n
                              \Zu TO EXIT TYPE  "q" \n \n \n \n \n " 27 88
read number
case "$number" 
in
	1) cd input/ ;./input.sh ;cd .. ;;
	2) make -C ./src/ 
		cd bin/
		nohup ./3d_dft.x>../info/screen.out &
		cd ..  
		dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --msgbox "THE PROGRAM IS STARTED :::::" 10 40 ;;
	3) cd info/ ; ./info.sh ;cd .. ;;
	4) make clean all -C ./src/ ;dialog --msgbox "THE PROGRAM COMPILED SUCESSFULY :::::" 5 50 ;;
	5) cd output_for_vmd/ ; ./eigen_vec_for_vmd.sh; cd .. ;;
	6) cd automate/ ; ./automate.sh ; cd .. ;;
	7) dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --infobox "ARE YOU SURE YOU WANT TO DELETE ALL ??(yes/n)" 5 50 
		read choice
		if [ $choice = yes ]
 		then		
		make clean -C ./src/;cd bin/;rm *.x ; cd ../info/;rm -f *.info *.out;cd ../output_for_vmd/
		rm -f *.dat *.VASP_CHGCAR ; cd ../back_up/ ; rm -f *.dat;cd ../output/ ; rm -f *.dat; cd .. 
		dialog --msgbox ":::::::::::::: CLEANED SUCCESSFULY:::::::::::::::::::" 5 50
		fi ;;
	8) dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --infobox "ENTER THE ARCHIVE NAME" 5 50
	   read tarname
	   day=`date |cut -c 9-10`
	   month=`date |cut -c 5-7`
	   year=`date |cut -c 25-28`
	   tar -cvvjf $tarname.$day$month$year.tar.bz2 info/ back_up/ output/ output_for_vmd/ chargeden_for_vmd/  input/
	   dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --msgbox ":::::::::THE ARCHIVE CREATED SUCESSFULY:::::::::::" 10 70 ;;
	q) inf=2;dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --colors --timeout 1 --msgbox " \n \Zb\Z1\Zu EXITING THE PROGRAM NOW " 10 40;;	
esac
done





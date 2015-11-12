#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#script by sandip de 19.02.09
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#!/bin/sh
zenity --question --text "THIS PROGRAM CAN BE SEEN ONLY IN TERMINAL \n IF YOU ARE NOT RUNNING THE PROGRAM FROM TERMINAL \n press CANCEL press OK TO PROCEED" 
case $? 
in 
0)inf=1
while [ "$inf" = 1 ]
do
dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "AUTOMATIC MESH INPUT GENERATER" --colors --infobox "\n \Zb\Z1
   ===================================================================================\n
	 This script will generate different configuration files for mesh, run the program \n
	               and store the results in different names \n
 ====================================================================================\n \Z0\n\n

	                which parameter you want to change ?\n\n
				\n
    		=======================================================\n
    		|	(1)---ALL 3 DIMENSIONS (ALL SAME)                    |\n
    		|	(2)---ALL 3 DELTAS  (ALL SAME)                       |\n
    		|	(3)---ONLY Z_DIMENSION                               |\n 
    		|	(4)---ONLY DELTA_Z                                   |\n
    		 =======================================================\n \n
                        ENTER YOUR CHOICE ::\Zb\Z1\Zu(will not be displayed)\ZU \n\n\n
                              \Zu TO EXIT TYPE  "q" \n \n \n \n \n " 25 90
read choice
day=`date |cut -c 9-10`
month=`date |cut -c 5-7`
year=`date |cut -c 25-28`
case "$choice"
in
   1)
	dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!" --infobox  "give the starting dimension ,end dimension , step dim :: " 10 50
	read startdim enddim stepdim
#	printf "give the end dimension :: "
#	read enddim
#	printf "give the increase in each step :: "
#	read stepdim
	dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!" --infobox  "give the fixed delta ::" 10 50
	read delta
 	sed s/0.1/$delta/g mesh.in >temp_mesh.in
	for i in `seq $startdim  $stepdim $enddim`
	do
	 sed  s/30/$i/g temp_mesh.in >../input/mesh.in
	dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!" --infobox  "input file genarated .....\n 
						executing the program now >>>>>>>>\n" 10 50
	 cd ../bin/
	 nohup ./3d_dft.x>../info/screen.out
	 cd ..
	  tar -cvvjf dim=$i.$day$month$year.tar.bz2 info/ output_for_vmd/ output/ input/
	cd automate/
	done
	dialog --msgbox "AUTOMATED WORK COMPLEATED" 10 30;;
   2)
  	dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!" --infobox "give starting delta ,end delta ,steps in delta ::" 10 50
        read startdelta enddelta stepdelta
#	printf "give end delta :: "
#	read enddelta
#	printf "give increase step ::"
#	read stepdelta
        dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!" --infobox "give the fixed dimension ::" 10 50
	read dim
	sed s/30/$dim/g mesh.in >temp_mesh.in
	for i in `seq $startdelta $stepdelta $enddelta`
	do
	 sed  s/0.1/$i/g temp_mesh.in >../input/mesh.in
   	 dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!" --infobox  "input file genarated .....
					\n executing the program now >>>>>>>>\n" 10 50
	 cd ../bin/
	 nohup ./3d_dft.x>../info/screen.out
	 cd ..
	 tar -cvvjf delta=$i.$day$month$year.tar.bz2 info/ output_for_vmd/ output/ input/
	cd automate/
	done
	dialog --msgbox "AUTOMATED WORK COMPLEATED" 10 30;;
   3)   dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!!" --msgbox "EDIT THIS FILE WITHOUT CHANGING THE Z_DIM VALUE" 10 50
	vi zdim_mesh.in
        dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!" --infobox "give starting z_dim ,end z_dim ,step in zdim ::" 10 50
        read start_zdim end_zdim stepzdim
#	dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!" --msgbox "give end z_dim ::" 10 50
#	read end_zdim
#	dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!" --msgbox "give increase step in zdim ::" 10 50
#	read stepzdim
	for i in `seq $start_zdim  $stepzdim $end_zdim`
	do
	 sed  s/300/$i/g zdim_mesh.in >../input/mesh.in
	printf "input file genarated .....\n executing the program now >>>>>>>>\n"
	  cd ../bin/
	 nohup ./3d_dft.x>../info/screen.out
	 cd ..
	  tar -cvvjf z_dim=$i.$day$month$year.tar.bz2 info/ output_for_vmd/ output/ input/
	cd automate/
	done
	dialog --msgbox "AUTOMATED WORK COMPLEATED" 10 30;;
   4)  dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!" --msgbox "EDIT THIS FILE WITHOUT CHANGING THE DELTA_Z VALUE" 10 50
	vi deltaz_mesh.in
        dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!" --msgbox "give starting delta_z, end delta_z, step in delta_z ::" 10 50
        read start_deltaz end_deltaz stepdeltaz
#	dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!" --msgbox "give end delta_z ::" 10 50
#	read end_deltaz
#	dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "IMPORTANT !!!" --msgbox "give increase step in delta_z ::" 10 50
#	read stepdeltaz
	for i in `seq $start_deltaz  $stepdeltaz $end_deltaz`
	do
	 sed  s/0.001/$i/g deltaz_mesh.in >../input/mesh.in
	printf "input file genarated .....\n executing the program now >>>>>>>>\n"
	  
	  cd ../bin/
	 nohup ./3d_dft.x >>../info/screen.out
	 cd ..
	  tar -cvvjf delta_z=$i.$day$month$year.tar.bz2 info/ output_for_vmd/ output/ input/
	cd automate/
	done
	dialog --msgbox "AUTOMATED WORK COMPLEATED" 10 30;;
  q) inf=2;dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --colors --timeout 1 --msgbox " \n \Zb\Z1\Zu GOING BACK TO AUTOMATED MENU " 10 40;;
esac
done;;
esac
clear


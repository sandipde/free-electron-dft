#!/bin/sh
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT BY SANDIP DE 21.02.09
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#!/bin/sh
inf=1
while [ "$inf" = 1 ]
do
ans=$(zenity  --list --title " 3D_DFT PROGRAM MAIN MENU " --text "\n ===================================================\n
                         WELCOME TO THE 3D-DFT PROGRAM \n
                     THIS IS THE MAIN MENU OF THE PROGRAM \n
 =====================================================\n\n\n" --radiolist --column "choose" --column "                         TASKS" TRUE "README" FALSE "EDIT INPUTS" FALSE "COMPILE & RUN" FALSE "VIEW THE OUTPUTS" FALSE "VIEW THE OUTPUTS LIVE !!!" FALSE "COMPILE ONLY" FALSE "SEPERATE MOLECULAR ORBITALS" FALSE "SEE THE 3D WAVE FUNCTION IN VMD" FALSE "LOAD MULTIPLE FILES IN VMD" FALSE "AUTOMATED RUN !!!" FALSE "CLEAN UP(ALL RESULTS WILL BE DELETED)" FALSE "ARCHIVE THE RESULTS" FALSE "KILL THE RUNNING PROGRAM"  --height=470 --width=400 )
#echo $ans    		 
case "$ans" 
in
	"README") zenity --title "README" --text-info --filename="readme" --height=500 --width=400 ;;
	"EDIT INPUTS") cd input/ ;./ginput.sh ;cd .. ;;
	"COMPILE & RUN")(echo "10"; echo "# COMPILING" ; 
			make -s -C ./src/ 
			case $?
			in
			   0)	cd bin/
				echo "50" ; echo "# EXECUTING THE PROGRAM....YOU CAN SEE THE OUTPUTS BY CHOOSING THE OPTION FROM MAIN MENU .."
				zenity --info --text ":::::THE PROGRAM IS STARTED :::::" --timeout=1
				nohup ./3d_dft.x>../info/screen.out ;echo "100"
				cd .. ;;
			   *) zenity --error --text "COMPILATION FAILED !!! ";;
			esac)|(zenity --progress --pulsate --title "PLEASE WAIT" --text "RUNNING  >>>" --auto-close --width=500 --height=50
				case $?
				in
				1) zenity --question --text "ARE YOU SURE ?? \n THIS WILL KILL THE PROGRAM "
					case $?
					in
						0)pkill "3d_dft.x";iter=2
						  zenity --info --text "YOU HAVE KILLED THE PROGRAM !! ";;
						   
					esac;;
				0)zenity --notification --text "PROGRAM COMPLEATED SUCESSFULLY ..." --timeout=3 ;iter=2
				  zenity --info --text "PROGRAM COMPLEATED SUCESSFULLY ..." ;;
				esac ) & ;;

	"VIEW THE OUTPUTS") cd info/ ; ./gfinalinfo.sh ;cd .. ;;
	"VIEW THE OUTPUTS LIVE !!!")cd info/ ; ./ginfo.sh ;cd .. ;;

	"COMPILE ONLY") (echo "10" echo "# cleaning ..." ;sleep 1 ;make  clean  -C ./src/ 
			echo "50" echo "# COMPILING ....";sleep 1
			make -s -C ./src/ 
			echo "98";echo "# COMPLETING >>>>")|zenity --progress --pulsate --text "COMPILING >>>" --auto-close
			case $?
			in
			   0) zenity --info --text "THE PROGRAM COMPILED SUCESSFULY :::::" --timeout=3 ;;
			   *) zenity --error --text "COMPILATION FAILED !!! ";;
			esac;;

	"SEPERATE MOLECULAR ORBITALS") cd output_for_vmd/ ; ./geigen_vec_for_vmd.sh; cd .. ;;
	"AUTOMATED RUN !!!") cd automate/ ;./gautomate.sh ; cd .. ;;
	"CLEAN UP(ALL RESULTS WILL BE DELETED)")zenity --question --text  "ARE YOU SURE YOU WANT TO DELETE ALL ??" 
		case $?
		in
		0)
		make clean -C ./src/;cd bin/;rm *.x ; cd ../info/;rm -f *.info *.out;cd ../output_for_vmd/
		rm -f *.dat *.VASP_CHGCAR ; cd ../back_up/ ; rm -f *.dat;cd ../output/ ; rm -f *.dat; cd .. 
		zenity --info --text ":::::::::::::: CLEANED SUCCESSFULY:::::::::::::::::::" --timeout=3 ;;
		1)zenity --info --text ":::::CLEANUP CANCELLED :::::" --timeout 2;;
		esac;;
	"ARCHIVE THE RESULTS") tarname=$(zenity --file-selection --title "GIVE ARCHIVE NAME" --save --confirm-overwrite)
	   name=`echo "$tarname"|wc -c`
#	   echo $tarname
#	   echo $name
	   case $name
	   in
		1) iter=1 ;;
		*)  
		   day=`date |cut -c 9-10`
		   month=`date |cut -c 5-7`
		   year=`date |cut -c 25-28`
		   (echo "# CREATING ARCHIVE::::"; echo "50" ; tar -cvzf "$tarname.$day$month$year.tar.gz" info/ back_up/ output/ output_for_vmd/ chargeden_for_vmd/ input/ ;echo "100" )|zenity --progress --pulsate --text "ARCHIVING ..." --auto-close
		   zenity --info --text ":::::::::THE ARCHIVE CREATED SUCESSFULY:::::::::::";; 
	   esac;;
	"KILL THE RUNNING PROGRAM")b=$(`a=`top -n 1 |grep "3d_dft.x"|awk ' { print $1 }'``)
					status=` echo $b | wc -c `
					case $status
					in
						1)zenity --error --text "YOU FOOL !!! \n THE PROGRAM IS NOT RUNNING !!";;
						*)zenity --question --text "ARE YOU SURE ?? \n THIS WILL KILL THE PROGRAM "
							case $?
							in	
								0)pkill "3d_dft.x";iter=2
						     		  zenity --info --text "YOU HAVE KILLED THE PROGRAM !! " --timeout 3;;
							esac;;
						   
					esac;;
	"SEE THE 3D WAVE FUNCTION IN VMD")filename=$(zenity --file-selection --title "SELECT THE FILE")
						case $?
						in
						0) vmd -CHGCAR "$filename" ;;
						esac ;;
	"LOAD MULTIPLE FILES IN VMD")vecmin=$(zenity --scale --text "pick your initial vector " --min-value=1 --max-value=30 --value=1 --step 2)
			vecmax=$(zenity --scale --text "pick your final vector " --min-value=1 --max-value=30 --value=1 --step 2)
			cd output_for_vmd
			printf "vmd -m" > vmd.sh
			for vec in `seq ${vecmin} 1 ${vecmax}`;do printf " -CHGCAR eigen_vec_spinup$vec.VASP_CHGCAR" ;done >>vmd.sh
			for vec in `seq ${vecmin} 1 ${vecmax}`;do printf " -CHGCAR eigen_vec_spindn$vec.VASP_CHGCAR" ;done >>vmd.sh
			chmod 777 vmd.sh
			./vmd.sh
			cd .. ;;
						  
	*)zenity --question --text "ARE YOU SURE ,YOU WANT TO EXIT ?? " 
		case $?
		in
			0) inf=2;;
		esac;;
esac
done




